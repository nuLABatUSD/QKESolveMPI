#include "../code/include.hh"
#include "mpi.h"


#include "../solvingProblems/WorkerHeaderCreation/headers/collision_distribution.hh"


#include <iostream>
#include <string>

using std::cerr;
using std::cout;
using std::string;

int main(int argc, char* argv[])
{
    int myid = 0;
    int numprocs = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    /*
        argv[1] = output filename
        argv[2] = collision mode supplied by the shell script

        mode 0 = original collisions
        mode 1 = optimized collisions
    */
    if (argc != 3)
    {
        if (myid == 0)
        {
            cerr
                << "usage: coll output_filename collision_mode\n"
                << "collision_mode:\n"
                << "    0 = original collisions\n"
                << "    1 = optimized collisions\n";
        }

        MPI_Finalize();
        return 1;
    }

    const string output_filename = argv[1];
    const int collision_mode = std::stoi(argv[2]);

    if (collision_mode != 0 &&
        collision_mode != 1)
    {
        if (myid == 0)
        {
            cerr
                << "Error: collision mode must be 0 or 1.\n";
        }

        MPI_Finalize();
        return 1;
    }

    const bool use_optimized_collisions =
        collision_mode == 1;

    int** core_job_rows = nullptr;

    if (use_optimized_collisions)
    {
        const int distribution_core_count =
            sizeof(core_jobs) / sizeof(core_jobs[0]);

        if (distribution_core_count != numprocs)
        {
            if (myid == 0)
            {
                cerr
                    << "Error: collision distribution contains "
                    << distribution_core_count
                    << " core rows, but MPI started "
                    << numprocs
                    << " ranks.\n";
            }

            MPI_Finalize();
            return 1;
        }

        core_job_rows =
            new int*[numprocs];

        for (int rank = 0;
             rank < numprocs;
             ++rank)
        {
            core_job_rows[rank] =
                core_jobs[rank];
        }
    }

    linspace_and_gl* eps =
        new linspace_and_gl(
            0.,
            EPS_MAX_LINSPACE,
            EPS_LINSPACE_POINTS,
            5
        );

    density* ics =
        new density(
            eps,
            IC_TCM,
            IC_NU_E,
            IC_NU_MU,
            IC_NUBAR_E,
            IC_NUBAR_MU,
            IC_MAX_DISTFUN
        );

    ics->set_T_Tcm(
        IC_TEMP,
        IC_TCM
    );

    QKEMPI* sim =
        new QKEMPI(
            myid,
            numprocs,
            PARAM_SIN_2THETA,
            PARAM_DELTA_M_SQUARED,
            0.,
            PARAM_DT_INIT,
            eps,
            ics,
            use_optimized_collisions,
            core_job_rows,
            max_jobs
        );

    /*
        The optimized collisions constructor copies the job indexes into
        its own worker_result_indexes arrays during construction.
    */
    delete[] core_job_rows;

    if (myid == 0)
    {
        cout
            << "Running QKEMPI with "
            << (
                use_optimized_collisions
                    ? "optimized"
                    : "original"
            )
            << " collisions.\n";
    }

    sim->run(
        PARAM_N_STEPS,
        PARAM_DN,
        PARAM_T_FINAL,
        output_filename,
        true
    );

    delete sim;
    delete ics;
    delete eps;

    MPI_Finalize();

    return 0;
}
