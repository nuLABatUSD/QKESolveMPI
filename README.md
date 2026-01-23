**Requirements**
- clone ``base_code``. Both ``QKESolveMPI`` and ``base_code`` should be at the same level
- file ``run_params.hh`` in the folder

**Scripts**
- ``bash script/run_coherent.sh`` Tcm k_e k_mu k_ebar k_mubar dm2 file_name : runs coherent oscillations
- ``bash script/test_QKEMPI_Coll.sh`` N_cores : calculates R with N_cores
- ``bash script/run_collisionsQKE_MPI.sh`` N_cores : calculates C with N_cores
- ``bash script/run_QKEMPI.sh`` file_name : solves QKEMPI (as defined in ``run_params.hh``) with output file file_name
- ``bash script/compile_QKEMPI.sh`` file_name : compiles QKEMPI (as defined in ``run_params.hh``) to prepare for ``sbatch`` on Anvil
