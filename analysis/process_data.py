import numpy as np
import matplotlib.pyplot as plt
import subprocess

hbar = 6.582e-22
GF = 1.166e-11

NU_E = 0
NU_MU = 1
ANTI_NU_E = 2
ANTI_NU_MU = 3

def make_data_dictionary(data_file, eps_file):
    results = dict()
    
    data = np.loadtxt(data_file, delimiter=',')
    eps_data = np.loadtxt(eps_file, delimiter=',')
    num_bins = len(eps_data[0])

    results['time'] = data[:,0] * hbar * 1e6
    results['N_bins'] = num_bins

    rho = data[:,2:2+4*num_bins]
    rhobar = data[:,2+num_bins*4:-2]
    results['rho'] = rho
    results['rhobar'] = rhobar

    P0 = rho[:,::4]
    Pz = rho[:,3::4]
    P0bar = rhobar[:,::4]
    Pzbar = rhobar[:,3::4]

    f_e = 0.5*P0*(1+Pz)
    f_m = 0.5*P0*(1-Pz)
    f_ebar = 0.5*P0bar*(1+Pzbar)
    f_mbar = 0.5*P0bar*(1-Pzbar)

    results['f'] = [f_e, f_m, f_ebar, f_mbar]
    
    results['Tcm'] = data[-1,-1]

    results['eps'] = eps_data[0]
    results['w'] = eps_data[1]

    dnde = []
    for i in range(4):
        dnde.append(results['f'][i] * results['eps']**2 / (2 * np.pi**2))
    results['dnde'] = dnde

    return results

def make_P(rho):
    return rho[:,::4], rho[:,1::4], rho[:,2::4], rho[:,3::4]

def V_mat(data):
    eps = data['eps']
    w = data['w']

    P0, Px, Py, Pz = make_P(data['rho'])
    P0bar, Pxbar, Pybar, Pzbar = make_P(data['rhobar'])

    Vx = 0.5 * P0 * Px - 0.5 * P0bar * Pxbar
    Vy = 0.5 * P0 * Py - 0.5 * P0bar * Pybar
    Vz = 0.5 * P0 * Pz - 0.5 * P0bar * Pzbar

    z = np.zeros(len(Vz))
    y = np.zeros(len(Vz))
    x = np.zeros(len(Vz))

    res = np.zeros((3, len(Vz)))
    for i in range(len(z)):
        zz = Vz[i,:] * eps**2
        yy = Vy[i,:] * eps**2
        xx = Vx[i,:] * eps**2
        res[2,i] = np.sum(w * zz)
        res[1,i] = np.sum(w * yy)
        res[0,i] = np.sum(w * xx)

    return np.sqrt(2) * GF * data['Tcm']**3 * res
    
def s_classical(data, time_index):
    eps = data['eps']
    w = data['w']
    se = np.zeros(len(eps))
    for i in range(len(se)):
        for j in range(4):
            fe = data['f'][j][time_index,:]
            if fe[i] > 0 and fe[i] < 1:
                se[i] += fe[i] * np.log(fe[i]) + (1-fe[i]) * np.log(1-fe[i])
    se *= -eps**2 / (2 * np.pi**2)

    return np.sum(se * w) * data['Tcm']**3

def s_quantum(data, time_index):
    eps = data['eps']
    w = data['w']

    density = [data['rho'][time_index,:], data['rhobar'][time_index,:]]
    se = np.zeros(len(eps))
    for i in range(len(se)):
        for r in density:
            den = r[4*i:4*i+4]
            lam1 = 0.5 * den[0] + 0.5 * den[0] * np.sqrt(np.sum(den[1:]**2))
            lam2 = 0.5 * den[0] - 0.5 * den[0] * np.sqrt(np.sum(den[1:]**2))
    
            if lam1 > 0 and lam1 < 1:
                se[i] += lam1 * np.log(lam1) + (1-lam1) * np.log(1-lam1)
            if lam2 > 0 and lam2 < 1:
                se[i] += lam2 * np.log(lam2) + (1-lam2) * np.log(1-lam2)

    se *= -eps**2 / (2 * np.pi**2)
    return np.sum(se * w) * data['Tcm']**3

def entropy(data):
    s = np.zeros_like(data['time'])
    sq = np.zeros_like(s)

    for i in range(len(s)):
        s[i] = s_classical(data, i)
        sq[i] = s_quantum(data, i)

    return s, sq

def num_density(data):
    eps = data['eps']
    w = data['w']

    dn = np.zeros((4,len(data['time']),data['N_bins']))
    for i in range(4):
        for j in range(len(data['time'])):
            dn[i,j,:] = data['f'][i][j,:] * eps**2 / (2 * np.pi**2)

    n = np.zeros((4, len(data['time'])))
    for i in range(4):
        for j in range(len(data['time'])):
            n[i,j] = np.sum(dn[i,j,:] * w)

    return n * data['Tcm']**3

def run_coherentsolve(outfile_name, Tcm=32, dm2=1.e-18, ke=0.9, km=1.8, kebar=0.9, kmbar=1.8):
    res = subprocess.run("cd .. && bash script/run_coherent.sh {} {} {} {} {} {} analysis/{}".format(Tcm, ke, km, kebar, kmbar, dm2, outfile_name), shell=True, capture_output=True)

def find_coherent_frequency_MHz(Tcm=32, dm2=1.e-18, ke=0.9, km=1.8, kebar=0.9, kmbar=1.8):
    run_coherentsolve("asym", Tcm, dm2, ke, km, kebar, kmbar)

    asym = make_data_dictionary("asym_run.csv", "asym_eps.csv")

    V_asym = V_mat(asym)
    t = asym['time']
    dv = np.diff(V_asym[1])
    t_max = []
    for i in range(len(dv)-1):
        if dv[i] >= 0 and dv[i+1] < 0:
            t_zero = t[i] - dv[i] * (t[i+1]-t[i])/(dv[i+1]-dv[i])
            t_max.append(t_zero)
    
    f = (1/np.diff(t_max))[1:]

    if np.std(f) < 0.01 * np.mean(f):
        subprocess.run("rm asym_run.csv", shell=True, capture_output=True)
        subprocess.run("rm asym_eps.csv", shell=True, capture_output=True)
        return np.mean(f)
    else:
        plt.plot(t, V_asym[1])
        plt.xlabel("t")
        plt.ylabel(r"$V_y$")
        print("Error: significant variance in frequency calculation.")
        return np.mean(f)
        
        