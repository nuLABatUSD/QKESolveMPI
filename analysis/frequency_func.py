import numpy as np
import matplotlib.pyplot as plt
import process_data as pd

def return_data(outfile_name):
    run = outfile_name+"_run.csv"
    epsilon = outfile_name+"_eps.csv"
    data = pd.make_data_dictionary(run, epsilon)
    return data

def calc_Vs(data,dm2):
    P0, Px, Py, Pz = pd.make_P(data['rho'])
    P0bar, Pxbar, Pybar, Pzbar = pd.make_P(data['rhobar'])
    eps = data['eps']
    w=data['w']
    time=data['time']
    Tcm=data['Tcm']
    
    S=np.array([(P0*Px)+(P0bar*Pxbar),(P0*Py)+(P0bar*Pybar),(P0*Pz)+(P0bar*Pzbar)])
    A=np.array([(P0*Px)-(P0bar*Pxbar),(P0*Py)-(P0bar*Pybar),(P0*Pz)-(P0bar*Pzbar)])
    B=np.array([0.8,0,-0.6]) #sin2theta=0.8
    
    Ntrap=201
    Ngl=5
    Emax=20
    G=np.sqrt(2) * pd.GF * (Tcm**3)*0.5/(np.pi**2)
    V=pd.V_mat(data)
    
    timeInd=len(P0)
    Vdd1=np.empty((3,timeInd),dtype=float)
    Vdd2=np.empty((3,timeInd),dtype=float)
    Vdd3=np.empty((3,timeInd),dtype=float)
    Vd=np.empty((3,timeInd),dtype=float)
    for t in range(timeInd):
      v1=0
      v2=0
      v3=0
      v=0
      for e in range(len(eps)):
        v1+= w[e]*np.dot(B,A[:,t,e])*B
        v2+= -w[e]*A[:,t,e]
        v3+= w[e]*(eps[e])*np.dot(B,S[:,t,e])
        v+= w[e]*(eps[e])*np.cross(B,S[:,t,e])
      Vdd1[:,t]=v1*G*(dm2**2)*0.25/(Tcm**2) # substituted p=E=eps*Tcm
      Vdd2[:,t]=v2*G*(dm2**2)*0.25/(Tcm**2)
      Vdd3[:,t]=v3*V[:,t]*G*dm2*0.5/Tcm
      Vd[:,t]=v*G*dm2*0.5/Tcm
    
    return V,Vd,Vdd1,Vdd2,Vdd3
    
def frequency(data,dm2):
    Vs=calc_Vs(data,dm2)
    V=Vs[0]
    Vd=Vs[1]
    Vdd=Vs[2]+Vs[3]+Vs[4]
    timeInd=len(V[0])
    time=data['time']
    
    co=np.empty((3,timeInd),dtype=float)
    for i in range(timeInd):
        for j in range(len(V)):
            if V[j,i]==0:
                co[j,i]=0
            else:
                co[j,i]=Vdd[j,i]/V[j,i]
                
    Dco = np.diff(co[1])
    coSign = np.sign(Dco)
    switch=np.diff(coSign)
    co0=np.array([])
    for i in range(len(coSign)-1):
      if switch[i]!=0:
          co0=np.append(co0,i)
    co0=co0.astype(int)
    cycle=co0[::4] #2 sign changes per cycle
    
    period = (time[cycle[-1]]-time[cycle[0]])/(len(cycle)-1)/pd.hbar*1e-6 #1e-6 for microsecond to second conversion
    frequency = 1/period #MeV
    af = frequency*2*np.pi
    min=np.min(co[1])
    
    return co,co0,min,af


def full(outfile_name, Tcm=32, dm2=1e-18, graph=True, ke=0.9, km=1.8, kebar=0.9, kmbar=1.8):
    # graph='yes' or 'no'
    pd.run_coherentsolve(outfile_name, Tcm, dm2, ke, km, kebar, kmbar)
    data=return_data(outfile_name)
    freq=frequency(data,dm2)
    
    time=data['time']
    co=freq[0]
    cycle=freq[1]
    min=freq[2]
    af=freq[3]
    
    if graph==True:
        c=cycle[2]
        plt.plot(time[:c],[0]*c, '--') #x axis is blue
        plt.plot(time[:c],[-(af**2)]*c, label='angular frequency squared (-ω²)') #angular frequency squared
        plt.plot(time[:c],co[1,:c], label='V**/V') #integral
        plt.plot(time[:c],[min]*c,'--')
        plt.title('Components of V** vs. Time(μs)')
        plt.legend()
        
    return min,af