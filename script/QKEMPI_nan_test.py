import numpy as np
import sys

data = np.loadtxt(sys.argv[1]+".csv", delimiter=',')

d = True

if np.sum(np.isnan(data[-1,:]))==0:
    rho = data[-1,2:]
elif np.sum(np.isnan(data[-2,:]))==0:
    rho = data[-2,2:]
else:
    print("Error, final two lines of data file are nan")
    d = False

if(d):
    f = open("density_test.hh", 'w')
    f.write("double test_density[] = {\n")
    for i in range(len(rho)):
        f.write("    ")
        f.write(str(rho[i]))
        if i != len(rho) - 1:
            f.write(",\n")
    f.write("};")
    f.close()