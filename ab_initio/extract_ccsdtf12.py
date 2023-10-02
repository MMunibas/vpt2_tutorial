#!/usr/bin/env python3


from ase.units import *
import re
from ase.io import read, write
import os
import numpy as np


file_list = []
for path, dirs, files in os.walk("ch2o_ccf12_inp_molpro"):
    for file in files:
        if file.endswith(".out"):
            file_list.append(path+"/"+file)

file_list = np.sort(file_list)
#print(file_list)
#quit()
Nmax = 6

num = len(file_list)
N = np.zeros([num], dtype=int)
E = np.zeros([num], dtype=float)
Q = np.zeros([num], dtype=float)
D = np.zeros([num, 3], dtype=float)
Z = np.zeros([num, Nmax], dtype=int)
R = np.zeros([num, Nmax, 3], dtype=float)
F = np.zeros([num, Nmax, 3], dtype=float)

#atom energies from uccsd(t)-f12b,ansatz=3C(fix) with avtz-f12 basis (molpro)
#in hartree
Eref = np.zeros([10], dtype=float)
Eref[1] = -0.499946213283
Eref[6] = -37.788204984713
Eref[8] = -75.000839553994


index = 0
for file in file_list:
    # open file and read contents
    with open(file, "r") as f:
        contents = f.read().splitlines()

    #search for CARTESIAN COORDINATES:
    linenumber = [i for i,line in enumerate(contents) if re.search(' geometry={', line)][0]
    linenumber1 = [i for i,line in enumerate(contents) if re.search(' }', line)][0]
    Ntmp = linenumber1-linenumber-3
    #print(Ntmp)
    Ztmp = []
    Rtmp = []
    for line in contents[linenumber+3:linenumber+3+Ntmp]:
        l, x, y, z = line.split()
        if l == 'H':
            Ztmp.append(1)
        elif l == 'C':
            Ztmp.append(6)
        elif l == 'N':
            Ztmp.append(7)
        elif l == 'O':
            Ztmp.append(8)
        else:
            print("UNKNOWN LABEL", l)
            quit()
        Rtmp.append([float(x), float(y), float(z)])
    #print(Ztmp)
    #print(Rtmp)
    #quit()
    #search for forces (CARTESIAN GRADIENT):
    linenumberF = [i for i,line in enumerate(contents) if re.search('MP2 GRADIENT FOR STATE 1.1', line)][0]# MP2 gradient is output from molpro even if its ccsdt (its a bug, asked Max)
    #print(linenumberF)
    Ftmp = []
    for line in contents[linenumberF+4:linenumberF+4+Ntmp]:
        _,x, y, z = line.split()
        Ftmp.append([-float(x), -float(y), -float(z)])
    #print(Ftmp)
    
    #search for DIPOLE MOMENT:
    linenumberD = [i for i,line in enumerate(contents) if re.search('Dipole moment /Debye', line)][1]
    #print(linenumberD)
    Dx = float(contents[linenumberD].split()[3])
    Dy = float(contents[linenumberD].split()[4])
    Dz = float(contents[linenumberD].split()[5])
    Dtmp = [Dx, Dy, Dz]
    #print(Dtmp)
    #quit()


    #search for FINAL SINGLE POINT ENERGY:
    linenumberE = [i for i,line in enumerate(contents) if re.search('aug-cc-pVTZ-F12 energy=', line)][0]

    Etmp = float(contents[linenumberE].split()[2])
    #print(Etmp)
    #quit()
    #subtract asymptotics
    for z in Ztmp:
        Etmp -= Eref[z]

    #search for TOTAL CHARGE
    #linenumberQ = [i for i,line in enumerate(contents) if re.search('symmetry,nosym', line)][0]

    Qtmp = 0.0
    #print(Qtmp)

    N[index] = Ntmp
    E[index] = Etmp
    Q[index] = Qtmp
    D[index,:] = np.asarray(Dtmp)
    Z[index,:Ntmp] = np.asarray(Ztmp)
    R[index,:Ntmp,:] = np.asarray(Rtmp)
    F[index,:Ntmp,:] = np.asarray(Ftmp)

    index += 1
    if index%100 == 0:
        print(index/num*100,"%")


#unit conversion
E *= Hartree  # *27.211386024367243 conversion to eV
D *= Debye # *0.20819433442462576 conversion to e*angstrom
F *= Hartree/Bohr # *27.211386024367243/0.5291772105638411 conversion to eV/angstrom

np.savez_compressed("dataset.npz", N=N, E=E, Q=Q, D=D, Z=Z, R=R, F=F)


