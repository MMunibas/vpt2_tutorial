#!/usr/bin/env python

import sys
import numpy as np
from goptimizer import parse_ifile, write_ofile
from ase import Atoms
from HessianNNCalculator.NNCalculator import *
import numpy as np
import time
#sys.stderr = open('errorlog.txt', 'w')


if __name__ == "__main__":

    # The input and output filenames from gaussian (don't change).
    ifile = sys.argv[2]
    ofile = sys.argv[3]

    # What you get out from Gaussian.
    # Cartesian coordinates in angstrom, use 0.52917721092 to get Bohr.
    (natoms, deriv, charge, spin, atomtypes, coordinates) = parse_ifile(ifile)

#===============================================================================
#NN part, get nn energy and force and whatever needed.

    atoms = Atoms(atomtypes, coordinates)
    
    calc = NNCalculator(
        checkpoint="/home/kaeser/teaching/vpt2_tutorial/tutorial_vpt2_sena/vpt2_hcooh/tl_models/tl_866_fad_ccsdt", #load the model you want to used
        atoms=atoms,
        charge=0.0,
        F=128,
        K=64,
        num_blocks=5,
        num_residual_atomic=2,
        num_residual_interaction=3,
        num_residual_output=1,
        sr_cut=10.0,
        use_electrostatic=True,
        use_dispersion=True,
        s6=1.0000,                    #s6 coefficient for d3 dispersion, by default is learned
        s8=2.3550,                    #s8 coefficient for d3 dispersion, by default is learned
        a1=0.5238,                    #a1 coefficient for d3 dispersion, by default is learned
        a2=3.5016)                   #a2 coefficient for d3 dispersion, by default is learned)





    atoms.set_calculator(calc)

    #get energy prediction from NN
    energy = atoms.get_potential_energy() #output in eV
    energy /= 27.211386024367243 #conversion to hartree
    
    gradient = -atoms.get_forces() #output in eV/angstrom
    gradient = (gradient*0.5291772105638411)/27.211386024367243 #converstion to hartree/bohr

    print(gradient)
    hessian = np.reshape(calc.get_hessian(atoms), (3*natoms, 3*natoms)) #output in ev/angstr**2

    #gaussian needs only lower triangle of hessian, thus reformat
    hessiantmp = []
    
    count=1
    for i in range(3*natoms):
        for j in range(count):
            hessiantmp.append(hessian[i,j])
        count+=1


    hessian = np.reshape(np.array(hessiantmp), (-1,3))
        
    hessian = (hessian*0.5291772105638411**2)/27.211386024367243 #conversion to ha/bohr**2

    #mu = np.dot(atoms.get_charges(), atoms.get_positions())/0.20819433442462576 #conversion to debye 
    mu = np.dot(atoms.get_charges(), atoms.get_positions())/0.5291772105638411 #conversion to e*bohr 

    dipder = np.array(calc.get_dipder(atoms))#being in e
    dipdertmp = []
    for i in range(3*len(atoms)):
        dipdertmp.append(dipder[0].reshape(-1)[i])
        dipdertmp.append(dipder[1].reshape(-1)[i])
        dipdertmp.append(dipder[2].reshape(-1)[i])
    dipder = np.array(dipdertmp).reshape(-1, 3)
    
#===============================================================================
    # Produce the Gaussian input file with the supplied data
    write_ofile(ofile, energy, natoms, gradient=gradient, hessian=hessian, dipole=mu, dipole_derivative=dipder)
    #write_ofile(ofile, energy, natoms, gradient=gradient, hessian=hessian, dipole=mu)

#sys.stderr.close()
#sys.stderr = sys.__stderr__



