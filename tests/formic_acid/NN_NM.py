#!/usr/bin/env python3


import argparse
import numpy as np
from ase import io
from ase.io.trajectory import Trajectory
from ase.io import read, write
from HessianNNCalculator.NNCalculator import *
from os.path import splitext
import ase.units as units

'''
Usage: python NN_NM.py -i structure.xyz
'''

#parse command line arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
required.add_argument("-i", "--input",   type=str,   help="input traj",  required=True)
optional = parser.add_argument_group("optional arguments")
optional.add_argument("--charge",  type=float, help="total charge", default=0.0)


args = parser.parse_args()
filename, extension = splitext(args.input)

atoms = read(args.input)
n_atom = len(atoms)

calc = NNCalculator(
    checkpoint="tl_models/tl_866_fad_ccsdt", #load the model you want to used
    atoms=atoms,
    charge=args.charge,
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


#attach calc to atoms object
atoms.set_calculator(calc)

#get hessian
hessian = np.reshape(calc.get_hessian(atoms), (3*n_atom, 3*n_atom))

m = atoms.get_masses()
indices = range(len(atoms))
indices = np.asarray(indices)

#calculate normal mode freq.
im = np.repeat(m[indices]**-0.5, 3)
omega2, modes = np.linalg.eigh(im[:, None] * hessian * im)
modes = modes.T.copy()

# Conversion factor:
s = units._hbar * 1e10 / np.sqrt(units._e * units._amu)
hnu = s * omega2.astype(complex)**0.5


print('---------------------\n')
print('  #    meV     cm^-1\n')
print('---------------------\n')
s = 0.01 * units._e / units._c / units._hplanck
for n, e in enumerate(hnu):
    if e.imag != 0:
        c = 'i'
        e = e.imag
    else:
        c = ' '
        e = e.real
    print('%3d %6.2f%s  %7.2f%s\n' % (n, 1000 * e, c, s * e, c))
print('---------------------\n')












