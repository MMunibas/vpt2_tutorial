#!/usr/bin/env python3


import argparse
import numpy as np
from ase.io import read, write
from NNCalculator.NNCalculator import *
import argparse
from os.path import splitext
from ase.vibrations import Vibrations

'''
Usage: python ase_vibrations.py -i structure.xyz
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

#setup calculator object, which in this case is the NN calculator
#it is important that it is setup with the settings as used in the
#training procedure.
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

#setup calculator (which will be used to describe the atomic interactions)
atoms.set_calculator(calc)

vib = Vibrations(atoms)
vib.run()

#print frequencies
vib.summary()

#delete written files
vib.clean()









