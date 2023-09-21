#!/usr/bin/env python
#
# MIT License
#
# Copyright (c) 2018 Anders Steen Christensen
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import sys
import numpy as np
from goptimizer import parse_ifile, write_ofile
from ase import Atoms
from NNCalculator.NNCalculator import *
import numpy as np

if __name__ == "__main__":
    # The input and output filenames from gaussian (don't change).
    ifile = sys.argv[2]
    ofile = sys.argv[3]
    # What you get out from Gaussian.
    # Cartesian coordinates in angstrom, use 0.52917721092 to get Bohr.
    (natoms, deriv, charge, spin, atomtypes, coordinates) = parse_ifile(ifile)

#===============================================================================
#NN part, get nn energy and force
    #quit()
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
    energy = atoms.get_potential_energy() #output in eV
    energy /= 27.211386246 #conversion to hartree
    gradient = -atoms.get_forces() #output in eV/angstrom
    gradient = (gradient*0.52917721092)/27.211386245
    mu = np.dot(atoms.get_charges(), atoms.get_positions())/0.5291772105638411 #conversion to e*bohr 
#===============================================================================

    # Produce the Gaussian input file with the supplied data
    write_ofile(ofile, energy, natoms, gradient=gradient, dipole=mu)
