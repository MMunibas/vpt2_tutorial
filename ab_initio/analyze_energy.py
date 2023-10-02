#!/usr/bin/env python3


import argparse
import numpy as np
from ase import io
from ase.io.trajectory import Trajectory
from ase.io import read, write
from os.path import splitext
import ase.units as units

import sys

'''
Usage: python analyze_energy.py dataset.npz
'''


import matplotlib.pyplot as plt
fig, (ax1, ax2) = plt.subplots(2)
files = sys.argv
for f in files[1:]:
    data = np.load(f)
    E = data["E"]*23.06054
    Z = data["Z"]
    Z_unique = np.unique(Z, axis=0)

    ax1.hist(E,int(len(E)/50))

    ax1.set_xlabel('$E$ [Kcal/mol]')
    ax1.set_ylabel('Count')

    ax2.scatter(range(len(E)),E)
    
    
    
    
ax2.set_ylabel('$E$ [Kcal/mol]')
ax2.set_xlabel('Index')

plt.savefig("energy_analysis.pdf",bbox_inches='tight',dpi=250)
plt.show()

