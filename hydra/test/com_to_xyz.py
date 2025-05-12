#!/usr/bin/env python3
'''
file to convert gview data into xyz datausage with bash command: 
for filename in ./*.com; do ./com_to_xyz.py -i "$filename"; done
CAUTION WITH USING IT, ADAPT LINENUMBERS !!!
'''

from ase.units import *
import re
from ase.io import read, write
import os
import numpy as np
from os.path import splitext
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
required.add_argument("-i", "--input",   type=str,   help="input traj",  required=True)
args = parser.parse_args()


filenametot = args.input
filename, extension = splitext(filenametot)

with open(filenametot, "r") as f:
    contents = f.read().splitlines()

linenumber = [i for i,line in enumerate(contents) if re.search('0 ', line)][0]
linenumber1 = [i for i,line in enumerate(contents) if re.search(' 1 ', line)][0]
N = linenumber1-linenumber-2
with open(filename + ".xyz", "w") as g:
    g.write(str(N) +"\n\n")
    for line in contents[linenumber+1:linenumber+1+N]:
        g.write(line + "\n")

    


print(linenumber1)
print(N)

