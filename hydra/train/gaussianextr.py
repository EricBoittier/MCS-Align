#!/usr/bin/env python3

import argparse
from ase.units import *
import re
from ase.io import read, write
import os
import numpy as np


#parse command line arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
required.add_argument("-i", "--input",   type=str,   help="input outfile",  required=True)
args = parser.parse_args()
from os.path import splitext
filename, extension = splitext(args.input)

with open(args.input, "r") as f:
    with open(filename + ".xyz", "w") as g:
        contents = f.read().splitlines()

        ###search for natom
        linenumber0 = [i for i,line in enumerate(contents) if re.search(' Charge =  0 Multiplicity', line)][0]
        linenumber1 = [i for i,line in enumerate(contents) if re.search(' GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad', line)][0]
        natoms = linenumber1-linenumber0-4


        ###search for optimized geom
        linenumber2 = [i for i,line in enumerate(contents) if re.search(' Optimization completed.', line)][0]
        linenumber3 = np.array([i for i,line in enumerate(contents) if re.search('Standard orientation:', line)])

        #print(linenumber2)
        #print(linenumber3[np.where(linenumber3>linenumber2)])
        linenumberR = linenumber3[np.where(linenumber3>linenumber2)][0] + 5
        #print(contents[linenumberR])

        #search for frequencies
        freqtmp = []
        linenumberfreq = [i for i,line in enumerate(contents) if re.search(' Frequencies -- ', line)]
        for i in linenumberfreq:
            _, _, a, b, c, = contents[i].split()
            freqtmp.append(float(a))
            freqtmp.append(float(b))
            freqtmp.append(float(c))
        freqtmp = np.array(freqtmp)

        #search OH bonded symm OH frequency (somewhere between 3400 and 3700 or 3780 for H-pi)
        found = False
        freq = freqtmp[freqtmp<3800]
        freq = freq[freq>3400]
        if len(freq) > 1:
            print("ATTENTION: Frequency not clearly assigned")
        g.write(str(natoms) + '\n')
        for i in range(len(freq)):
            g.write(str(freq[i]) + '    ')
        g.write("\n")
        for line in contents[linenumberR:linenumberR+natoms]:
            _, l, _, x, y, z = line.split()

            if l == '1':
                l = 'H'
            elif l == '6':
                l = 'C'
            elif l == '7':
                l = 'N'
            elif l == '8':
                l = 'O'
            elif l == '9':
                l = 'F'
            else:
                print("UNKNOWN LABEL", l)
                quit()
            g.write(l + '    ' + str(x) + '    ' + str(y)+ '    ' + str(z) + '\n')



