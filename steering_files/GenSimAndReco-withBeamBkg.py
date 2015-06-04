#!/usr/bin/env python
# -*- coding: utf-8 -*-

#######################################################                                                                                     #
#
# Adapted from B2A104 tutorial
#
# Author: A.Gaz
#
######################################################

from basf2 import *
from modularAnalysis import generateY4S
from modularAnalysis import analysis_main
from simulation import add_simulation
from reconstruction import add_reconstruction
from reconstruction import add_mdst_output
from ROOT import Belle2
import glob

# generate signal MC
generateY4S(100, Belle2.FileSystem.findFile('../dec_files/B0_phi-K+K-_KS-pi+pi-.dec'))

# simulation
add_simulation(analysis_main, None, bkgfiles=glob.glob('/sw/belle2/bkg/*.root'))

# reconstruction
add_reconstruction(analysis_main)

# dump in MDST format
add_mdst_output(analysis_main, True,
                'B0_Phi-KK_Ks-pi+pi-_gsim-BKGx1.root')

# Process the events
process(analysis_main)

# print out the summary
print statistics
