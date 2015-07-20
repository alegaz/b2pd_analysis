#!/usr/bin/env python
# -*- coding: utf-8 -*-

from basf2 import *
from modularAnalysis import inputMdst
from modularAnalysis import inputMdstList
from modularAnalysis import analysis_main
from simulation import add_simulation
from reconstruction import add_reconstruction
from reconstruction import add_mdst_output

import sys

# load input ROOT file

inputName = ['testSkimming_2-7.root','testSkimming_8-10.root']
outputName = 'reconstructedContinuum.root'

inputMdstList(inputName)

bkgdir = '/b2data/gaz/local_installation/bkg/bkg/'
bkgFiles = [bkgdir + '*.root']

# simulation
add_simulation(analysis_main)

# reconstruction
add_reconstruction(analysis_main)

# dump in MDST format
add_mdst_output(analysis_main, True,outputName)

# Show progress of processing
progress = register_module('Progress')
analysis_main.add_module(progress)

# Process the events
process(analysis_main)

# print out the summary
print statistics
