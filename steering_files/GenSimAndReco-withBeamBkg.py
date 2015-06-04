#!/usr/bin/env python
# -*- coding: utf-8 -*-

#######################################################                                                                                     #
#
# Adapted from B2A104 tutorial
#
# Author: A.Gaz
#
######################################################

# gbasf2 setting
gb2_setuprel = 'build-2015-04-19'
gb2_input_sandboxfiles = '../dec_files/B0_phi-K+K-_KS-pi+pi-.dec'

from basf2 import *
from modularAnalysis import generateY4S
from modularAnalysis import analysis_main
from simulation import add_simulation
from reconstruction import add_reconstruction
from reconstruction import add_mdst_output
from ROOT import Belle2
import glob

# specify number of events to be generated (this is mostly for the grid)
eventinfosetter = register_module('EventInfoSetter')
eventinfosetter.param('evtNumList', [10])  
eventinfosetter.param('runList', [1])  # from run number 1
eventinfosetter.param('expList', [1])  # and experiment number 1
#analysis_main.add_module(eventinfosetter)

# generate signal MC
generateY4S(10, Belle2.FileSystem.findFile('analysis/ale_stuff/dec_files/B0_phi-K+K-_KS-pi+pi-.dec'))
#generateY4S(10, glob.glob('B0_phi-K+K-_KL.dec'))

# simulation
#add_simulation(analysis_main, None, bkgFiles)
add_simulation(analysis_main, None, bkgfiles=glob.glob('/sw/belle2/bkg/*.root'))

# reconstruction
add_reconstruction(analysis_main)

# dump in MDST format
add_mdst_output(analysis_main, True,
                'B0_Phi,KK_Ks-pi+pi-_gsim-BKGx1-100.root')

# Process the events
process(analysis_main)

# print out the summary
print statistics
