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
from modularAnalysis import loadGearbox
from modularAnalysis import analysis_main
from reconstruction import add_mdst_output
from ROOT import Belle2

# generate signal MC
generateY4S(1000000, Belle2.FileSystem.findFile('../dec_files/genericBBbar.dec'))

# dump in MDST format
add_mdst_output(analysis_main, True,
                'Gen_bbar_gsim-BKGx1-1M.root')

# if simulation/reconstruction scripts are not added than one needs to load gearbox
loadGearbox()

# Process the events
process(analysis_main)

# print out the summary
print statistics
