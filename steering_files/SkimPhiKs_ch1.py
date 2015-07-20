#!/usr/bin/env python
# -*- coding: utf-8 -*-

# load all necessary scripts
from basf2 import *
from reconstruction import *

from modularAnalysis import inputMdstList
from modularAnalysis import fillParticleListsFromMC
from modularAnalysis import applyCuts
from modularAnalysis import summaryOfLists
from modularAnalysis import reconstructDecay
from modularAnalysis import matchMCTruth
from modularAnalysis import analysis_main
from modularAnalysis import ntupleFile
from modularAnalysis import ntupleTree
from modularAnalysis import copyLists
from modularAnalysis import cutAndCopyList

import sys

inputName = ['Gen_bbar_gsim-BKGx1-1M-1.root',
             'Gen_bbar_gsim-BKGx1-1M-2.root',
             'Gen_bbar_gsim-BKGx1-1M-3.root',
             'Gen_bbar_gsim-BKGx1-1M-4.root',
             'Gen_bbar_gsim-BKGx1-1M-5.root',
             'Gen_bbar_gsim-BKGx1-1M-6.root',
             'Gen_bbar_gsim-BKGx1-1M-7.root',
             'Gen_bbar_gsim-BKGx1-1M-8.root',
             'Gen_bbar_gsim-BKGx1-1M-9.root',
             'Gen_bbar_gsim-BKGx1-1M-10.root',
             ]
outputName = 'testSkimming_BBbar-10M-1.root'

# load files
inputMdstList(inputName)

# create lists of FSPs
pions = ('pi+:all', '')
kaons = ('K+:all', '')

fillParticleListsFromMC([pions, kaons], True)

# reconstruct phi and K_S0
reconstructDecay('phi:all -> K+:all K-:all', 'M < 1.3')
reconstructDecay('K_S0:all -> pi+:all pi-:all', '0.4 < M < 0.6')

# recostruct a B candidate out of those
reconstructDecay('B0:all -> phi:all K_S0:all', '5.0 < M < 5.5')


# add SkimFilter module to set condition variable based on the number of reconstructed B-tag mesons
skim = register_module('SkimFilter')
skim.param('particleLists', ['B0:all'])
analysis_main.add_module(skim)

# save event to a new mdst file if it contained at least one B-tag meson
empty_path = create_path()
skim.if_false(empty_path)
add_mdst_output(analysis_main, True, outputName)

# print out some further info
summaryOfLists(['B0:all'])

# process the events
process(analysis_main)

# print out the summary
print statistics
