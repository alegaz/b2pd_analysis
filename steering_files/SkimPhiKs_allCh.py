#!/usr/bin/env python
# -*- coding: utf-8 -*-

from basf2 import *
from reconstruction import *
from modularAnalysis import *
from stdCharged import *
from stdPhotons import *
from stdPi0s import *
import sys


outputName = 'testSkim.root'

fileNames = "/chai/sgt3/users/gaz/Belle2/PhiK0_rootfiles/MC9/SkimSet/mdst_sample_ch2.root"
fileList = [fileNames]
inputMdstList('default', fileList)


loadStdCharged()
stdPhotons('pi0')


# reconstruct candidates

reconstructDecay('pi0:veryLoose -> gamma:pi0 gamma:pi0', '0.09 < M < 0.165')
applyCuts('pi0:veryLoose', 'E > 0.2')

fillParticleList('K_S0:mdst','0.4 < M < 0.6')
reconstructDecay('K_S0:neu -> pi0:veryLoose pi0:veryLoose', '0.4 < M < 0.6')

reconstructDecay('phi:KK -> K+:all K-:all', 'M < 1.3')
reconstructDecay('phi:3pi -> pi0:veryLoose pi+:all pi-:all', 'M < 1.3')

reconstructDecay('B0:ch1 -> phi:KK K_S0:mdst', 'Mbc > 5.2 and abs(deltaE) < 0.2', 1)
reconstructDecay('B0:ch2 -> phi:KK K_S0:neu', 'Mbc > 5.0 and abs(deltaE) < 0.3', 2)
reconstructDecay('B0:ch3 -> phi:3pi K_S0:mdst', 'Mbc > 5.0 and abs(deltaE) < 0.3', 3)

copyLists('B0:all', ['B0:ch1', 'B0:ch2', 'B0:ch3'])


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
print (statistics)
