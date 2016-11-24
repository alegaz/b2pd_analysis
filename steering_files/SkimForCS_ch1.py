#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################################
#
# Skim candidate events for
#
# B0 -> phi[K+K-] K0[pi+pi-]
#
# to be used for continuum suppression training
#
# Author: A. Gaz
#
######################################################

import sys

if len(sys.argv) != 2:
    sys.exit('\n Usage: basf2 SkimForCS_ch1.py `sample`.\n\n Where `sample` = `training` or `test`\n'
    )

import os
from basf2 import *
from modularAnalysis import *

sample = str(sys.argv[1])


# set the input files

filenameTrain = "/chai/sgt3/users/gaz/Belle2/PhiK0_rootfiles/TrainingSet/Ch1/*.root"
filelistTrain = [filenameTrain];

filenameTest = "/chai/sgt3/users/gaz/Belle2/PhiK0_rootfiles/MiniTestSet/Ch1/*.root"
filelistTest = [filenameTest];

outfile = 'SkimmedCS_ch1_'

if sample == 'training':
    inputMdstList('MC5', filelistTrain)
    outfile += 'train.root'

elif sample == 'test':
    inputMdstList('MC5', filelistTest)
    outfile += 'test.root'

else:
    sys.exit('Input sample does not match any of the availible samples: `training`, or `test`.'
             )


photons   = ('gamma:all',   '')
electrons = ('e-:all',      '')
muons     = ('mu-:all',     '')
pions     = ('pi-:all',     '')
kaons     = ('K-:all',      '')
protons   = ('anti-p-:all', '')

fillParticleLists([photons, electrons, muons, pions, kaons, protons])


# start reconstructing the candidates

fillParticleList('K_S0:mdst','0.3 < M < 0.7')
applyCuts('K_S0:mdst','0.475 < M < 0.520')

cutAndCopyList('gamma:good','gamma:all','0.060 < E < 6.000 and -150 < clusterTiming < 0 and clusterE9E25 > 0.75')
calibratePhotonEnergy('gamma:good', 0.0064)

reconstructDecay('phi:all -> K-:all K+:all','M < 1.1')
reconstructDecay('B0:ch1 -> phi:all K_S0:mdst','Mbc > 5.2 and abs(deltaE) < 0.2')
matchMCTruth('B0:ch1')


# train only on correctly reconstructed signal and continuum events:
applyCuts('B0:ch1','formula(isContinuumEvent+isSignal)>0.5')

# get the rest of the event:
buildRestOfEvent('B0:ch1')

# get tag vertex ('breco' is the type of MC association)
TagV('B0:ch1', 'breco')

# get continuum suppression (needed for flavor tagging)
buildContinuumSuppression('B0:ch1')


# variables for continuum suppression training:

from variables import variables

trainVars = [
    'R2',
    'cosTBTO',
    'KSFWVariables(hso02)',
    'KSFWVariables(hso12)',
    'cosTBz',
    'CleoCone(9)',
    'thrustBm',
    'thrustOm',
    'KSFWVariables(et)',
    'KSFWVariables(mm2)',
    'KSFWVariables(hso00)',
    'KSFWVariables(hso04)',
    'KSFWVariables(hso10)',
    'KSFWVariables(hso14)',
    'KSFWVariables(hso20)',
    'KSFWVariables(hso22)',
    'KSFWVariables(hso24)',
    'KSFWVariables(hoo0)',
    'KSFWVariables(hoo1)',
    'KSFWVariables(hoo2)',
    'KSFWVariables(hoo3)',
    'KSFWVariables(hoo4)',
    'CleoCone(1)',
    'CleoCone(2)',
    'CleoCone(3)',
    'CleoCone(4)',
    'CleoCone(5)',
    'CleoCone(6)',
    'CleoCone(7)',
    'CleoCone(8)',
    'CThetaB',
    ]

# Save target variable necessary for training.
targetVar = ['isNotContinuumEvent']
sigVar = ['isSignal']

# Create output file.
variablesToNTuple('B0:ch1', trainVars + targetVar + sigVar, treename='B0_ch1', filename=outfile)


# Process the events
process(analysis_main)

# print out the summary
print(statistics)
