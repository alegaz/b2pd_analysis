#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################################
#
# Skim candidate events for
#
# B0 -> phi[pi+pi-pi0] K0[pi+pi-]
#
# to be used for continuum suppression training
#
# Author: A. Gaz
#
######################################################

import sys

if len(sys.argv) != 2:
    sys.exit('\n Usage: basf2 SkimForCS_ch3.py `sample`.\n\n Where `sample` = `training` or `test`\n'
    )

import os
from basf2 import *
from modularAnalysis import *
from flavorTagger import *
from stdCharged import *
from stdPhotons import *
from stdPi0s import *

sample = str(sys.argv[1])


# set the input files

filenameTrain = "/chai/sgt3/users/gaz/Belle2/PhiK0_rootfiles/TrainingSet/Ch3/*.root"
filelistTrain = [filenameTrain];

filenameTest = "/chai/sgt3/users/gaz/Belle2/PhiK0_rootfiles/MiniTestSet/Ch3/*.root"
filelistTest = [filenameTest];

outfile = 'SkimmedCS_ch3_'

if sample == 'training':
    inputMdstList('MC5', filelistTrain)
    outfile += 'train.root'

elif sample == 'test':
    inputMdstList('MC5', filelistTest)
    outfile += 'test.root'

else:
    sys.exit('Input sample does not match any of the availible samples: `training`, or `test`.'
             )


loadStdCharged()
loadStdPhoton()
loadStdPi0()


# start reconstructing the candidates

fillParticleList('K_S0:mdst','0.3 < M < 0.7')
vertexKFit('K_S0:mdst', 0.0)
applyCuts('K_S0:mdst','0.475 < M < 0.520')
matchMCTruth('K_S0:mdst')

fillParticleList('gamma:pi0', '')
applyCuts('gamma:pi0', 'E > 0.05')

reconstructDecay('pi0:veryLooseFit -> gamma:pi0 gamma:pi0', '0.09 < M < 0.165')
massKFit('pi0:veryLooseFit', 0.0)
applyCuts('pi0:veryLooseFit', 'E > 0.2')
matchMCTruth('pi0:veryLooseFit')

reconstructDecay('phi:all -> pi0:veryLooseFit pi+:all pi-:all','0.95< M < 1.07')
vertexKFit('phi:all', 0.0)
matchMCTruth('phi:all')

reconstructDecay('B0:ch3 -> phi:all K_S0:mdst','Mbc > 5.2 and abs(deltaE) < 0.2')
vertexRave('B0:ch3', 0.0, 'B0:ch3 -> [phi -> ^pi+ ^pi- pi0] ^K_S0', 'iptube')
matchMCTruth('B0:ch3')


# train only on correctly reconstructed signal and continuum events:
applyCuts('B0:ch3','formula(isContinuumEvent+isSignal)>0.5')

# get the rest of the event:
buildRestOfEvent('B0:ch3')

# get tag vertex ('breco' is the type of MC association)
TagV('B0:ch3', 'breco')

# get continuum suppression (needed for flavor tagging)
buildContinuumSuppression('B0:ch3')


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
variablesToNTuple('B0:ch3', trainVars + targetVar + sigVar, treename='B0_ch3', filename=outfile)


# Process the events
process(analysis_main)

# print out the summary
print(statistics)
