#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################################
#
# Author: A. Gaz
#
######################################################

import sys

if len(sys.argv) != 5:
    sys.exit('\n Usage: basf2 SelectPhiK0_ch3.py `action` `sample` `file_name_stub` `output file name`.\n\n Where `action` = `simple`, or `expert`\n and `sample` = `signal`, `BBbar` or `continuum`\n'
    )

import os
from basf2 import *
from modularAnalysis import *
from flavorTagger import *
from stdCharged import *
from stdPhotons import *
from stdPi0s import *

action = str(sys.argv[1])
sample = str(sys.argv[2])
stub   = str(sys.argv[3])
outFile = str(sys.argv[4])


# set the input files

filenameSIG = "/chai/sgt3/users/gaz/Belle2/PhiK0_rootfiles/MC7/Ch3/"
filenameSIG += stub
filenameSIG += "*.root"
filelistSIG = [filenameSIG];

filenameCC = "/chai/sgt3/users/gaz/Belle2/PhiK0_rootfiles/GenericMC/continuum/"
filenameCC += stub
filenameCC += "*.root"
filelistCC = [filenameCC];

filenameBBbar = "/chai/sgt3/users/gaz/Belle2/PhiK0_rootfiles/GenericMC/BBbar/"
filenameBBbar += stub
filenameBBbar += "*.root"
filelistBBbar = [filenameBBbar];


if sample == 'signal':
    inputMdstList('default', filelistSIG)
            
elif sample == 'BBbar':
    inputMdstList('MC5', filelistBBbar)
            
elif sample == 'continuum':
    inputMdstList('MC5', filelistCC)
        
else:
    sys.exit('Input sample does not match any of the availible samples: `signal`, `BBBar` or `continuum`.'
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

# get the rest of the event:
buildRestOfEvent('B0:ch3')

# flavor tagging
flavorTagger(
    mode='Expert',
    particleList='B0:ch3',
    combinerMethods=['TMVA-FBDT', 'FANN-MLP'],
    workingDirectory="./JpsiMuMu_ftTrain/BGx1",
    belleOrBelle2='Belle2')

# get tag vertex ('breco' is the type of MC association)
TagV('B0:ch3', 'breco')

# get continuum suppression (needed for flavor tagging)
buildContinuumSuppression('B0:ch3')


if action == 'expert':

    # Variables for training.
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

    # Target variable used in training.
    targetVar = ['isNotContinuumEvent']

    # MVAExpert
    analysis_main.add_module('MVAExpert', listNames=['B0:ch3'], extraInfoName='FastBDT', identifier='MVAFastBDT_ch3.root')

    # Variables from MVAExpert.
    expertVars = ['extraInfo(FastBDT)', 'transformedNetworkOutput(FastBDT,0.,1.0)']



# variable aliases

from variables import variables

variables.addAlias('flLenSig', 'significanceOfDistance')
variables.addAlias('CosTHel', 'decayAngle(0)')
variables.addAlias('out_FastBDT',  'extraInfo(FastBDT)')
variables.addAlias('csv_FastBDT',  'transformedNetworkOutput(FastBDT,0.,1.0)')


# define Ntuple tools 

# with tracks, neutrals, pi0's and intermediate resonances
# are meant only for debugging
#
# for the real analysis, anything included in B0_ch3
# should be sufficient

toolsTrk = ['EventMetaData', 'pi+:all']
toolsTrk += ['Kinematics', '^pi+:all']
toolsTrk += ['MCKinematics', '^pi+:all']
toolsTrk += ['MCTruth','^pi+:all']
toolsTrk += ['MCHierarchy', '^pi+:all']
toolsTrk += ['CustomFloats[isContinuumEvent]', '^pi+:all']
toolsTrk += ['PID', '^pi+:all']
toolsTrk += ['Track', '^pi+:all']

toolsNeu = ['EventMetaData', 'gamma:all']
toolsNeu += ['Kinematics', '^gamma:all']
toolsNeu += ['MCKinematics', '^gamma:all']
toolsNeu += ['MCTruth', '^gamma:all']
toolsNeu += ['MCHierarchy', '^gamma:all']

toolsPi0 = ['EventMetaData', 'pi0:veryLooseFit']
toolsPi0 += ['InvMass', '^pi0:veryLooseFit -> gamma gamma']
toolsPi0 += ['Kinematics', '^pi0:veryLooseFit -> ^gamma ^gamma']
toolsPi0 += ['MCKinematics', '^pi0:veryLooseFit -> ^gamma ^gamma']
toolsPi0 += ['MCTruth', '^pi0:veryLooseFit -> ^gamma ^gamma']
toolsPi0 += ['MCHierarchy', '^pi0:veryLooseFit -> ^gamma ^gamma']

toolsKsChg = ['InvMass', '^K_S0:mdst -> pi+ pi-']
toolsKsChg += ['Vertex', '^K_S0:mdst -> pi+ pi-']
toolsKsChg += ['EventMetaData', '^K_S0:mdst']

toolsPhi3Pi = ['InvMass', '^phi -> ^pi0 pi+ pi-']
toolsPhi3Pi += ['Vertex', '^phi -> pi0 pi+ pi-']
toolsPhi3Pi += ['EventMetaData', '^phi']
toolsPhi3Pi += ['Kinematics',               '^phi -> ^pi0 ^pi+ ^pi-']
toolsPhi3Pi += ['MCTruth',                  '^phi']
toolsPhi3Pi += ['Vertex',                   '^phi']
toolsPhi3Pi += ['MCReconstructible',        'phi -> ^pi0 ^pi+ ^pi-']
toolsPhi3Pi += ['DeltaLogL',                'phi -> ^pi0 ^pi+ ^pi-']
toolsPhi3Pi += ['Track',                    'phi -> ^pi0 ^pi+ ^pi-']
toolsPhi3Pi += ['CustomFloats[isSignal:M]', '^phi']

toolsBsigCh3 = ['EventMetaData', '^B0:ch3']
toolsBsigCh3 += ['InvMass', '^B0:ch3 -> [^phi -> pi0 pi+ pi-] ^K_S0']
toolsBsigCh3 += ['InvMass[BeforeFit]', 'B0:ch3 -> [phi -> ^pi0 pi+ pi-] K_S0']
toolsBsigCh3 += ['Charge', 'B0:ch3 -> [phi -> pi0 ^pi+ ^pi-] [K_S0 -> ^pi+ ^pi-]']
toolsBsigCh3 += ['Kinematics', '^B0:ch3 -> [^phi -> [^pi0 -> ^gamma ^gamma] ^pi+ ^pi-] [^K_S0 -> ^pi+ ^pi-]']
toolsBsigCh3 += ['CMSKinematics', '^B0:ch3 -> [^phi -> [^pi0 -> ^gamma ^gamma] ^pi+ ^pi-] [^K_S0 -> ^pi+ ^pi-]']
toolsBsigCh3 += ['DeltaEMbc', '^B0:ch3']
toolsBsigCh3 += ['PID', 'B0:ch3 -> [phi -> pi0 ^pi+ ^pi-] [K_S0 -> ^pi+ ^pi-]']
toolsBsigCh3 += ['Track', 'B0:ch3 -> [phi -> pi0 ^pi+ ^pi-] [K_S0 -> ^pi+ ^pi-]']
toolsBsigCh3 += ['TrackHits', 'B0:ch3 -> [phi -> pi0 ^pi+ ^pi-] [K_S0 -> ^pi+ ^pi-]']
toolsBsigCh3 += ['MCTruth', '^B0:ch3 -> [^phi -> ^pi0 ^pi+ ^pi-] [^K_S0 -> ^pi+ ^pi-]']
toolsBsigCh3 += ['MCHierarchy', '^B0:ch3 -> [^phi -> [^pi0 -> ^gamma ^gamma] ^pi+ ^pi-] [^K_S0 -> ^pi+ ^pi-]']
toolsBsigCh3 += ['MCKinematics', '^B0:ch3 -> [^phi -> [^pi0 -> ^gamma ^gamma] ^pi+ ^pi-] [^K_S0 -> ^pi+ ^pi-]']
toolsBsigCh3 += ['Vertex', '^B0:ch3 -> ^phi ^K_S0']
toolsBsigCh3 += ['MCVertex', '^B0:ch3 -> ^phi ^K_S0']
toolsBsigCh3 += ['TagVertex', '^B0:ch3']
toolsBsigCh3 += ['MCTagVertex', '^B0:ch3']
toolsBsigCh3 += ['DeltaT', '^B0:ch3']
toolsBsigCh3 += ['MCDeltaT', '^B0:ch3']
toolsBsigCh3 += ['CustomFloats[isSignal]', '^B0:ch3 -> ^phi ^K_S0']
toolsBsigCh3 += ['CustomFloats[isContinuumEvent]', '^B0:ch3']
toolsBsigCh3 += ['CustomFloats[flLenSig]', 'B0:ch3 -> phi ^K_S0']
toolsBsigCh3 += ['CustomFloats[CThetaB]', '^B0:ch3']
toolsBsigCh3 += ['CustomFloats[CThHelVto3P]', 'B0:ch3 -> ^phi K_S0']
toolsBsigCh3 += ['ContinuumSuppression', '^B0:ch3']
toolsBsigCh3 += ['ContinuumSuppression[FS1]', '^B0:ch3']
toolsBsigCh3 += ['FlavorTagging[TMVA-FBDT, FANN-MLP, qrCategories]', '^B0:ch3']

if action == 'expert':

    toolsBsigCh3 += ['CustomFloats[out_FastBDT]', '^B0:ch3']
    toolsBsigCh3 += ['CustomFloats[csv_FastBDT]', '^B0:ch3']



toolsRS = ['RecoStats', '^B0:ch3']


# save stuff to root file
ntupleFile(outFile)

#ntupleTree('Trks', 'pi+:all', toolsTrk)
#ntupleTree('Neu', 'gamma:all', toolsNeu)
#ntupleTree('Pi0', 'pi0:veryLooseFit', toolsPi0)
#ntupleTree('KsChg', 'K_S0:mdst', toolsKsChg)
#ntupleTree('Phi3Pi', 'phi:all', toolsPhi3Pi)
ntupleTree('B0_ch3', 'B0:ch3', toolsBsigCh3)
ntupleTree('RecoStats', 'B0:ch3', toolsRS)

# Process the events
process(analysis_main)

# print out the summary
print(statistics)
