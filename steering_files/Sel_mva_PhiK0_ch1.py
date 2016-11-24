#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################################
#
# Author: A. Gaz
#
######################################################

import sys

if len(sys.argv) != 5:
    sys.exit('\n Usage: basf2 SelectPhiK0_ch1.py `action` `sample` `file_name_stub` `output file name`.\n\n Where `action` = `simple`, or `expert`\n and `sample` = `signal`, `BBbar` or `continuum`\n'
    )

import os
from basf2 import *
from modularAnalysis import *

action = str(sys.argv[1])
sample = str(sys.argv[2])
stub   = str(sys.argv[3])
outFile = str(sys.argv[4])


# set the input files

filenameSIG = "/chai/sgt3/users/gaz/Belle2/PhiK0_rootfiles/SignalMC/Ch1/private/"
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
    inputMdstList('MC5', filelistSIG)
            
elif sample == 'BBbar':
    inputMdstList('MC5', filelistBBbar)
            
elif sample == 'continuum':
    inputMdstList('MC5', filelistCC)
        
else:
    sys.exit('Input sample does not match any of the availible samples: `signal`, `BBBar` or `continuum`.'
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
vertexKFit('K_S0:mdst', 0.0)
applyCuts('K_S0:mdst','0.475 < M < 0.520')
matchMCTruth('K_S0:mdst')

cutAndCopyList('gamma:good','gamma:all','0.060 < E < 6.000 and -150 < clusterTiming < 0 and clusterE9E25 > 0.75')
calibratePhotonEnergy('gamma:good', 0.0064)

reconstructDecay('pi0:all -> gamma:good gamma:good',' 0.080 < M < 0.200')
massKFit('pi0:all',0.0)
matchMCTruth('pi0:all')

reconstructDecay('phi:all -> K-:all K+:all','M < 1.1')
vertexKFit('phi:all', 0.0)
matchMCTruth('phi:all')

reconstructDecay('B0:ch1 -> phi:all K_S0:mdst','Mbc > 5.2 and abs(deltaE) < 0.2')
vertexRave('B0:ch1', 0.0, 'B0:ch1 -> [phi -> ^K+ ^K-] ^K_S0', 'iptube')
matchMCTruth('B0:ch1')

# get the rest of the event:
buildRestOfEvent('B0:ch1')

# get tag vertex ('breco' is the type of MC association)
TagV('B0:ch1', 'breco')

# get continuum suppression (needed for flavor tagging)
buildContinuumSuppression('B0:ch1')


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
    analysis_main.add_module('MVAExpert', listNames=['B0:ch1'], extraInfoName='FastBDT', identifier='MVAFastBDT_ch1.root')

    # Variables from MVAExpert.
    #expertVars = ['extraInfo(FastBDT)', 'transformedNetworkOutput(FastBDT,0.,1.0)']



# flavor tagging
#FlavorTagger(
#    mode='Expert',
#    weightFiles='B2JpsiKs_mu',
#    categories=['Electron', 'Muon', 'KinLepton', 'Kaon', 'SlowPion', 'FastPion', 'Lambda', 'FSC', 'MaximumP*', 'KaonPion'])


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
# for the real analysis, anything included in B0_ch1
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

toolsPi0 = ['EventMetaData', 'pi0:all']
toolsPi0 += ['InvMass', '^pi0:all -> gamma gamma']
toolsPi0 += ['Kinematics', '^pi0:all -> ^gamma ^gamma']
toolsPi0 += ['MCKinematics', '^pi0:all -> ^gamma ^gamma']
toolsPi0 += ['MCTruth', '^pi0:all -> ^gamma ^gamma']
toolsPi0 += ['MCHierarchy', '^pi0:all -> ^gamma ^gamma']

toolsKsChg = ['InvMass', '^K_S0:mdst -> pi+ pi-']
toolsKsChg += ['Vertex', '^K_S0:mdst -> pi+ pi-']
toolsKsChg += ['EventMetaData', '^K_S0:mdst']

toolsPhiKK = ['InvMass', '^phi -> K+ K-']
toolsPhiKK += ['Vertex', '^phi -> K+ K-']
toolsPhiKK += ['EventMetaData', '^phi']
toolsPhiKK += ['Kinematics',               '^phi -> ^K- ^K+']
toolsPhiKK += ['MCTruth',                  '^phi']
toolsPhiKK += ['Vertex',                   '^phi']
toolsPhiKK += ['MCReconstructible',        'phi -> ^K- ^K+']
toolsPhiKK += ['DeltaLogL',                'phi -> ^K- ^K+']
toolsPhiKK += ['Track',                    'phi -> ^K- ^K+']
toolsPhiKK += ['CustomFloats[isSignal:M]', '^phi']

toolsBsigCh1 = ['EventMetaData', '^B0:ch1']
toolsBsigCh1 += ['InvMass', '^B0:ch1 -> [^phi -> K+ K-] ^K_S0']
toolsBsigCh1 += ['Charge', 'B0:ch1 -> [phi -> ^K+ ^K-] [K_S0 -> ^pi+ ^pi-]']
toolsBsigCh1 += ['Kinematics', '^B0:ch1 -> [^phi -> ^K+ ^K-] [^K_S0 -> ^pi+ ^pi-]']
toolsBsigCh1 += ['CMSKinematics', '^B0:ch1 -> [^phi -> ^K+ ^K-] [^K_S0 -> ^pi+ ^pi-]']
toolsBsigCh1 += ['DeltaEMbc', '^B0:ch1']
toolsBsigCh1 += ['PID', 'B0:ch1 -> [phi -> ^K+ ^K-] [K_S0 -> ^pi+ ^pi-]']
toolsBsigCh1 += ['Track', 'B0:ch1 -> [phi -> ^K+ ^K-] [K_S0 -> ^pi+ ^pi-]']
toolsBsigCh1 += ['TrackHits', 'B0:ch1 -> [phi -> ^K+ ^K-] [K_S0 -> ^pi+ ^pi-]']
toolsBsigCh1 += ['MCTruth', '^B0:ch1 -> [^phi -> ^K+ ^K-] [^K_S0 -> ^pi+ ^pi-]']
toolsBsigCh1 += ['MCHierarchy', '^B0:ch1 -> [^phi -> ^K+ ^K-] [^K_S0 -> ^pi+ ^pi-]']
toolsBsigCh1 += ['MCKinematics', '^B0:ch1 -> [^phi -> ^K+ ^K-] [^K_S0 -> ^pi+ ^pi-]']
toolsBsigCh1 += ['Vertex', '^B0:ch1 -> ^phi ^K_S0']
toolsBsigCh1 += ['MCVertex', '^B0:ch1 -> ^phi ^K_S0']
toolsBsigCh1 += ['TagVertex', '^B0:ch1']
toolsBsigCh1 += ['MCTagVertex', '^B0:ch1']
toolsBsigCh1 += ['DeltaT', '^B0:ch1']
toolsBsigCh1 += ['MCDeltaT', '^B0:ch1']
toolsBsigCh1 += ['CustomFloats[isSignal:isContinuumEvent]', '^B0:ch1']
toolsBsigCh1 += ['CustomFloats[flLenSig]', 'B0:ch1 -> phi ^K_S0']
toolsBsigCh1 += ['CustomFloats[CThetaB]', '^B0:ch1']
toolsBsigCh1 += ['CustomFloats[CosTHel]', 'B0:ch1 -> ^phi K_S0']
toolsBsigCh1 += ['ContinuumSuppression', '^B0:ch1']
toolsBsigCh1 += ['ContinuumSuppression[FS1]', '^B0:ch1']

if action == 'expert':

    toolsBsigCh1 += ['CustomFloats[out_FastBDT]', '^B0:ch1']
    toolsBsigCh1 += ['CustomFloats[csv_FastBDT]', '^B0:ch1']


#toolsBsigCh1 += ['FlavorTagging', '^B0:ch1']

toolsRS = ['RecoStats', '^B0:ch1']


# save stuff to root file
ntupleFile(outFile)

#ntupleTree('Trks', 'pi+:all', toolsTrk)
#ntupleTree('Neu', 'gamma:all', toolsNeu)
#ntupleTree('Pi0', 'pi0:all', toolsPi0)
#ntupleTree('KsChg', 'K_S0:mdst', toolsKsChg)
#ntupleTree('PhiKK', 'phi:all', toolsPhiKK)
ntupleTree('B0_ch1', 'B0:ch1', toolsBsigCh1)
ntupleTree('RecoStats', 'B0:ch1', toolsRS)

# Process the events
process(analysis_main)

# print out the summary
print(statistics)
