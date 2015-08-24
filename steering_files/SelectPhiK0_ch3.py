#!/usr/bin/env python
# -*- coding: utf-8 -*-

######################################################
#
# Get started with B0 -> phi K0
#
# These are the channels that are going to be 
# investigated:
#
# 1: phi[KK] KS[pi+pi-]
# 2: phi[KK] KS[pi0pi0]
# 3: phi[pi+pi-pi0] KS[pi+pi-]
# 4: phi[KK] KL
#
# Author: A. Gaz
#
######################################################

import sys

if len(sys.argv) != 3:
    sys.exit('\n Usage: basf2 SelectPhiK0_ch3.py `action` `sample`.\n\n Where `action` = `simple`, `train`, or `expert`\n and `sample` = `signal`, `BBbar` or `continuum`\n'
    )

from basf2 import *
from modularAnalysis import *


action = str(sys.argv[1])
sample = str(sys.argv[2])
outFile = ''


# set the input files

filelistSIG = ['../saved_rootfiles/B0_phi-pi+pi-pi0_KS-pi+pi-_gsim-BKGx1-50000-*.root']
#filelistSIG = ['B0_phi-pi+pi-pi0_KS-pi+pi-_gsim-BKGx1-50000.root']
filelistBBbar = ['../saved_rootfiles/reconstructedBBbar-10M.root']
filelistCC = ['../saved_rootfiles/reconstructedContinuum_ch3_1M_*.root']


# define input and output files

outFile = 'B0_PhiK0_ch3_output.root'

if action == 'training':
    inputMdstList(filelistSIG+filelistCC)

else:
    if sample == 'signal':
        inputMdstList(filelistSIG)
        if action == 'simple':
            outFile = 'B0_PhiK0_ch3_output_signal.root'
        elif action == 'expert':
            outFile = 'B0_PhiK0_ch3_output_cs_signal.root'
            
    elif sample == 'BBbar':
        inputMdstList(filelistBBbar)
        if action == 'simple':
            outFile = 'B0_PhiK0_ch3_output_BBbar.root'
        elif action == 'expert':
            outFile = 'B0_PhiK0_ch3_output_cs_BBbar.root'
            
    elif sample == 'continuum':
        inputMdstList(filelistCC)
        if action == 'simple':
            outFile = 'B0_PhiK0_ch3_output_continuum.root'
        elif action == 'expert':
            outFile = 'B0_PhiK0_ch3_output_cs_continuum.root'
        
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

reconstructDecay('phi:all -> pi0:all pi+:all pi-:all','0.95< M < 1.07')
vertexKFit('phi:all', 0.0)
matchMCTruth('phi:all')

reconstructDecay('B0:ch3 -> phi:all K_S0:mdst','Mbc > 5.2 and abs(deltaE) < 0.2')
vertexRave('B0:ch3', 0.0, 'B0:ch3 -> [phi -> pi0 ^pi+ ^pi-] K_S0')
matchMCTruth('B0:ch3')

# get the rest of the event:
buildRestOfEvent('B0:ch3')

# get tag vertex ('breco' is the type of MC association)
TagV('B0:ch3', 'breco')

# get continuum suppression (needed for flavor tagging)
buildContinuumSuppression('B0:ch3')


if action == 'training':

    # Define the input variables.
    variables = [
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
        ]

    # Define the methods.
    methods = [('FastBDT', 'Plugin',
                'H:V:CreateMVAPdfs:NTrees=100:Shrinkage=0.10:RandRatio=0.5:NCutLevel=8:NTreeLayers=3')]

    # TMVA training/testing
    teacher = register_module('TMVATeacher')
    teacher.param('prefix', 'B0_PhiKs_ch3_TMVA')
    teacher.param('methods', methods)
    teacher.param('variables', variables)
    teacher.param('target', 'isContinuumEvent')
    teacher.param('listNames', ['B0:ch3'])
    teacher.param('workingDirectory', 'training')
    analysis_main.add_module(teacher)


if action == 'expert':
    # run the expert mode
    methods = ['FastBDT']

    for method in methods:
        expert = register_module('TMVAExpert')
        expert.param('prefix', 'B0_PhiKs_ch3_TMVA')
        expert.param('method', method)
        expert.param('listNames', ['B0:ch3'])
        expert.param('expertOutputName', method + '_Probability')
        expert.param('workingDirectory', 'training')
        analysis_main.add_module(expert)

    # Network output
    networkOutput = ['extraInfo({method}_Probability)'.format(method=method)
                     for method in methods]
    #transformedNetworkOutputNB = \
    #    ['transformedNetworkOutput(NeuroBayes_Probability,-0.9,1.0)']
    transformedNetworkOutputFBDT = \
                                 ['transformedNetworkOutput(FastBDT_Probability,0.0,1.0)']





# flavor tagging
#FlavorTagger(
#    mode='Expert',
#    weightFiles='B2JpsiKs_mu',
#    categories=['Electron', 'Muon', 'KinLepton', 'Kaon', 'SlowPion', 'FastPion', 'Lambda', 'FSC', 'MaximumP*', 'KaonPion'])

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
toolsBsigCh3 += ['InvMass', '^B0:ch3 -> [^phi -> ^pi0 pi+ pi-] ^K_S0']
toolsBsigCh3 += ['Charge', 'B0:ch3 -> [phi -> pi0 ^pi+ ^pi-] [K_S0 -> ^pi+ ^pi-]']
toolsBsigCh3 += ['Kinematics', '^B0:ch3 -> [^phi -> ^pi0 ^pi+ ^pi-] [^K_S0 -> ^pi+ ^pi-]']
toolsBsigCh3 += ['CMSKinematics', '^B0:ch3 -> [^phi -> ^pi0 ^pi+ ^pi-] [^K_S0 -> ^pi+ ^pi-]']
toolsBsigCh3 += ['DeltaEMbc', '^B0:ch3']
toolsBsigCh3 += ['PID', 'B0:ch3 -> [phi -> pi0 ^pi+ ^pi-] [K_S0 -> ^pi+ ^pi-]']
toolsBsigCh3 += ['Track', 'B0:ch3 -> [phi -> pi0 ^pi+ ^pi-] [K_S0 -> ^pi+ ^pi-]']
toolsBsigCh3 += ['TrackHits', 'B0:ch3 -> [phi -> pi0 ^pi+ ^pi-] [K_S0 -> ^pi+ ^pi-]']
toolsBsigCh3 += ['MCTruth', '^B0:ch3 -> [^phi -> ^pi0 ^pi+ ^pi-] [^K_S0 -> ^pi+ ^pi-]']
toolsBsigCh3 += ['MCHierarchy', '^B0:ch3 -> [^phi -> ^pi0 ^pi+ ^pi-] [^K_S0 -> ^pi+ ^pi-]']
toolsBsigCh3 += ['MCKinematics', '^B0:ch3 -> [^phi -> ^pi0 ^pi+ ^pi-] [^K_S0 -> ^pi+ ^pi-]']
toolsBsigCh3 += ['Vertex', '^B0:ch3 -> ^phi ^K_S0']
toolsBsigCh3 += ['MCVertex', '^B0:ch3 -> ^phi ^K_S0']
toolsBsigCh3 += ['TagVertex', '^B0:ch3']
toolsBsigCh3 += ['MCTagVertex', '^B0:ch3']
toolsBsigCh3 += ['DeltaT', '^B0:ch3']
toolsBsigCh3 += ['MCDeltaT', '^B0:ch3']
toolsBsigCh3 += ['CustomFloats[isSignal:isContinuumEvent]', '^B0:ch3']

if action == 'expert':
    toolsBsigCh3 += ['ContinuumSuppression', '^B0:ch3']
    toolsBsigCh3 += ['ContinuumSuppression[FS1]', '^B0:ch3']
    toolsBsigCh3 += ['CustomFloats[' + ':'.join(networkOutput) + ']', '^B0:ch3']
    toolsBsigCh3 += ['CustomFloats[' + ':'.join(transformedNetworkOutputFBDT) + ']', '^B0:ch3']
    #toolsBsigCh3 += ['CustomFloats[' + ':'.join(transformedNetworkOutputNB) + ']', '^B0:ch3']

#toolsBsigCh3 += ['FlavorTagging', '^B0:ch3']

toolsRS = ['RecoStats', '^B0:ch3']


# save stuff to root file
ntupleFile(outFile)

ntupleTree('Trks', 'pi+:all', toolsTrk)
ntupleTree('Neu', 'gamma:all', toolsNeu)
ntupleTree('Pi0', 'pi0:all', toolsPi0)
ntupleTree('KsChg', 'K_S0:mdst', toolsKsChg)
ntupleTree('Phi3Pi', 'phi:all', toolsPhi3Pi)
ntupleTree('B0_ch3', 'B0:ch3', toolsBsigCh3)
ntupleTree('RecoStats', 'B0:ch3', toolsRS)

# Process the events
process(analysis_main)

# print out the summary
print statistics
