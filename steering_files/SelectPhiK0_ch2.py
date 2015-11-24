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
    sys.exit('\n Usage: basf2 SelectPhiK0_ch2.py `action` `sample`.\n\n Where `action` = `simple`, `train`, or `expert`\n and `sample` = `signal`, `BBbar` or `continuum`\n'
    )

from basf2 import *
from modularAnalysis import *


action = str(sys.argv[1])
sample = str(sys.argv[2])
outFile = ''


# set the input files

#filelistSIG = ['../saved_rootfiles/B0_Phi-KK_Ks-pi0pi0_gsim-BKGx1-50000-1.root']
filelistSIG = ['test_signal_ch2_v9.root']
filelistCC = ['/ghi/fs01/belle2/bdata/MC/fab/sim/release-00-05-03/DBxxxxxxxx/MC5/prod00000005/s00/e0001/4S/r00001/ssbar/sub00/mdst_00001*_prod00000005_task0000001*.root']
filelistBBbar = ['../saved_rootfiles/reconstructedBBbar-10M.root']


# define input and output files

outFile = 'B0_PhiK0_ch2_output.root'

if action == 'training':
    inputMdstList(filelistSIG+filelistCC)

else:
    if sample == 'signal':
        inputMdstList(filelistSIG)
        if action == 'simple':
            outFile = 'B0_PhiK0_ch2_output_signal.root'
        elif action == 'expert':
            outFile = 'B0_PhiK0_ch2_output_cs_signal.root'
            
    elif sample == 'BBbar':
        inputMdstList(filelistCC)
        if action == 'simple':
            outFile = 'B0_PhiK0_ch2_output_BBbar.root'
        elif action == 'expert':
            outFile = 'B0_PhiK0_ch2_output_cs_BBbar.root'
            
    elif sample == 'continuum':
        inputMdstList(filelistCC)
        if action == 'simple':
            outFile = 'B0_PhiK0_ch2_output_continuum.root'
        elif action == 'expert':
            outFile = 'B0_PhiK0_ch2_output_cs_continuum.root'
        
    else:
        sys.exit('Input sample does not match any of the availible samples: `signal`, `BBbar` or `continuum`.'
                 )


photons   = ('gamma:all',   '')
electrons = ('e-:all',      '')
muons     = ('mu-:all',     '')
pions     = ('pi-:all',     '')
kaons     = ('K-:all',      '')
protons   = ('anti-p-:all', '')

fillParticleLists([photons, electrons, muons, pions, kaons, protons])


# start reconstructing the candidates

cutAndCopyList('gamma:good','gamma:all','0.060 < E < 6.000 and -150 < clusterTiming < 0 and clusterE9E25 > 0.75')
calibratePhotonEnergy('gamma:good', 0.0064)

reconstructDecay('pi0:all -> gamma:good gamma:good','0.080 < M < 0.200')
massKFit('pi0:all',0.0)
matchMCTruth('pi0:all')

reconstructDecay('K_S0:neu -> pi0:all pi0:all', '0.4 < M < 0.6')
massKFit('K_S0:neu',0.0)
matchMCTruth('K_S0:neu')

reconstructDecay('phi:all -> K-:all K+:all','M < 1.1')
vertexKFit('phi:all', 0.0)
matchMCTruth('phi:all')

reconstructDecay('B0:ch2 -> phi:all K_S0:neu','Mbc > 5.0 and abs(deltaE) < 0.5')
vertexRave('B0:ch2', 0.0, 'B0:ch2 -> [phi -> ^K+ ^K-] K_S0')
matchMCTruth('B0:ch2')

# get the rest of the event:
buildRestOfEvent('B0:ch2')

# get tag vertex ('breco' is the type of MC association)
TagV('B0:ch2', 'breco')

# get continuum suppression (needed for flavor tagging)
buildContinuumSuppression('B0:ch2')


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
    methods = [('LPCA', 'Likelihood',
                'H:V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA')]

    # TMVA training/testing
    teacher = register_module('TMVATeacher')
    teacher.param('prefix', 'B0_PhiKs_ch2_TMVA')
    teacher.param('methods', methods)
    teacher.param('variables', variables)
    teacher.param('target', 'isContinuumEvent')
    teacher.param('listNames', ['B0:ch2'])
    teacher.param('workingDirectory', 'training')
    analysis_main.add_module(teacher)


if action == 'expert':
    # run the expert mode
    methods = ['LPCA']

    for method in methods:
        expert = register_module('TMVAExpert')
        expert.param('prefix', 'B0_PhiKs_ch2_TMVA')
        expert.param('method', method)
        expert.param('listNames', ['B0:ch2'])
        expert.param('expertOutputName', method + '_Probability')
        expert.param('workingDirectory', 'training')
        analysis_main.add_module(expert)

    # Network output
    networkOutput = ['extraInfo({method}_Probability)'.format(method=method)
                     for method in methods]
    transformedNetworkOutputLPCA = \
        ['transformedNetworkOutput(LPCA_Probability,0.,1.0)']
    #transformedNetworkOutputNB = \
    #    ['transformedNetworkOutput(NeuroBayes_Probability,-0.9,1.0)']
    #transformedNetworkOutputFBDT = \
    #                             ['transformedNetworkOutput(FastBDT_Probability,0.0,1.0)']





# flavor tagging
#FlavorTagger(
#    mode='Expert',
#    weightFiles='B2JpsiKs_mu',
#    categories=['Electron', 'Muon', 'KinLepton', 'Kaon', 'SlowPion', 'FastPion', 'Lambda', 'FSC', 'MaximumP*', 'KaonPion'])

# define Ntuple tools 

# with tracks, neutrals, pi0's and intermediate resonances
# are meant only for debugging
#
# for the real analysis, anything included in B0_ch2
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

toolsKsNeu = ['InvMass', '^K_S0:neu -> ^pi0 ^pi0']
toolsKsNeu += ['Vertex', '^K_S0:neu -> ^pi0 ^pi0']
toolsKsNeu += ['EventMetaData', '^K_S0']

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

toolsBsigCh2 = ['EventMetaData', '^B0:ch2']
toolsBsigCh2 += ['InvMass', '^B0:ch2 -> [^phi -> K+ K-] K_S0']
toolsBsigCh2 += ['InvMass[BeforeFit]', 'B0:ch2 -> [phi -> K+ K-] [^K_S0 -> ^pi0 ^pi0]']
toolsBsigCh2 += ['Charge', 'B0:ch2 -> [phi -> ^K+ ^K-] [K_S0 -> pi0 pi0]']
toolsBsigCh2 += ['Kinematics', '^B0:ch2 -> [^phi -> ^K+ ^K-] [^K_S0 -> ^pi0 ^pi0]']
toolsBsigCh2 += ['CMSKinematics', '^B0:ch2 -> [^phi -> ^K+ ^K-] [^K_S0 -> ^pi0 ^pi0]']
toolsBsigCh2 += ['DeltaEMbc', '^B0:ch2']
toolsBsigCh2 += ['PID', 'B0:ch2 -> [phi -> ^K+ ^K-] [K_S0 -> pi0 pi0]']
toolsBsigCh2 += ['Track', 'B0:ch2 -> [phi -> ^K+ ^K-] [K_S0 -> pi0 pi0]']
toolsBsigCh2 += ['TrackHits', 'B0:ch2 -> [phi -> ^K+ ^K-] [K_S0 -> pi0 pi0]']
toolsBsigCh2 += ['MCTruth', '^B0:ch2 -> [^phi -> ^K+ ^K-] [^K_S0 -> ^pi0 ^pi0]']
toolsBsigCh2 += ['MCHierarchy', '^B0:ch2 -> [^phi -> ^K+ ^K-] [^K_S0 -> ^pi0 ^pi0]']
toolsBsigCh2 += ['MCKinematics', '^B0:ch2 -> [^phi -> ^K+ ^K-] [^K_S0 -> ^pi0 ^pi0]']
toolsBsigCh2 += ['Vertex', '^B0:ch2 -> ^phi ^K_S0']
toolsBsigCh2 += ['MCVertex', '^B0:ch2 -> ^phi ^K_S0']
toolsBsigCh2 += ['TagVertex', '^B0:ch2']
toolsBsigCh2 += ['MCTagVertex', '^B0:ch2']
toolsBsigCh2 += ['DeltaT', '^B0:ch2']
toolsBsigCh2 += ['MCDeltaT', '^B0:ch2']
toolsBsigCh2 += ['CustomFloats[isSignal:isContinuumEvent]', '^B0:ch2']

if action == 'expert':
    toolsBsigCh2 += ['ContinuumSuppression', '^B0:ch2']
    toolsBsigCh2 += ['ContinuumSuppression[FS1]', '^B0:ch2']
    toolsBsigCh2 += ['CustomFloats[' + ':'.join(networkOutput) + ']', '^B0:ch2']
    toolsBsigCh2 += ['CustomFloats[' + ':'.join(transformedNetworkOutputLPCA) + ']', '^B0:ch2']
    #toolsBsigCh2 += ['CustomFloats[' + ':'.join(transformedNetworkOutputNB) + ']', '^B0:ch2']

#toolsBsigCh2 += ['FlavorTagging', '^B0:ch2']

toolsRS = ['RecoStats', '^B0:ch2']


# save stuff to root file
ntupleFile(outFile)

#ntupleTree('Trks', 'pi+:all', toolsTrk)
#ntupleTree('Neu', 'gamma:all', toolsNeu)
#ntupleTree('Pi0', 'pi0:all', toolsPi0)
#ntupleTree('KsNeu', 'K_S0:neu', toolsKsNeu)
#ntupleTree('PhiKK', 'phi:all', toolsPhiKK)
ntupleTree('B0_ch2', 'B0:ch2', toolsBsigCh2)
ntupleTree('RecoStats', 'B0:ch2', toolsRS)

# Process the events
process(analysis_main)

# print out the summary
print statistics
