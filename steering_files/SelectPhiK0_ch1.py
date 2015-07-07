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

from basf2 import *
from modularAnalysis import *


# set the input files
filelistSIG = ['../../../../rel-00-05-00/analysis/b2pd_analysis/saved_rootfiles/B0_Phi-KK_Ks-pi+pi-_gsim-BKGx1-50000-1.root']

# load input ROOT file
inputMdstList(filelistSIG)

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
vertexRave('B0:ch1', 0.0, 'B0:ch1 -> [phi -> ^K+ ^K-] K_S0')
matchMCTruth('B0:ch1')

# get the rest of the event:
buildRestOfEvent('B0:ch1')

# get tag vertex ('breco' is the type of MC association)
TagV('B0:ch1', 'breco')

# get continuum suppression (needed for flavor tagging)
buildContinuumSuppression('B0:ch1')

# flavor tagging
#FlavorTagger(
#    mode='Expert',
#    weightFiles='B2JpsiKs_mu',
#    categories=['Electron', 'Muon', 'KinLepton', 'Kaon', 'SlowPion', 'FastPion', 'Lambda', 'FSC', 'MaximumP*', 'KaonPion'])

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

toolsKsChg = ['InvMass', '^K_S0:all -> pi+ pi-']
toolsKsChg += ['Vertex', '^K_S0:all -> pi+ pi-']
toolsKsChg += ['EventMetaData', '^K_S0']

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
toolsBsigCh1 += ['InvMass', '^B0:ch1 -> [^phi -> K+ K-] ^K_S0:all']
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
toolsBsigCh1 += ['CustomFloats[isSignal]', '^B0:ch1']
#toolsBsigCh1 += ['FlavorTagging', '^B0:ch1']

toolsRS = ['RecoStats', '^B0:ch1']


# save stuff to root file
ntupleFile('reco_test2.root')

ntupleTree('Trks', 'pi+:all', toolsTrk)
ntupleTree('Neu', 'gamma:all', toolsNeu)
ntupleTree('Pi0', 'pi0:all', toolsPi0)
ntupleTree('KsChg', 'K_S0:mdst', toolsKsChg)
ntupleTree('PhiKK', 'phi:all', toolsPhiKK)
ntupleTree('B0_ch1', 'B0:ch1', toolsBsigCh1)
ntupleTree('RecoStats', 'B0:ch1', toolsRS)

# Process the events
process(analysis_main)

# print out the summary
print statistics
