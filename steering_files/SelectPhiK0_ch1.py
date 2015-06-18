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
from FlavorTagger import *

# standard lists of particles
from stdFSParticles import stdPi0
from stdFSParticles import stdPhoton
from stdLooseFSParticles import stdLoosePi
from stdLooseFSParticles import stdVeryLoosePi
from stdLooseFSParticles import stdLooseK
from stdV0s import stdKshorts


# set the input files
import sys
filelistSIG = ['../saved_rootfiles/B0_phi-KK_KS-pi+pi-_gsim-BKGx1-50000-*.root']

# load input ROOT file
inputMdstList(filelistSIG)

# print contents of the DataStore before loading Particles
printDataStore()

# use the default functions for all other particles
stdPhoton()
stdPi0()
stdLoosePi()
stdVeryLoosePi()
stdLooseK()
stdKshorts()

# print contents of the DataStore after loading Particles
printDataStore()

# start reconstructing the candidates
# two lists of Kshorts can be defined (those from V0's are recommended)

# Ks from V0's
fillParticleList('K_S0:all', '0.3 < M < 0.7', True)
vertexKFit('K_S0:all', 0.0)

# Ks -> pi+pi-
reconstructDecay('K_S0:chg -> pi+:loose pi-:loose', '0.4 < M < 0.6')
vertexKFit('K_S0:chg', 0.0)

# phi -> K+K-
reconstructDecay('phi:KK -> K+:loose K-:loose', '0.1 < M < 4.1')
vertexKFit('phi:KK', 0.0)


# ch1: B -> phi[KK] Ks[pi+pi-]
#reconstructDecay('B0:ch1 -> phi:KK K_S0:chg', '5.0 < M < 5.5')
reconstructDecay('B0:ch1 -> phi:KK K_S0:all', '4.5 < M < 5.5')
vertexRave('B0:ch1', 0.0, 'B0:ch1 -> [phi -> ^K+ ^K-] K_S0')

# get the rest of the event:
buildRestOfEvent('B0:ch1')

# perform MC truth matching
matchMCTruth('B0:ch1')

# get tag vertex ('breco' is the type of MC association)
TagV('B0:ch1', 'breco')

# get continuum suppression (needed for flavor tagging)
#buildContinuumSuppression('B0:ch1')

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
#toolsBsigCh1 += ['FlavorTagging', '^B0:ch1']

toolsRS = ['RecoStats', '^B0:ch1']


# save stuff to root file
ntupleFile('B0_Phi-KK_K0-pi+pi-.root')

ntupleTree('Trks', 'pi+:all', toolsTrk)
ntupleTree('Neu', 'gamma:all', toolsNeu)
ntupleTree('Pi0', 'pi0:all', toolsPi0)
ntupleTree('KsChg', 'K_S0:chg', toolsKsChg)
ntupleTree('PhiKK', 'phi:KK', toolsPhiKK)
ntupleTree('B0_ch1', 'B0:ch1', toolsBsigCh1)
ntupleTree('RecoStats', 'B0:ch1', toolsRS)

# Process the events
process(analysis_main)

# print out the summary
print statistics
