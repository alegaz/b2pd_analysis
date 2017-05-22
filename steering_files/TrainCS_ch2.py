#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################################
#
# Apply continuum suppression training to:
#
# B0 -> phi[K+K-] K0[pi0pi0]
#
# Author: A. Gaz
#
######################################################

from basf2 import *
from modularAnalysis import *
import basf2_mva
import subprocess

train_data = 'SkimmedCS_ch2_train.root'
test_data = 'SkimmedCS_ch2_test.root'


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

general_options = basf2_mva.GeneralOptions()
general_options.m_datafiles = basf2_mva.vector(train_data)
general_options.m_treename = "B0_ch2"
general_options.m_identifier = "MVAFastBDT_ch2.root"
general_options.m_variables = basf2_mva.vector(*trainVars)
general_options.m_target_variable = "isNotContinuumEvent"
general_options.m_transformToProbability = False

fastbdt_options = basf2_mva.FastBDTOptions()
# more options here? nTrees, etc...


# Train a MVA method and store the weightfile (MVAFastBDT_ch2.root) locally.
basf2_mva.teacher(general_options, fastbdt_options)
