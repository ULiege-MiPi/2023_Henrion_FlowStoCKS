#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 16:13:39 2022

@author: lucas
"""

import numpy as np
import pandas as pd

PATHGlC3Seg = '...S.cerevisiae glc3/Segregostat/S.cerevisiae_Segrego_glc3_1/continuous/Reactor1/Summary.csv'

# set time array
SimDuration = 100
time=np.arange(0,SimDuration+60/3600,60/3600)
len_t=len(time)

#initialize empty vector
Feed_c2 = np.zeros((len_t))
Pulse_c2 = np.zeros((len_t))

#extract pulses time an get closest match with simulation time array
def PulseInducer(PATH,CPulse):
    ExpData = pd.read_csv(PATH,sep=';')
    Tpulse_data = np.array(ExpData['Time'][np.where(ExpData['Regulation'] != 'None')[0]])
    for i in range(len(Tpulse_data)):
        difference_array = np.absolute(time-Tpulse_data[i])
        index = difference_array.argmin()
        Pulse_c2[index] = CPulse
    return Pulse_c2

#Define a Process dictionary with all the process values
Process = {'Feed_glucose':5,'D':0.1,'InitBiomass':2.25,'InitC2':0,
           'InitCglucose':0.01,'PulseC2':0,'Nb_cell':10000,
           'time':time,'len_t':len_t,'Feed_c2':0,'Pulse_glu':0}

'''
List of all process conditions considered

#Process Chemostat
Process = {'Feed_glucose':5,'D':0.1,'InitBiomass':2.25,'InitC2':0,
           'InitCglucose':0.01,'PulseC2':0,'Nb_cell':10000,
           'time':time,'len_t':len_t,'Feed_c2':0,'Pulse_glu':0}

#Process Segregostat
Process = {'Feed_glucose':5,'D':0.1,'InitBiomass':2.25,'InitC2':0,
           'InitCglucose':0.01,'PulseC2':0,'Nb_cell':10000,
           'time':time,'len_t':len_t,'Feed_c2':0,'Pulse_glu':PulseInducer(PATHGlC3Seg,0.2)}

'''
