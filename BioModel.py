#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 13:46:49 2022

@author: lucas
µ_max_glu, Ks_glu, y_glu: Monod growth parameters on glucose
n, k: Repons function parameters
delay_time: delay in h between cellular switching decision and actual switching (GFP prod)
GfpProd: Production of GFP per switching event
KI: Inhibition constant if apply
Fu_intial: Mean fluorescence value at the start of the simularion
Inducer: define the inducer type, alternative carbon source (c2) or glucose
"""

PNoBurden = {'µ_max_glu': 0.54, 'Ks_glu': 0.034, 'y_glu':0.5,
    'n':-2,'k':0.05 ,'Pbasal':0,'delay_time':0.4,'GfpProd':1e5,'KI':None,'Fu_intial':0,'Inducer':'glu'}


PLowBurden = {'µ_max_glu': 0.54, 'Ks_glu': 0.034, 'y_glu':0.5,
    'n':-2,'k':0.05 ,'Pbasal':0,'delay_time':0.4,'GfpProd':1e5,'KI':1e5/20,'Fu_intial':0,'Inducer':'glu'}


PHighBurden = {'µ_max_glu': 0.54, 'Ks_glu': 0.034, 'y_glu':0.5,
    'n':-2,'k':0.05 ,'Pbasal':0,'delay_time':0.4,'GfpProd':1e5,'KI':1e5/200,'Fu_intial':0,'Inducer':'glu'}

