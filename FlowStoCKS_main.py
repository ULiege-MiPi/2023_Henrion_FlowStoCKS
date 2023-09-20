"""
Created on Mon Oct 11 13:56:30 2021

@author: lucas
This code aims at simulating the phenotype distribution within a Bioreactor.
The population is initially split into 10,000 cells that, at each time step, 
grow and consume their substrate according to their phenotype
"""
import numpy as np
import pandas as pd
import random


def FlowStoCKS(BioModel,ProcessParam):
    Bio = BioModel
    Process = ProcessParam
    Nb_cell = Process['Nb_cell']
    len_t=Process['len_t']
    time=Process['time']
    space = int(Nb_cell*1.5)

    # initialize empty vectors and define parameters
    X        = np.zeros((len_t,Nb_cell+space))
    f_u      = np.zeros((len_t,Nb_cell+space))
    delay    = np.zeros((len_t,Nb_cell+space))

    # set every variable to a list of zeros
    S_glu,Pulse_glu,Inducer= (np.zeros((len_t)) for _ in range(3))

    #initial biomass, glucose and arabinose concentration (g/l)
    for i in range(Nb_cell):
        X[0,i]   = abs(np.random.normal(1,1/2))*(Process['InitBiomass']/Nb_cell) #cells start with a biomass following a guassian to avoid intial sync
    S_glu[0]    = Process['InitCglucose']

    # set inducer pulse concentration
    Pulse_glu[:] = Process['Pulse_glu']
    # initialize fluorescence
    f_u[0,:] = np.random.gamma(2,(Bio['Fu_intial']/2),Nb_cell+space)
    #select inducer
    if Bio['Inducer']=='glu':
        Inducer = S_glu
    # compute growth rate given substrate and inhibitor concentration 
    func_µ_glu = lambda S_glu,f_u: (Bio['µ_max_glu']*(S_glu/(Bio['Ks_glu'] + S_glu))
                            *(Bio['KI']/(Bio['KI'] + f_u)))

    #set t_step, it is always the same
    t_step = time[1]-time[0]
    #set probability of beeing flushed
    PFlushed = Process['D']*t_step
    #initialize growth correction and step (the correction aims at reducing the growth rate and re-do the step if the glucose concentration
    # is negative at the end of the simulation step)
    correction = 1
    step = 1
    while step < len_t:

        S_glu[step] = S_glu[step-1]- Process['D']*S_glu[step-1]*t_step + Process['D']*Process['Feed_glucose']*t_step + Pulse_glu[step]

        #compute flushing of cell by dilution
        CellOut = PFlushed>np.random.rand(np.sum(X[step-1,:]!=0))
        listCell = np.array(np.where(X[step-1,:]!=0)[0])  
        X[step-1,listCell[CellOut]] = 0 
        delay[step-1,listCell[CellOut]]=0
        f_u[step-1,listCell[CellOut]]=0

        #set probability of switching to alternative phenotype
        Pinduction   = ( (Inducer[step-1]**Bio['n'])/
                             ( (Bio['k']**Bio['n']) + (Inducer[step-1]**Bio['n']) ) )

        #re-set list of cells and no cells
        listCell = list(np.where(X[step-1,:]!=0)[0])
        listNoCell = list(np.where(X[step-1,:]==0)[0])
        randomlist = random.sample(listCell, len(listCell))
        #initialize index of possible first new cell
        cell_junior = 0

        for cell in range(0,len(randomlist)):
            
            µ_glu =  func_µ_glu(S_glu[step-1],f_u[step-1,randomlist[cell]])

            µ_glu = correction*µ_glu  
            
            # Biomass
            X[step,randomlist[cell]] = ( X[step-1,randomlist[cell]]
                                       + (µ_glu*X[step-1,randomlist[cell]]*t_step))
            
            # glucose concentration
            S_glu[step]  = ( S_glu[step]
                           - µ_glu*X[step-1,randomlist[cell]]*t_step/Bio['y_glu'])
               
            # make cell play against the probability of switching given by Pinduction, if a cell win, it accumuylates tstep it its memory
            if Pinduction > random.random():

                delay[step,randomlist[cell]] = delay[step-1,randomlist[cell]] + t_step
                
            else:
                delay[step,randomlist[cell]] = delay[step-1,randomlist[cell]]
            #once a cell has accumulated enought tstep, it switches and produces fluorescence units (gfp concentration goes up)
            if delay[step,randomlist[cell]]>=Bio['delay_time']:

                f_u[step,randomlist[cell]] = ( f_u[step-1,randomlist[cell]]
                                + np.random.gamma(2,Bio['GfpProd']/2)*t_step
                                - (µ_glu*f_u[step-1,randomlist[cell]]*t_step) )

                delay[step,randomlist[cell]] = 0

            else:

                f_u[step,randomlist[cell]] = (f_u[step-1,randomlist[cell]])
    
            #if biomass of a cell double, the cell splits it with all the phenotype into two cells
            if X[step,randomlist[cell]] >= (2*(Process['InitBiomass']/Nb_cell)):
                X[step,randomlist[cell]] = X[step,randomlist[cell]]/2
                delay[step,randomlist[cell]] = delay[step,randomlist[cell]]/2
                f_u[step,randomlist[cell]] = f_u[step,randomlist[cell]]/2

                X[step,listNoCell[cell_junior]] = X[step,randomlist[cell]]
                delay[step,listNoCell[cell_junior]] = delay[step,randomlist[cell]]
                f_u[step,listNoCell[cell_junior]] = f_u[step,randomlist[cell]]
                cell_junior+=1
            #fluorescence cannot go bellow autofluorescence
            if  f_u[step,randomlist[cell]]<np.random.gamma(2,100):
                f_u[step,randomlist[cell]]=np.random.gamma(2,100)

        #If glucose concentration at the end of time step is negative, re-do simulation of time step but with a lower growth rate (correction factor)
        if S_glu[step]<0:
            Available = (S_glu[step-1]- Process['D']*S_glu[step-1]*t_step + Process['D']*Process['Feed_glucose']*t_step - Pulse_glu[step])
            correction = (Available/(Available-S_glu[step]))
            step -= 1 
        else:
            correction = 1   

        step += 1

    return time,f_u,S_glu,X,delay
