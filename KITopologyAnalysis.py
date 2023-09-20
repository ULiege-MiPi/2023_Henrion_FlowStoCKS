"""
Created on Mon Oct 11 13:56:30 2021

@author: lucas and JA


"""
import numpy as np
from matplotlib import pyplot as plt
from BioModel import *
from ProcessParam import *
from FlowStoCKS_main import FlowStoCKS
import matplotlib.colors  #if scatterplot failes define bins as float
import matplotlib.cm as cm
import itertools


PBiomass = {'µ_max_glu': 0.54, 'Ks_glu': 0.034, 'y_glu':0.5,
    'n':-2,'k':0.05 ,'Pbasal':0,'delay_time':0.4,'GfpProd':1e5,'KI':0,'Fu_intial':0,'Inducer':'glu'}

#KI=[0,1e-10,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e10,1e99]
KI=[8e3,9e3,1e4,5e4,1e5,1e10,1e99]
#KI=[1e2,5e2,1e3,2e3,3e3,4e3,5e3,6e3,7e3,8e3,9e3,1e4,5e4]
#KI=[1e2,2e2,3e2,4e2,5e2,6e2,7e2,8e2,9e2]

for k in range(9):
    PBiomass["KI"]=KI[k]
    print(PBiomass["KI"])
    sim=FlowStoCKS(PBiomass,Process) #StoCKs takes biological and process values as input

    time,f_u,S_glu,X,delay = sim[0],sim[1],sim[2],sim[3],sim[4]

    # perform stats on the simulation results
    Biomass = np.sum(X,axis=1)
    mean_fl1 = np.zeros((len(f_u)))

    for i in range(len(f_u)):
        mean_fl1[i] = np.mean(f_u[i,np.where(f_u[i,:]!=0)[0]])

    median_fl1 = np.median(np.where(f_u[:-1,:]!=0),axis=1)

    lent = 100
    fu = []
    for i in range(len(f_u)):
        fu.append(f_u[i,np.where(f_u[i,:]!=0)[0]])

    t=np.arange(0,lent+60/3600,60/3600)

    timeSim=[]
    for i in range(len(t)):
        for j in range(len(fu[i])):
            timeSim.append(t[i])

    logbins = np.logspace(0,6,100)


    fluorescenceSim = list(itertools.chain.from_iterable(fu))

    normalisation=matplotlib.colors.LogNorm()
    c=cm.ScalarMappable(cmap='viridis')
    plt.figure(figsize=(10,6))
    plt.hexbin(timeSim,fluorescenceSim,gridsize=500, mincnt=1, edgecolors="none", cmap=c.get_cmap(),yscale='log',norm=normalisation)
    plt.xlabel('Time (h)',fontsize=18)
    plt.ylabel('FL1-A',fontsize=18)
    plt.savefig("/home/lucas/Desktop/PhD/Python/StockTime adder/Results/KI"+str(KI[k])+"_KI.png",dpi=500)

    f_uWithout0 = pd.DataFrame(fu)
    f_uWithout0.to_csv("/home/lucas/Desktop/PhD/Python/StockTime adder/Results/KI"+str(KI[k])+"_uWithout0.gz", compression='gzip') 

    X = pd.DataFrame(X)
    X.to_csv("/home/lucas/Desktop/PhD/Python/StockTime adder/Results/KI"+str(KI[k])+"_x.gz", compression='gzip') 
    

