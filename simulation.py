"""
Created on Mon Oct 11 13:56:30 2021

@author: lucas


"""
import numpy as np
from matplotlib import pyplot as plt
from BioModel import *
from ProcessParam import *
from FlowStoCKS_main import FlowStoCKS
import matplotlib.colors  
import matplotlib.cm as cm
import itertools



sim=FlowStoCKS(PHighBurden,Process) #FlowStoCKs takes biological and process values as input
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

#plot time scatter plot of population fluorescence
logbins = np.logspace(0,6,100)
fluorescenceSim = list(itertools.chain.from_iterable(fu))
normalisation=matplotlib.colors.LogNorm()
c=cm.ScalarMappable(cmap='viridis')
plt.figure(figsize=(10,6))
plt.hexbin(timeSim,fluorescenceSim,gridsize=500, mincnt=1, edgecolors="none", cmap=c.get_cmap(),yscale='log',norm=normalisation)
plt.xlabel('Time (h)',fontsize=18)
plt.ylabel('FL1-A',fontsize=18)
plt.show()
plt.savefig(".../TimeScatter.png",dpi=500)

#Save single cell biomasses and fluo
f_u = pd.DataFrame(f_u)
f_u.to_csv(".../Fu.csv") 

f_uWithout0 = pd.DataFrame(fu) 
f_uWithout0.to_csv(".../FuWithout0.csv") 

X = pd.DataFrame(X)
X.to_csv(".../Biomass.csv") 


