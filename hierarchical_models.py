#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 09:56:59 2018

@author: cfaust
"""

##
import pandas
import matplotlib.pyplot as plt
from scipy.stats import beta
## for mcmc package uses python script
import sampyl as smp
from sampyl import np

#scenario1
path1="/data/bayesian/malaria_clean_maintext.csv"
#scenario2; 3-unknowns as non vir or vir
#path1="/data/bayesian/malaria_clean_unknowns.csv"
#include all data
#path1="/data/bayesian/malaria_clean_wcountries.csv"
path2="/data/bayesian/phenotype_clean.csv"
#path2="/data/bayesian/phenotype_clean_vietnamascam.csv"
#path2="/data/bayesian/phenotype_clean_X.csv"
mal = pandas.read_csv(path1)
phenotypes = pandas.read_csv(path2)

data = pandas.merge(mal, phenotypes, on='region')

vir = np.array(data['virmalinf'])
## uncomment here to include unknowns
#vir= np.array(data['vir_unknown'])
mal_N = np.array(data['mal_N'])
pheno_N = np.array(data['animals_N'])
not_A = np.array(data['not_A'])
   
nonvirp = []
nonvirp2 = []

cutoff =0.2
for i in range(len(mal_N)):
    virulent = vir[i]
    N = mal_N[i]
    nonvirp.append(beta.cdf(cutoff,1+virulent,1+N-virulent))
    nonvirp2.append(beta.mean(1+virulent,1+N-virulent))



## Here we are defining a function to estimate two parameters; 
## theta 1 (non-vir) & theta 2 (vir)
def logp(t1, t2):
    model = smp.Model()
    print(t1,t2)
    #retVal = 0.0
    arr = np.arange(len(mal_N))
    for j in range(len(mal_N)):
        i = int(arr[j])
        nA = not_A[i]
        N = pheno_N[i]
        s = nonvirp[i]
        ## log-likelihood
        model.add(np.log((s*(t1**nA)*((1.0-t1)**(N-nA)))  + ((1.0-s)*((t2**nA)*((1.0-t2)**(N-nA)))))) 
    ## log-priors
    
    model.add(smp.uniform(t1,lower=0,upper=1),
              smp.uniform(t2, lower=0, upper=1))

    return model()   

start = smp.find_MAP(logp, {'t1': 0.5, 't2': 0.5},bounds={'t1':(0.0, 1.0), 't2':(0.0,1.0)})
sampler = smp.Metropolis(logp, start)
chain = sampler.sample(1000, burn=200, thin = 100)
plt.plot(chain.t1,color='g')
plt.plot(chain.t2,color='b')
plt.show()

import matplotlib.pyplot as plt
plt.hist(chain.t1, bins=100,color='g',range=[0,1])
plt.hist(chain.t2, bins=100,color='b',range=[0,1])
plt.show()  
plt.plot(vir/mal_N.astype(float),not_A/pheno_N.astype(float), '.')
plt.xlim(0,0.1)
start
