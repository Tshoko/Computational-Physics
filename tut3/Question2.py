# -*- coding: utf-8 -*-
"""
Created on Mon May  4 01:46:09 2020
code for the plots https://notes.quantecon.org/submission/5b3b102eb9eab00015b89f8e
@author: Boitshoko
"""

# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import scipy.optimize as opt
import json
import intergration

#Constants
hist_data=np.array([3,4,3,8,8,2,5,0,1,2])
print(hist_data)


#Define the narmal pdf
def gaussian (x, mu, sigma):
    
    return (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2) )
# Define function that generates values of a normal pdf
def norm_cdf(x, mu, std, S):
    
    #S is the upper bound sum in counts for the distribution

    #cummulative pdf 
    Cum_prob = stats.norm.cdf(S, loc=mu, scale=std)
            
    pdf=gaussian (x, mu, std)/Cum_prob
    
    return pdf

def Likelihood(x, mu,std,S,B):
    """Calculates likelihood for data giventhe 
    given mean(mu) and standard deviation(std).
    S is the upper bound sum in counts for the distribution
    B Backround
    """ 
    x=x-B
    
    vals=norm_cdf(x, mu, std, S)
    
    #return the likelihood
    return [np.sum(np.log(vals)),np.product(vals)]
def Min(params, *args):
    """
     Use to minimize the log_likelihood from given parameters 
     by negating it
    """
    mu, std ,S,B= params
    x= args
    log_lik_val = Likelihood(x, mu,std,S,B)[0]

    return -log_lik_val

results = intergration.Json()
results["2"]={"a":Likelihood(hist_data, 5, 2, 20,1)[1]}
with open("tut3/MTSBOI022/tut3.json","w") as f:
        json.dump(results,f)

print('likelihood: ', Likelihood(hist_data, 5, 2, 20,1)[1])

params_init = np.array([5, 2,20,1])
results = opt.minimize(Min, params_init, hist_data)
m,st,S,B=results.x
print(results)


print(results)

x=[]
for i in range(10):x.append(i+0.5)
#plt.plot(x,norm_cdf(x, 3, 2, 20))
plt.bar(x,hist_data,width=1,edgecolor="black")
x=np.linspace(0,10,100)
plt.plot(x,norm_cdf(x, m, st, S)*40,color="r")
plt.show()
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
cmap1 = matplotlib.cm.get_cmap('summer')

mu_vals = np.linspace(0,10, 50)
sig_vals = np.linspace(1, 8, 50)
lnlik_vals = np.zeros((50, 50))
for mu_ind in range(50):
    for sig_ind in range(50):
        lnlik_vals[mu_ind, sig_ind] = Likelihood(hist_data, mu_vals[mu_ind],
                                                   sig_vals[sig_ind], 20,0)[1]
mu_mesh, sig_mesh = np.meshgrid(mu_vals, sig_vals)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(sig_mesh, mu_mesh, lnlik_vals, rstride=2,
                cstride=2, cmap=cmap1)
ax.set_title('likelihood for values of mu and sigma')
ax.set_xlabel(r'$\sigma$')
ax.set_ylabel(r'$\mu$')
ax.set_zlabel(r'likelihood')
'''






print(results)

x=[]
for i in range(10):x.append(i+0.5)
#plt.plot(x,norm_cdf(x, 3, 2, 20))
plt.bar(x,hist_data,width=1,edgecolor="black")
x=np.linspace(1,10,10)
plt.plot(x,norm_cdf(x, m, st, S)*40+B)
'''