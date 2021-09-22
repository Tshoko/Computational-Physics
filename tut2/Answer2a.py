# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 21:11:40 2020

@author: Boitshoko
"""


import json
import transmition_line



with open("input.json","r") as f:
    results = json.load(f)

print(results)
n=500  #the number ofcells in the circuit
#constants

w=628318530         #frequency

Ri=50                #internal resistance

for q,r in results.items():
    #print(q)
    if(q=="l"):
        RL=float(r)
    elif(q=='C'):
        C=float(r)*10e-2
        print(C)
    elif(q=='L'):
        L=float(r)*10e-2
    elif(q=='R'):
        R=float(r)*10e-2
    elif(q=='G'):
        G=float(r)*10e-2
    elif(q=='G'):
        N=int(r)
n=2*N+1
L=((L)*1j*w)   #inductance
R=R            #resistance
G=G              #ohmic losses
C=(1/(1j*w*C))     #capacitance
R_L=L+R             #Series impedence
RL=RL              #Terminal resistor
V_I=100
GC=1/(1/G+1/C)       #parallel impedance for C and G
Matrix=np.zeros((n,n), dtype=complex) #initialize the N*N MATRIX THAT will


def make(Matrix,w):
    #Make the matrix needed for calculations
    
    for i in range(1,n):
        
       if i%2==0 and i<n-1:
           Matrix[i][i]=L+R 
           Matrix[i][i+1]= 1/(1/G+1/C)
           Matrix[i+1][i+1]=-1
           Matrix[i+1][i+2]=-1
           Matrix[i][i-1]=-1/(1/G+1/C)
           Matrix[i+1][i]=1
    Matrix[0][0]=L+R+Ri
    Matrix[0][1]=1/(1/G+1/C)
    Matrix[1][2]=-1
    Matrix[n-1][n-2]=-1/(1/G+1/C)
    Matrix[n-1][n-1]=RL
    Matrix[1][0]=1
    Matrix[1][1]=-1
    return Matrix
Matrix=make(Matrix,w)

with open("solutions.json","r") as f:
    results = json.load(f)

#edite the solutions jason file using different values for RL
Matrix[n-1][n-1]=0
voltages=Calculate(Matrix)[0]
N=len(voltages)-1
#
results['a']['open']['middle']=voltages[int(N/2)]
results['a']['open']['end']=voltages[N]
Matrix[n-1][n-1]=10e100
voltages=Calculate(Matrix)[0]
results['a']['short']['middle']=voltages[int(N/2)]
results['a']['short']['end']=voltages[N]
Matrix[n-1][n-1]=100
voltages=Calculate(Matrix)[0]
results['a']['100ohm']['middle']=voltages[int(N/2)]
results['a']['100ohm']['end']=voltages[N]
results['c']=RLmin(Matrix)
results['f']=np.absolute(500/((3*10e8)*(L*C)**(0.5)))

print(results)
#write results to a Json file
with open("tut1.json","w") as f:
    json.dump(results,f)
