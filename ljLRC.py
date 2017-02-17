import os
import sys
import numpy as np 
import pandas as pandas

def readLJfromPRM(filename):
    infi=open(filename).readlines()
    NBS=0
    for n in range(len(infi)):
	if 'NONBONDED' in infi[n]:
		NBS=n+2
		break
    LJpar={}
    atno = 0
    for line in infi[NBS:]: 
	LJpar[atno]={'AT':line.split()[0],'EPS':-1.0*float(line.split()[2]),'SIG':float(line.split()[3])/0.561231}
        atno+=1 
    return LJpar 

def CalcLRC(solvpar,solupar,Nsolv,cutoff,denSolv):
    LRC = 0 
    for i in  solvpar.keys():
	for j in solupar.keys():
	    sigma = np.sqrt(solvpar[i]['SIG']*solupar[j]['SIG'])
	    epsil = np.sqrt(solvpar[i]['EPS']*solupar[j]['EPS'])
            ratio = sigma/cutoff
            print i,j, (pow(ratio,9)/3.0 - pow(ratio,3))*epsil*pow(sigma,3)
    	    LRC = LRC + (pow(ratio,9)/3.0 - pow(ratio,3))*epsil*pow(sigma,3)
    print LRC
    LRC = LRC*8*np.pi*Nsolv*denSolv/3.0
    return LRC

solv=readLJfromPRM('CYH.prm')
solu=readLJfromPRM('MET.prm')
print CalcLRC(solv,solu,Nsolv=140,cutoff=10.0,denSolv=0.774)
