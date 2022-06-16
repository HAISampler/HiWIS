__author__ = "LDIOP"

import time
import sys
from mpmath import *
from random import random


###########################################################################
#                        Elementary functions                             #
###########################################################################

class wItem():
    def __init__(self, item, wpos, wneg):
        self.item = item
        self.wpos = wpos
        self.wneg = wneg
        self.wt = ""
        


def combin(n, k):
    if k>n: return 0
    if k==0 and n==0: return 1
    if k > n//2:
        k = n-k
    x = 1
    y = 1
    i = n-k+1
    while i <= n:
        x = (x*i)//y
        y += 1
        i += 1
    return x

    
def find(tab,i,j,x):
    m=int((i+j)/2)
    if m==0 or (tab[m-1]<x and x<=tab[m]):
        return m
    if tab[m]<x:
        return find(tab,m+1,j,x)
    return find(tab,i,m,x)


def computeCnk(j, l, tabCnk, M):
    if j+1>= len(tabCnk):
        for n in range(len(tabCnk), j+2):
            tabCnk.append([combin(n, k) for k in range(M+2)])
    return tabCnk[j][l]
    

###########################################################################
#                        Weighting algorithm                              #
###########################################################################

def weightDatasetSum(dataset, delim, M, tabCnk):
    wDatabase = []
    weights = []
    z = 0
    with open(dataset, 'r') as base:
        line=base.readline()
        while line:
            trans = line.replace(" \n","").replace("\n","").split(" ")
            wTrans = []
            i = len(trans)-1
            j=1
            while i >= 0:
                info = trans[i].split(delim)
                if j==1:
                    wTrans.append(wItem(info[0], [float(info[1])], [0]))
                else:
                    wpos, wneg = [], []
                    for l in range(1, min(j,M)+1):
                        if l==1:
                            wpos.append(float(info[1]))
                        else:
                            wpos.append(float(info[1])*computeCnk(j-1, l-1, tabCnk, M) + wTrans[0].wpos[l-2] + wTrans[0].wneg[l-2])
                        if l < len(wTrans[0].wpos)+1:
                            wneg.append(wTrans[0].wpos[l-1] + wTrans[0].wneg[l-1])
                    wTrans = [wItem(info[0], wpos, wneg+[0]*(len(wpos)-len(wneg)))] + wTrans            
                if i==0:
                    wt, x = [], 0
                    for l in range(len(wTrans[0].wpos)):
                        x += (wTrans[0].wpos[l] + wTrans[0].wneg[l])/(l+1)
                        wt.append(x)
                    wTrans[0].wt = wt  
                    z += wt[-1] 
                i -= 1
                j += 1
            wDatabase.append(wTrans)
            weights.append(z)
            line=base.readline()
        base.close()
        del base
    return wDatabase, weights


def weightDatasetProd(dataset, delim, M):
    wDatabase = []
    weights = []
    z = 0
    with open(dataset, 'r') as base:
        line=base.readline()
        while line:
            trans = line.replace(" \n","").replace("\n","").split(" ")
            wTrans = []
            i = len(trans)-1
            j=1
            sumUtil = 0
            while i >= 0:
                info = trans[i].split(delim)
                if j==1:
                    wTrans.append(wItem(info[0], [float(info[1])], [0]))
                else:
                    wpos, wneg = [], []
                    for l in range(1, min(j,M)+1):
                        if l==1:
                            wpos.append(float(info[1]))
                        else:
                            wpos.append(float(info[1])*(wTrans[0].wpos[l-2] + wTrans[0].wneg[l-2]))
                        if l < len(wTrans[0].wpos)+1:
                            wneg.append(wTrans[0].wpos[l-1] + wTrans[0].wneg[l-1])
                    wTrans = [wItem(info[0], wpos, wneg+[0]*(len(wpos)-len(wneg)))] + wTrans   
                sumUtil += float(info[1])         
                if i==0:
                    wt, x = [], 0.
                    for l in range(len(wTrans[0].wpos)):
                        try:
                            x += (wTrans[0].wpos[l] + wTrans[0].wneg[l])/(sumUtil**l) #computeCnk(len(trans), l+1, tabCnk, M)
                        except OverflowError:
                            continue
                        wt.append(x)
                    wTrans[0].wt = wt  
                    z += wt[-1] 
                i -= 1
                j += 1
            wDatabase.append(wTrans)
            weights.append(z)
            line=base.readline()
        base.close()
        del base
    return wDatabase, weights

###########################################################################
#                        Sampling algorithm                               #
###########################################################################

def HAISampler(database, weights):
    xt = random()*weights[-1]
    i = find(weights,0,len(weights),xt)
    trans = database[i]
    xl = random()*trans[0].wt[-1]
    l = find(trans[0].wt,0,len(trans[0].wt),xl)+1
    j = len(trans)
    p=1
    y = 0
    pattern = []
    l1=l
    while l>0:
        d = random()*(y*tabCnk[j][l] + trans[p-1].wpos[l-1]+ trans[p-1].wneg[l-1])
        b = y*tabCnk[j-1][l-1] + trans[p-1].wpos[l-1]
        if d <= b:
            pattern.append(trans[p-1].item)#+":"+trans[p-1].wpos[0])
            y += trans[p-1].wpos[0]
            l -= 1
        p += 1
        j -= 1  
    return str(pattern)


def HPISampler(database, weights):
    xt = random()*weights[-1]
    i = find(weights,0,len(weights),xt)
    trans = database[i]
    xl = random()*trans[0].wt[-1]
    l = find(trans[0].wt,0,len(trans[0].wt),xl)+1
    j = len(trans)
    p=1
    y = 1
    pattern = []
    l1=l
    while l>0:
        d = random()*(y*(trans[p-1].wpos[l-1]+ trans[p-1].wneg[l-1]))
        b = y*trans[p-1].wpos[l-1]
        if d <= b:
            pattern.append(trans[p-1].item)
            y *= trans[p-1].wpos[0]
            l -= 1
        p += 1
        j -= 1  
    return str(pattern)
    

###########################################################################
#                           Record the sample                             #
###########################################################################

def printSample(sampledPatterns, name, opUtil):
    print("*********************** Recording... **************************")
    with open("Output/"+opUtil+"/"+name+"_"+str(M)+".txt", 'w') as output:
        output.write("Pattern\t&\tFrequency in the sample\n")
        for pattern in sampledPatterns:
            output.write(pattern+"\t"+str(sampledPatterns[pattern])+"\n")
    output.close()
    print("********************* End recording! **************************\n")
    
###########################################################################
#                            HiWIS function                               #
###########################################################################

def HiWIS(database, weights, N, opUtil):
    sampledPatterns = {}
    if opUtil=="Sum":
        for k in range(N):
            pattern = HAISampler(database, weights)
            if pattern in sampledPatterns.keys():
                sampledPatterns[pattern] += 1
            else:
                sampledPatterns[pattern] = 1
    elif opUtil=="Prod":
        for k in range(N):
            pattern = HPISampler(database, weights)
            if pattern in sampledPatterns.keys():
                sampledPatterns[pattern] += 1
            else:
                sampledPatterns[pattern] = 1
    return sampledPatterns
    
    
##########################################################################
#                            main function                               #
##########################################################################

if __name__ == '__main__':
    print("#############################################################################")
    print("#         Welcome to the HiWIS tool for high weighted itemset sampling!     #")
    print("#############################################################################\n")
    datasets    =  ["adult"]#, "BMS", "mushroom", "retail", "chess", "foodmart"] 
    delim   =   ":" # Delimits an item with its utility in a transaction
    N    =   10000     # The desired sample size
    M   =   5         # Maximum length constraint
    opUtil = "Prod" # Sum or Prod
    for name in datasets:
        dataset = "Datasets/"+name+".num"
        print("Dataset :",name)
        database, weights =[],[]
        print("M =",str(M))            
        # Preprocessing phase  
        beginTime = time.time()       
        if   opUtil == "Sum":
            tabCnk = []  
            database, weights = weightDatasetSum(dataset, delim, M, tabCnk)
        elif   opUtil == "Prod":
            database, weights = weightDatasetProd(dataset, delim, M)
        endTime = time.time() - beginTime
        print("\tPreprocessing time (s)   :",endTime)
        # Sampling phase 
        beginTime = time.time() 
        sampledPatterns = HiWIS(database, weights, N, opUtil)
        endTime = time.time() - beginTime
        print("\tSampling time (s)        :",endTime)
        print("\tDistinct sampled patterns:",len(sampledPatterns))
        printSample(sampledPatterns, name, opUtil)
        del sampledPatterns
        print("#############################################################################")
    print("#############################################################################")
    print("#                     Thank for using HiWIS with",opUtil+"!                      #")
    print("#############################################################################\n")


