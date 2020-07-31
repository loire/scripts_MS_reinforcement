#!/usr/bin/python
#########################################
# File Name: parseC2_and_simul.py
#
# Purpose:
#
# Usage:
#
# Creation date: 31-10-2017
#
# Last modified:
#
# Created by:
#
###########################################

import sys
import numpy as np

C2file = sys.argv[1]
C2genefile = sys.argv[2]

outfile = open(C2genefile+"_mean_C2_simul_results_nsimtot10000.txt","w")

# commande line
# parseC2_and_simul_group_Cmean.py res.ana4pops_allchr_C2.bed res.ana4pops_allchr_collapsedC2pergene.bed



## First add all C2values in a list
#C2values = []
#with open(C2genefile,'r') as f:
#    line = f.readline()
#    for line in f:
#        c = line.split()
#        if c[5]!=".":
#            j=map(float,c[5].split(','))
#            C2values.extend(j)
#

## Store all values with genomic coordinates in order.
## Store in dic for each chromosome

listchrom = []
for i in xrange(1,20):
    listchrom.append(str(i))

listchrom.append("X")
dicochrom = {k:[] for k in listchrom}

#C2genefile = "res.ana4pops_allchr_C2.bed"

with open(C2file,'r') as f:
    line = f.readline()
    for line in f:
        c = line.split()
        chrom = c[0]
        pos = int(c[1])
        C2 = float(c[3])
        dicochrom[chrom].append((C2,pos))

import random

# Function to pick a random chrom weighted by its size
def select_weighted(d):
   offset = random.randint(0, sum(d.itervalues())-1)
   for k, v in d.iteritems():
      if offset < v:
         return k
      offset -= v

dicosize = {}

for i in dicochrom:
    dicosize[i]=len(dicochrom[i])

#select_weighted(dicosize)


# To extract mean of a group of SNP in the all genome
# We want the same number of SNP as in the gene
# and roughly the same length
# num = number of  SNP in gene
# leng = length of gene
# chrom = chromosome name (a string like "1" or "X")
# stdlength variance around the length
## exemple of values:
#stdlength = 0.1
#chrom = "1"
#num= 10
#leng = 10000
# Returned values:
# Val = position in the chrom distribution
# moy = mean C2 values from val position to val + num position if leng is ok
def getC2s(num,leng,stdlength,chrom):
    ratio = 0
    cpt = 0
    minlengthratio = 1-stdlength
    maxlengthratio = 1+stdlength
    while ratio > maxlengthratio or ratio < minlengthratio:
        val = random.randint(0,len(dicochrom[select_weighted(dicosize)]))
        if val > len(dicochrom[chrom])-num:
                     continue
        sample = dicochrom[chrom][val:val+num]
        taille = sample[-1][1] - sample[0][1]
        ratio = taille/float(leng)
        cpt +=1
    #print val,taille,ratio
    #print "Essai : ",cpt
    moy = 0
    for i in sample:
        moy+=i[0]
    moy/=len(sample)
    #print moy
    return((moy,val))


#   exemple of paramters values
#Zval = 0.99 # number of values before threshold in sim.
#stdlength = 0.1
#chrom = "1"
#num= 10
#leng = 500
## nsimtot: number of desired values for the permut
#nsimtot = 1000
#maxtry = 10 # max number of trials (times simtots before giving away)
#
import bisect


def get_gene_pvals(chrom,num,leng,stdlength,nsimtot,maxtry,Zval,C2meanobs):
    sim_ok={} # dic to store mean values for position (keys are positions)
    ntry=0
    nsim=0
    while nsim < nsimtot and ntry< maxtry*nsimtot:
        ntry+=1
        res = getC2s(num, leng, stdlength , chrom)
        if not res[1] in sim_ok.keys():
            sim_ok[res[1]]=res[0]
            nsim+=1
#            print "new val in results, ",nsim," done so far"
    print "Number of trials: ",ntry
    #print sim_ok
    if ntry >= maxtry*nsimtot:
        return(("NA","NA"))
    else:
        th =int(Zval * nsimtot)
        pvalmax = bisect.bisect_left(sorted(sim_ok.values()),C2meanobs)
        threshold = sorted(sim_ok.values())[th]
        return((1-(pvalmax/float(nsimtot)),threshold))



#get_gene_pvals("X",2,500,0.1,10000,5,0.99,3.26)


#get_gene_pvals("X",14,2056,0.1,10000,5,0.99,3.26)



count = 0
with open(C2genefile,'r') as f:
    for line in f:
        c=line.split()
        if c[5]==".":
            continue
        chrom = c[0]
        start = int(c[1])
        stop = int(c[2])
        name = c[3]
        ID = c[4]
        j=map(float,c[5].split(','))
        loc = len(j)
        if loc>1:
            C2meanobs = np.mean(np.array(j))
            print "gene to analyze:"
            print name,ID,chrom,stop-start,loc,C2meanobs
            ressim = get_gene_pvals(chrom,loc,stop-start,0.2,10000,100,0.99,C2meanobs)
            count+=1
            print str(count)+" genes analyzed"
        #if C2maxobs > smax:
        #    outlier = True
        #print name,ID,chrom,stop-start,loc,C2meanobs,smean,C2maxobs,smax,outlier
            outfile.write("\t".join(map(str,[name,ID,chrom,stop-start,loc,C2meanobs,ressim[1],ressim[0]]))+"\n")
        else:
            outfile.write("\t".join(map(str,[name,ID,chrom,stop-start,loc,C2meanobs,"NA","NA"]))+"\n")











