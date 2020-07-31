#!/usr/bin/python
#########################################
# File Name: parseC2_and_simul_max_onlycoding.py
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

# usage ./parseC2_and_simul_max_onlycoding.py res.ana4pops_allchr_collapsedC2pergene.bed > Max_C2_per_gene_analysis_with_codingSNP_1000000permut.txt
#
###########################################

import sys
import numpy as np

# here we need to choose in which distribution we are going to sample.
# We can either keep only SNP in coding regions
# Or use genomic SNP but outside repeats regions (masked)
# Let's do genomic SNP first

C2genefile = sys.argv[1]



# First add all C2values in a list
C2values = []
with open(C2genefile,'r') as f:
    line = f.readline()
    for line in f:
        c = line.split()
        if c[5]!=".":
            j=map(float,c[5].split(','))
            C2values.extend(j)

arr = np.array(C2values)
import random
import bisect
random.sample(arr,3)
numsim=1000000
thres = int(0.95*numsim)
print >> sys.stderr, "Read all C2 values, now going into file again..."

count = 0
#with open(C2genefile,'r') as f:
with open(sys.argv[1],'r') as f:
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
        C2maxobs = max(j)
        simmax=[]
        for sim in range(numsim):
            tmp = random.sample(arr,loc)
            simmax.append(max(tmp))
        smax = sorted(simmax)[thres]
        pvalmax = bisect.bisect_left(sorted(simmax),C2maxobs)
        pvalmax = pvalmax/float(numsim)
        count+=1
        print >> sys.stderr, str(count)+" genes analyzed"
        #if C2maxobs > smax:
        #    outlier = True
        #print name,ID,chrom,stop-start,loc,C2meanobs,smean,C2maxobs,smax,outlier
        if C2maxobs > smax:
            outlier = True
        print name,ID,chrom,stop-start,loc,C2maxobs,smax,pvalmax



