#-------------------------------------------------------------------------------
# Name:        Bayes A/B Test
# Purpose:     Analyze testing continuously
#
# Author:      Zach Riddle
#
# Created:     29/07/2014
#-------------------------------------------------------------------------------

from matplotlib import use
use('wx')
from pylab import *
from scipy.stats import beta, norm, uniform
from random import random
from numpy import *
import numpy as np
import os

##################################################
# Input data
#Default Obs and Test Obs
N = [2492,2463]
#Default conversions and Test Conversions
s = [1415,1407]
###################################################


def prob_winner(y,n,ndraws):
    #number of arms
    k=len(y)
    if k==2:
        ans = prob_B_wins(y[0],n[0]-y[0],y[1],n[1]-y[1])
        ans = [1-ans,ans]
    else:
        #number of failures
        b=[x-z for x,z in zip(n,y)]
        post=np.zeros((ndraws,k))
        sample = np.zeros(k)
        #loop to sample from posterior distributions
        for i in range(0,k):
            post[:,i]=beta.rvs(y[i]+1,b[i]+1,size=ndraws)
        #Count number of wins for each distribution
        for i in range(0,ndraws):
            sample[np.argmax(post[i,:])]+=1

        ans=sample/np.sum(sample)
    return ans

def lbeta(a,b):
    beta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a+b)
    return beta

def prob_B_wins(s0, b0, s1, b1):
    total = 0.0
    for i in range(0,(s1)):
        total += exp(lbeta(s0+i+1, b1+b0+2) - log(b1+i+1) - lbeta(1+i, b1+1) - lbeta(s0+1, b0+1))
    return total

def expected_loss(N,s):
    if len(N)==2:
        s0=s[0]
        s1=s[1]
        N0=N[0]
        N1=N[1]
        if(s[1] / float(N[1])) < (s[0] / float(N[0])):
            s0=s[1]
            s1=s[0]
            N0=N[1]
            N1=N[0]
        prob1 = prob_B_wins(s0,N0-s0,s1+1,N1-s1)
        prob2 = prob_B_wins(s0+1,N0-s0,s1,N1-s1)
        loss = abs((s0+1)/float(N0+2)*(1-prob1) - (s1+1)/float(N1+2)*(1-prob2))/(s0/float(N0))
    else:
        conv = [(x+1)/float(y+2) for x,y in zip(s,N)]
        probs = prob_winner(s,N,500000)
        loss = 0.0
        for i in range(0,len(N)):
            st1 = s
            st1[conv.index(max(conv))]+=1
            probs1=prob_winner(st1,N,500000)
            if(conv.index(max(conv))!=i):
                loss += probs1[i]*conv[i]-(probs1[i]*max(conv)*1.024)
        loss = abs(loss)/min([(x)/float(y) for x,y in zip(s,N)])
    return loss

e_loss = expected_loss(N,s)
probs = prob_winner(s,N,500000)

print 'Expected Loss from pushing out the current winner = '+str(e_loss)
print 'Probability Default will win = ' + str(probs[0])
print 'Probability Test Group will win = ' + str(probs[1]) #+ ',' + str(probs[2])











