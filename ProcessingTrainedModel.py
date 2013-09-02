1# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 08:29:35 2013

@author: xinghua
"""

import networkx as nx

# this script process the trained PySBVSigNet models and identify a best model
# for rat and human respectively.  
import os, re, cPickle
from PySBVSigNetClass import PySBVSigNet
from SigNetNode import SigNetNode

pickledir = "/home/xinghua/projects/src/PySBVSigNet/pickles/"
os.chdir(pickledir)
files = os.listdir(pickledir)
models = list()
bestLikelihood = float("-inf")
bestModelFileName = ""
for f in files:
    #if re.search('bestHumanAugmentedNet-alpha-03.pickle$', f):
    if re.search('bestRatAugmented.+?\.pickle$', f):
        print "Loading " + f
        m = cPickle.load(open(f, 'rb'))
        likelihood = m.getLikelihoodArray()[-1]
        if  likelihood > bestLikelihood:
            bestLikelihood = likelihood
            bestModel = m
            bestModelFileName = f
        
print  pickledir + bestModelFileName + ".  \nAlpha: " + str(bestModel.alpha) + " Likelihood: " + str(bestLikelihood)
print "The number of nodes : " + str(len(bestModel.network.nodes())) + "; the number of edges: " + str(len(bestModel.network.edges()))
print "Start trimming edges"
#
#bestModel.trimNetwork()
bestModel.trimNetworkByConsensus(0.95)
print "After trimming, the graph has " + str(len(bestModel.network.edges()))
outf = open (bestModelFileName+"network_edges.csv", 'w')
seenNodes = set()
for e in bestModel.network.edges():
    seenNodes.add(e[0])
    seenNodes.add(e[1])
    print str(e)
    outf.write  (e[0] + "," + e[1] + "\n")

for n in bestModel.network:
    if n not in seenNodes:
        outf.write(n + "\n")
        
outf.close()

#nx.write_graphml(bestModel.network, pickledir + "human.best.model.xgmml")


        
