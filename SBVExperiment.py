# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 21:25:19 2013

@author: xinghualu
"""

#import pickle

#from SigNetNode import SigNetNode
from PySBVSigNetClass import PySBVSigNet
import numpy as np

edgeFile = "Data/Reference_Network.augmented-08-29-13.csv" 
#edgeFile = "Data/Reference_Network.csv"
nodeFile = "Data/Refnet.node.type.csv"

net = PySBVSigNet(nodeFile, edgeFile)

net.assocData("Data/Human.matrix.csv")
#net.assocData("Data/Rat.data.matrix.csv")
alphas = np.array(range(1, 9)) * .1
for a in alphas:
    net.gibbsUpdate(pickleDumpFile="pickles/bestHumanAugmentedNet-09-01-13.alpha-"+ str(a) +".pickle", alpha = a, nChains = 20)


