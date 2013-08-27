# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 21:25:19 2013

@author: xinghualu
"""

import pickle

#from SigNetNode import SigNetNode
from PySBVSigNetClass import PySBVSigNet


edgeFile = "Data/Reference_Network.mody20130826.csv"
#edgeFile = "Data/Reference_Network.csv"
nodeFile = "Data/Refnet.node.type.csv"

net = PySBVSigNet(nodeFile, edgeFile)

net.assocData("Data/Rat.data.matrix.csv")

trainedNetwork = net.gibbsUpdate(pickleDumpFile="tmpHumanAugmentedSBVSigNet-alpha-07.pickle", alpha = .7)

pickle.dump(trainedNetwork, open("convergedHumanAugmentedSBVSigNet-alpha-07.pickle", 'wb'))