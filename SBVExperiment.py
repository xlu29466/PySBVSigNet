# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 21:25:19 2013

@author: xinghualu
"""

import pickle

from SigNetNode import SigNetNode
from PySBVSigNetClass import PySBVSigNet


edgeFile = "Data/Reference_Network.csv"
nodeFile = "Data/Refnet.node.type.csv"

net = PySBVSigNet(nodeFile, edgeFile)

net.assocData("Data/Rat.data.matrix.csv")

trainedSigNet = net.gibbsUpdate()

pickle.dump(trainedSigNet, open("trainedSBVSigNet.pickle", 'wb'))