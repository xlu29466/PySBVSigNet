""" PySBVSigNet

This class represent the data from SBV Improver challenge as a Bayesian network 
and infer the network structure.

The signaling network consists of 3 types of node: stimuli, proteins and genes 
expression.  Gene expression data, stimuli phosphorylation state of some 
proteins are given as observed variables.  The task is to infer the activation 
state of unobserved proteins and infer the network among proteins and genes.

Starting from the provided reference network, we will create a cannonical 
netwrok for both rat and human. This implementation apply Gibbs sampling to 
instanatiate the activation state of proteins.  For proteins with observed
phosphorylation state, we assume phosphorylated protein is in an active state, 
although it impact on downstream moleclule can be activation or inhibition.   
After instantiate latent variables, the network structure contains software 
for modeling causal relationship between proetin kinases and their target 
proteins.

"""

import networkx as nx
import numpy as np
from numpy import matlib
from rpy2 import robjects 
import rpy2.robjects.numpy2ri
from StringIO import StringIO
from SigNetNode import SigNetNode  # import the node class


rpy2.robjects.numpy2ri.activate()   # enable directly pass numpy arrary or matrix as arguments to rpy object
R = robjects.r                      # load R instance
R.library("glmnet")
glmnet = R('glmnet')                # make glmnet a callable python object
lm_fit = R('lm.fit')                # make "lm.fit" a callable Pyhton object

class PySBVSigNet:
    ## Constructor.  Create an empty instanc.
    def __init__(self, nodeFileName = None, edgeFileName = None, dataMatrix = None):
        """Constructor
        
           Args:
              nodeFileName  A string of pathname of the file containing the nodes
                            of the network
              edgeFileName  A string of pathname of the file containing the edges
                            of the network
        """
        self.network = None
        self.data = None
        
        if nodeFileName and edgeFileName:
            self.initNetwork(nodeFileName, edgeFileName)
            
        if self.network and dataMatrix:
            self.assocData(dataMatrix)

    def getNetwork(self):
        return self.network
        

    def initNetwork(self, nodeFileName, edgeFileName):
        """  Create an instance of network using networkx by parsing definition 
             files of nodes and edges
        
            Args:
              nodeFileName   CSV text file; each line describe one node 
                             (gene, protein, or stimuli, whether measured )
              edgeFileName   CSV text file; each line provide source and sink 
                             node of an edge
        """  
        
        if self.network:
            self.network.clear()

        try:
            nf = open(nodeFileName, "r")
            nodeLines = nf.readlines()
            nf.close()
        except IOError:
            print "Failed to open the file containing nodes"
            return
            
        if len(nodeLines) == 1:  # Mac files end a line with \r instead of \n
            nodeLines = nodeLines[0].split("\r")
           
        try:
            ef = open(edgeFileName, 'r')
            edgeLines = ef.readlines()
            ef.close()
        except IOError:
            print "Fail to open file " + nodeFileName
            return
        if len(edgeLines) == 1:
            edgeLines = edgeLines[0].split("\r")  

        # creat a networkx graph 
        self.network = nx.DiGraph()
        print "Creating network"

        # parse nodes
        for line in nodeLines:
            #print line
            name, nodeType, bMeasured = line.split(',')
            if bMeasured == '1':
                bMeasured = True
            elif bMeasured == '0':
                bMeasured = False
            else:
                print "unknown type for bMeasured variable"
                continue
            # create a node and add to network
            self.network.add_node(name, nodeObj=SigNetNode(name, nodeType, bMeasured))
        print "Added " + str(len(nodeLines)) + " nodes"

        # parse edges
        for line in edgeLines:
            #print "Current edge line: " + line
            source, sink = line.split(',')
            if source not in self.network or sink not in self.network:
                print "Cannot add edges between nodes not in the network"
                continue
            self.network.add_edge(source, sink)
        print "Added " + str(len(edgeLines)) + " edges"
            
        # randomly initialize the parameters of the network    
        self.randomInitParams()
        print "Done with creating network"


    def assocData(self, dataFileName):
        """  Associate data with the network of an instance.
        
            Args:
                dataFileName  Pathname of text file containing data matrix.  
                              Rows are cases; columns are values of variables

            The data file is an nCase * nNodes matrix, with latent nodes 
            randomly instantiated to 0/1 values, whereas observed
            variables are set to the observed values
        """
       
        try:
            f = open(dataFileName)
            lines = f.readlines()
        except IOError:
            print "Fail to read in file " + dataFileName
        
        colnames = None
        if len(lines) == 1:  # Mac version csv, with "\r" as return
            lines = lines[0].split("\r")
            colnames = lines.pop(0).split(',') # split  header and extract colnames
            map(lambda x: x.rstrip(), lines)  # remove the "\r"
            lines = "\n".join(lines)
        else:
            colnames = lines.pop(0).split(',')
            lines = "".join(lines)
            
        # read in data and generate a numpy data matrix
        self.data = np.genfromtxt(StringIO(lines), delimiter = ",", usecols=tuple(range(1, len(colnames)))) 
            
        #check in which column the data for a node in graph locates
        self.dictNode2MatrixIndx = dict()  
        for node in self.network:
            nodeIndex = colnames.index(node)
            if not nodeIndex:  # there is no data for the node
                raise Exception("The data for node " + node + " is missing.  Quit!")
            self.dictNode2MatrixIndx[node] = nodeIndex
        
        # find column indices for the predecessors of a node
        self.dictParentOfNodeToMatrixIndx = dict()  # a diction of list
        for node in self.network:
            preds = self.network.predecessors(node)
            if len(preds) > 0:
                self.dictParentOfNodeToMatrixIndx[node] = []          

                for p in preds:                
                    self.dictParentOfNodeToMatrixIndx[node].append(self.dictNode2MatrixIndx[p])  
                
        print "Done with associating data to network"
                
        
    def randomInitParams(self):
        """Initialize the parameter vector associated with each node 
           with a random vector sampled from standard Gaussian
        """
        print "Initialize parameters associated with each node"
        for nodeId in self.network:           
            parents = self.network.predecessors(nodeId)
            if len(parents) > 0:
                self.network.node[nodeId]['nodeObj'].params = np.random.randn(len(parents) + 1)
            else:
                self.network.node[nodeId]['nodeObj'].params = None
            
            # To do, change intial values for edges with high confidence
            
    def initEdgeParam(self, source, sink, priormu):
        """To initialize the parameter of associated with a specific edge
           based on prior knowledge.
           
           Args:
             source  The source node of the edge
             sink    The sink node of  the edge
             priormu The mean of a Guassiun distribution from which to sample a weight
        """
        pass
    
    
    def gibbsUpdate(self, nChains = 1, nSamples = 10):
        """ Sampling the states of hidden variables using Gibbs sampling.
            
            Each node take binary state.
            Conditional distribution of each node p(x | Pa(x)) is a logistic 
            function of its parents.  Update of each node is conditioning on 
            its Markov blanket.            
        """
        # set up Markov chains.  We replcate rows of data matrix "nChain" times
        # to create multiple chains 
        self.nodeStates = matlib.repmat(self.data, nChains, 1)
        
        # create a matrix to collect mulitple samples from a single chain
        # to collect expected states of nodes within a chain
        self.expectedStates = np.zeros(np.shape(self.nodeStates))

        bConverged = False
        nIter = 0
        maxIter = 5000
        burnIter = 100
        sampleCount = 0
        while not bConverged:
            nIter +=1
            if nIter > maxIter:
                break

            # E-step of EM
            self._updateStates()            
            if nIter > burnIter and nIter % 5 == 0:
                self.expectedStates = self.expectedStates + self.nodeStates
                sampleCount += 1
                
            # M-step of EM
            if sampleCount >= nSamples:
                sampleCount = 0
                 # take expectation of sample states
                self.expectedStates = self.expectedStates / nSamples
                self._updateParams()
                self.expectedStates = np.zeros(np.shape(self.nodeStates))
            
            bConverged = self._checkConvergence()
            
    def _checkConvergence(self):
        # To do, add convergence checking code
        return False
                        

    def _updateStates(self):
        #Performing Gibbs sampling of the node states

        # interate through all nodes.  Potentially parallelizable
        for nodeId in self.network:
            print "Gibbs updating node: " + nodeId + "; has " + str(len(self.network.predecessors(nodeId))) + " predecessors"
            
            nodeObj = self.network.node[nodeId]['nodeObj'] 
            nodeParams = nodeObj.params
            
            print "param vector: " + str(nodeObj.params)
            # skip observed node
            if nodeObj.bMeasured:
                continue
            
            nrow, ncol = np.shape(self.nodeStates)
            
            # collect the state of the parent nodes
            curNodeMatrixIndx = self.dictNode2MatrixIndx[nodeId]
            print "Index of current node: " + str(curNodeMatrixIndx)
            parentIndx = self.dictParentOfNodeToMatrixIndx[nodeId]
            print "Indices of predecessors: " + str(parentIndx)
            
            if len(parentIndx) > 0:  # if the node has parents
                parentStates = np.column_stack((np.ones(nrow), self.nodeStates[:,parentIndx]))  
                print str(parentStates)

                # calculate p(node = 1 | parents)
                pOneCondOnParents = 1 / (1 + np.exp( - np.dot(parentStates, nodeParams)))  
                logProbCondOnParents = np.log(pOneCondOnParents)
                logProbZeroCondOnParents = np.log(np.ones(nrow) - pOneCondOnParents)
            else:
                logProbCondOnParents = np.zeros(nrow)
                logProbZeroCondOnParents = np.zeros(nrow)

            # collect all children terms
            logProbDChildCondOne = 0
            logProdOfChildCondZeros = 0
            children = self.network.successors(nodeId)
            if len(children) > 0:
                for child in children:   
                    # collect data and parameters associated with the node
                    nodeParams = self.network.node[child]['nodeObj'].params 
                    curChildStates = self.nodeStates[:, self.dictNode2MatrixIndx[child]]                    
                    childParentIndx = self.dictParentOfNodeToMatrixIndx[child]
                    
                    # find out which column the data from the current node is in
                    curNodePositionInMatrix = childParentIndx.index(curNodeMatrixIndx)  
    
                    # Collect states of the predecessors of the child, condition,
                    # Set the state of current node to ones 
                    childParentStates = self.nodeStates[:, childParentIndx]
                    childParentStates[:, curNodePositionInMatrix] = np.ones(nrow)  # set value of current node to 1s
                    childParentStates = np.column_stack((np.ones(nrow), childParentStates)) # padding data with a column ones as bias
                    pChildCondCurNodeOnes = 1 / (1 + np.exp(-np.dot(childParentStates, nodeParams)))
                    logProbDChildCondOne += np.log ( curChildStates * pChildCondCurNodeOnes + (1 - curChildStates) * (1 - pChildCondCurNodeOnes))
    
                    # set the state of the current node (nodeId) to zeros and calculate the probability of observed values of the child
                    childParentStates [:, curNodePositionInMatrix] = np.zeros(nrow)
                    childParentStates = np.column_stack((np.ones(nrow), childParentStates))
                    pChildCondCurNodeZeros = 1 / (1 + np.exp(- np.dot(childParentStates, nodeParams))) 
                    logProdOfChildCondZeros += np.log(curChildStates * pChildCondCurNodeZeros + (np.ones(nrow) - curChildStates) * (np.ones(nrow) - pChildCondCurNodeZeros))

            # now we can calculate the marginal probability of current node by collecting statistics from Markov blanket
            curNodeMarginal = 1 / (1 + np.exp(logProbZeroCondOnParents + logProdOfChildCondZeros - logProbCondOnParents - logProbDChildCondOne))

            # generate samples based on the prob
            rProb = np.random.rand(nrow)
            sampleState = np.zeros(nrow)
            sampleState[curNodeMarginal >= rProb] = 1.
            self.nodeStates[:, curNodeMatrixIndx] = sampleState


    def _updteParams(self):
        # Update the parameter associated with each node, p(n | Pa(n)) using logistic regression,
        # using expected states of precessors as X and current node states acrss samples as y

        for nodeId in self.network:
            nrow, ncol = np.shape(self.nodeStates)
            ancestors = self.network.predecessors(nodeId)
            
            x = np.column_stack((np.ones(nrow), self.expectedStates[:, ancestors]))  # create data matrix containing data of predecessor nodes
            y = self.nodeStates[:, self.dictNode2MatrixIndx[nodeId]]

            # call logistic regression funciton from Rpy
            fit = lm_fit (x, y, family = "binomial")
            # extract coefficients from Rpy2 vector object
            self.network.node[nodeId]['nodeObj'].params = np.array(fit[0])   
                                
 
 