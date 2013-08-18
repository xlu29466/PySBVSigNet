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
from rpy2 import robjects 
import rpy2.robjects.numpy2ri
from StringIO import StringIO
from SigNetNode import SigNetNode  # import the node class


rpy2.robjects.numpy2ri.activate()   # enable directly pass numpy arrary or matrix as arguments to rpy object
R = robjects.r                      # load R instance
R.library("glmnet")
glmnet = R('glmnet')                # make glmnet from R a callable python object
lm_fit = R('lm.fit')                # make "lm.fit" from R a callable Pyhton object

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
        
        # initial look up table of parameters and mapping between network node 
        # and datamatrix column indices
        self.dictNodeParams = dict()
        self.dictNode2MatrixIndx = dict()
        self.dictParentOfNodeToMatrixIndx = dict()  # a diction of list indices
       
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
        print "Added " + str(len(edgeLines)) + " edges.  Done with creating network"


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
        
        colnames.pop(0)  # get rid of the name of first column        
        
        nodesWoData = set(self.network.nodes()) - set(colnames)
        if len(nodesWoData) > 0:
            for node in nodesWoData:
                self.network.remove_node(node)
                
        # read in data and generate a numpy data matrix
        self.data = np.genfromtxt(StringIO(lines), delimiter = ",", usecols=tuple(range(1, len(colnames) + 1)))
        print "Data matrix dimension " + str(np.shape(self.data))
            
        #check in which column the data for a node in graph locates and populate dictionaries
        for node in self.network:
            try:
                nodeIndex = colnames.index(node)  
            except ValueError:
                raise ValueError("The data for node " + node + " is not in data Matrix.  Quit!")
            self.dictNode2MatrixIndx[node] = nodeIndex            
        
            # find column indices for the predecessors 
            preds = self.network.predecessors(node)
            if len(preds) > 0:
                try: 
                    self.dictParentOfNodeToMatrixIndx[node] = [(colnames.index(p) ) for p in preds]
                except:
                    raise ValueError("The data for certain node is not in data matrix.  Quit!")
            else:
                self.dictParentOfNodeToMatrixIndx[node] = []
               
        print "Done with associating data to network"
        #print str(self.dictNode2MatrixIndx)
        #print str(self.dictParentOfNodeToMatrixIndx)
                
        
    def randomInitParams(self):
        """Initialize the parameter vector associated with each node 
           with a random vector sampled from standard Gaussian.
           
           The parameters associated with each node for each MCMC chain is 
           stored in two dimensional numpy array with nChain * nPredecessors 
        """
        print "Initialize parameters associated with each node in each MCMC chain"
        for nodeId in self.network:           
            parents = self.network.predecessors(nodeId)
            if len(parents) > 0:
                self.dictNodeParams[nodeId] = np.random.randn(self.nChains, len(parents) + 1)
            else:
                self.dictNodeParams[nodeId]  = None
            
            # To do, change intial values for edges with high confidence
            
   
    
    def gibbsUpdate(self, nChains = 10, nSamples = 10, maxIter = 5000, p = 0.2):
        """ Sampling the states of hidden variables using Gibbs sampling.
            
            Each node take binary state.
            Conditional distribution of each node p(x | Pa(x)) is a logistic 
            function of its parents.  Update of each node is conditioning on 
            its Markov blanket.            
        """
        nCases, nVariables = np.shape(self.data)
        self.nChains = nChains
        self.randomInitParams()
        print "Start Gibbs sampling update."
        
        # set up Markov chains. 
        self.nodeStates = list()
        self.expectedStates = list()
        evidenceNodes = [n  for n in self.network.nodes() if self.network.node[n]['nodeObj']] 
        hiddenNodes = list(set(self.network.nodes()) - set(evidenceNodes))
        for c in range(self.nChains):  
            # each Markov chain keeps a state matrix
            self.nodeStates.append(self.data)
            # random intialize the hidden states for each chain
            for nodeId in hiddenNodes:
                nodeIndx = self.dictNode2MatrixIndx[nodeId]
                self.nodeStates[c][:, nodeIndx] = 0
                rn = np.random.rand(nCases)
                self.nodeStates[c][rn <= p,nodeIndx] = 1
                
            # each chain collect expected statistics of nodes from samples along the chain
            self.expectedStates.append(np.zeros(np.shape(self.data)))

                        
        bConverged = False
        nIter = 0
        burnIter = 300
        sampleCount = 0
        while not bConverged:
            nIter +=1
            if nIter > maxIter:
                break

            # E-step of EM
            self._updateStates()            
            if  nIter % 2 == 0:
                sampleCount += 1
                for c in range(self.nChains):
                    self.expectedStates[c] = self.expectedStates[c] + self.nodeStates[c]
                
                
            # M-step of EM
            if sampleCount >= nSamples:
                sampleCount = 0
                 # take expectation of sample states
                self.expectedStates = map(lambda x: x / nSamples, self.expectedStates)
                self._updteParams()
                print "Iteration: " + str(nIter) + "; log marginal probability of observed variables: " + str(self.calcEvidenceMarginal())
                for c in range(self.nChains):
                    self.expectedStates[c] = np.zeros(np.shape(self.data))
            
            bConverged = self._checkConvergence()
             
            
            
    def _checkConvergence(self):
        # To do, add convergence checking code
        return False
                        

    def _updateStates(self):
        nCases, nVariables = np.shape(self.data)
        # interate through all nodes. 
        for c in range(self.nChains):
            for nodeId in self.network:
                curNodeIndx = self.dictNode2MatrixIndx[nodeId]
                # skip observed nodes
                if self.network.node[nodeId]['nodeObj'].bMeasured:
                    continue
                
                curNodeMarginal = self._calcNodeMarginal(nodeId, c)
                
                # sample states of current node based on the prob, and update 
                sampleState = np.zeros(nCases)
                sampleState[curNodeMarginal >= np.random.rand(nCases)] = 1.
                self.nodeStates[c][:, curNodeIndx] = sampleState


    def _calcNodeMarginal(self, nodeId, c):
        """
        Calculate the marginal probability of a node's state set to "1" 
        
        args:
             nodeId   A string id of the node of interest
             c        An integer indicate the chain from which the parameter 
                         vector to be used  
        """
        nCases, nVariables = np.shape(self.data)
        # collect the state of the predecessors of the node
        predIndices = self.dictParentOfNodeToMatrixIndx[nodeId]  
        logProbOneCondOnParents = 0
        logProbZeroCondOnParents = 0
        if len(predIndices) > 0:  # if the node has parents  
        # calculate p(node = 1 | parents);   
            nodeParams = self.dictNodeParams[nodeId][c,:] 
            predStates =  np.column_stack((np.ones(nCases), self.nodeStates[c][:, predIndices])) 
            pOneCondOnParents = 1 / (1 + np.exp( - np.dot(predStates, nodeParams)))  
            logProbOneCondOnParents  = np.log(pOneCondOnParents)
            logProbZeroCondOnParents = np.log(1 - pOneCondOnParents)

        # collect  evidence from all children 
        logProbDChildCondOne = 0
        logProdOfChildCondZeros = 0
        children = self.network.successors(nodeId)
        curNodeIndx = self.dictNode2MatrixIndx[nodeId]

        if len(children) > 0:
            for child in children:   
                # collect data and parameters associated with the node
                childNodeParams = self.dictNodeParams[child][c,:] 
                curChildStates = self.nodeStates[c][:, self.dictNode2MatrixIndx[child] ]                    
                        
                # find out predecessors and the position of curNode in the list
                childPredIndices = self.dictParentOfNodeToMatrixIndx[child] 
                curNodePosInPredList = childPredIndices.index(curNodeIndx)  
        
                # Collect states of the predecessors of the child
                # Set the state of current node to ones 
                childPredStatesWithCurSetOne = self.nodeStates[c][:, childPredIndices]    
                childPredStatesWithCurSetOne[:, curNodePosInPredList] = np.ones(nCases)  
                childPredStatesWithCurSetOne = np.column_stack((np.ones(nCases), childPredStatesWithCurSetOne)) # padding data with a column ones as bias
                pChildCondCurNodeOnes = 1 / (1 + np.exp(-np.dot(childPredStatesWithCurSetOne, childNodeParams)))
                logProbDChildCondOne += np.log ( curChildStates * pChildCondCurNodeOnes + (1 - curChildStates) * (1 - pChildCondCurNodeOnes))
        
                # set the state of the current node (nodeId) to zeros 
                childPredStatesWithCurSetZero = self.nodeStates[c][:,childPredIndices]
                childPredStatesWithCurSetZero [:, curNodePosInPredList] = np.zeros(nCases)
                childPredStatesWithCurSetZero = np.column_stack((np.ones(nCases), childPredStatesWithCurSetZero))
                pChildCondCurNodeZeros = 1 / (1 + np.exp(- np.dot(childPredStatesWithCurSetZero, childNodeParams))) 
                logProdOfChildCondZeros += np.log(curChildStates * pChildCondCurNodeZeros + (1 - curChildStates) * (1 - pChildCondCurNodeZeros))

        # now we can calculate the marginal probability of current node 
        curNodeMarginal = 1 / (1 + np.exp(logProbZeroCondOnParents + logProdOfChildCondZeros - logProbOneCondOnParents - logProbDChildCondOne))
        return curNodeMarginal
       
        
        
    def _updteParams(self):
        # Update the parameter associated with each node, p(n | Pa(n)) using logistic regression,
        # using expected states of precessors as X and current node states acrss samples as y
        nCases, nVariables = np.shape(self.data)
        for c in range(self.nChains):
            for nodeId in self.network:
                print "Estimate parameter for node" + nodeId
                predIndices = self.dictParentOfNodeToMatrixIndx[nodeId]
                nodeIdx = self.dictNode2MatrixIndx[nodeId]
                
                if len(predIndices) > 0: 
                    x = np.column_stack((np.ones(nCases), self.expectedStates[c][:, predIndices]))  # create data matrix containing data of predecessor nodes
                    y = self.nodeStates[c][:, nodeIdx]
                    #print str(x)
                    #print str(y)
                
                    # call logistic regression funciton from Rpy
                    fit = lm_fit (x, y, family = "binomial")
                    # extract coefficients from Rpy2 vector object
                    self.dictNodeParams[nodeId][c,:] = np.array(fit[0])  
                    print "params: " + str(self.dictNodeParams[nodeId][c,:])
                
                
    def calcEvidenceMarginal(self):
        # collect all evidence nodes
        evidenceNodes = [n  for n in self.network.nodes() if self.network.node[n]['nodeObj']]  
        # calculate the marginal by calcualte probability of instantiation of a
        # node in all samples and average chains and samples 
        logTotalMarginal = 0
        for c in range(self.nChains):
            for node in evidenceNodes:
                curNodeIndx = self.dictNode2MatrixIndx[node]
                curNodeStates = self.nodeStates[c][:,curNodeIndx]
                nodeProb = self._calcNodeMarginal(node, c)
                logProb = np.log(curNodeStates * nodeProb + (1 - curNodeStates) * (1 - nodeProb) )
                logTotalMarginal += sum(logProb)
                
        return logTotalMarginal / c 
        
        
    def trimEdgeWithLasso(self):
        pass
    
            
