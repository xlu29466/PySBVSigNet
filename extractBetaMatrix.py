# this module contain some utility funcitons for working wiht glmnet in RPy2


def extractBetaMatrix(glmnet_res):        
    """ Parse the 'beta' matrix returned by calling glmnet through RPy2.
        Return the first column of 'beta' matrix of the glmnet object 
        with all none zero values returned by the glmnet
        """
                
    # Read in lines of beta matrix txt, which is a nVariables * nLambda.
    # Since we call glmnet by padding x with a column of 1s, we only work
    # with the 'beta' matrix returned by fit
    betaLines = StringIO(str(glmnet_res.rx('beta'))).readlines()
    dimStr = re.search("\d+\s+x\s+\d+", betaLines[1]).group(0)
    if not dimStr:
        raise Exception("'parse_glmnet_res' could not determine the dims of beta")
    nVariables , nLambda = map(int, dimStr.split(' x ')) 
    betaMatrix = np.zeros( (nVariables, nLambda), dtype=np.float)
    
    # glmnet print beta matrix in mulitple blocks with 
    # nVariable * blockSize
    blockSize = len(betaLines[4].split()) - 1
    curBlockColStart = - blockSize
    for line in betaLines:  #read in blocks
        m = re.search('^V\d+', line)
        if not m:  # only find the lines begins with 'V\d'
            continue
        else:
            rowIndx = int(m.group(0)[1:len(m.group(0))]) 
        if rowIndx == 1:
            curBlockColStart += blockSize
            
        # set 'rowIndx' as start from 0
        rowIndx -= 1

        fields = line.rstrip().split()
        fields.pop(0)
        if len(fields) != blockSize:
            blockSize = len(fields)
        for j in range(blockSize):
            if fields[j] == '.':
                continue
            else:
                betaMatrix[rowIndx, curBlockColStart + j] = float(fields[j])

    return betaMatrix
