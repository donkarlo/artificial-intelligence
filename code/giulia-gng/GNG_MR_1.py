# GROWING NEURAL GAS ##########################################################


import matplotlib.pyplot as plt
import numpy as np
import sklearn.metrics as metrics
from sklearn.manifold import TSNE

import time
st=time.time()


# LOADING THE DATA TO CLUSTER #################################################

class ThreeDPosVelFile:
    def __init__(self, filePath: str):
        self.__filePath = filePath

    def getNpArr(self):
        file = open(self.__filePath, "r")
        arr = np.empty((0, 6), np.dtype('float64'))
        for lineCounter, curLine in enumerate(file):
            arr = np.append(arr, np.array([[r for r in curLine.split(",")]]), axis=0)
        opt = arr.astype(np.float)
        return np.asarray(opt)


t3dposVel = ThreeDPosVelFile("pos-vel-obs-from-gps.txt")
inputData = t3dposVel.getNpArr()

# Positional data
# inputData = 'inputData.mat'
#
# inputData = mat4py.loadmat(inputData)
# inputData = [ v for v in inputData.values() ]
# inputData = inputData[0]
#
# inputData = np.asarray(inputData)
# inputData = np.asarray(inputData)


# PARAMETERS DEFINITION #######################################################

N = 100;  # Number of nodes
MaxIt = 100;  # Iteration (repetition of input data)
L_growing = 1000;  # Growing rate
epsilon_b = 0.05;  # Movement of winner node
epsilon_n = 0.0006;  # Movement of all other nodes except winner
alpha = 0.5;  # Decaying global error and mUtility
delta = 0.9995;  # Decaying local error and mUtility
T = 100;  # It could be a function of params.L_growing, e.g., params.LDecay = 2*params.L_growing
L_decay = 1000;  # Decay rate should be faster than the growing then it will remove extra nodes
alpha_utility = 0.0005;
k = 1.5;
seedvector = 2;

PlotFlag = True

# CLUSTERING PROCEDURE ########################################################

# Normalization of the input data
minDataNorm = np.min(inputData, axis=0)
dataNorm = inputData - np.tile(minDataNorm, (inputData.shape[0], 1))  # tile(a, (m, n)) = repmat(a, m, n)
maxDataNorm = np.max(dataNorm, axis=0)
CV = dataNorm / np.tile(maxDataNorm, (inputData.shape[0], 1))
inputNormOrd = CV

# Size of input data (number of training samples)
nData = CV.shape[0]
# Dimension of input data
nDim = CV.shape[1]

# Permutation of the rows of the input data vector
permutationValues = np.random.RandomState(seed=seedvector).permutation(nData)
CV = CV[permutationValues, :]

CVmin = np.min(CV, axis=0)
CVmax = np.max(CV, axis=0)

# Initialization ##############################################################
# Initial 2 nodes for training the algorithm
Ni = 2

wNorm = np.zeros((Ni, nDim))
# It returns an array of random numbers generated from the continuous uniform 
# distributions with lower and upper endpoints specified by 'CVmin' and 'CVmax'.
for i in range(Ni):
    wNorm[i, :] = np.random.uniform(low=CVmin, high=CVmax, size=None)

# error
E = np.zeros(Ni)
# mUtility
utility = np.ones(Ni)
# Connections between nodes
C = np.zeros((Ni, Ni))
# Ages of the edges
t = np.zeros((Ni, Ni))

# Loop ########################################################################

nx = 0

# Loop over the number of iterations parameter
for it in range(MaxIt):

    print("Iteration " + str(it) + " out of " + str(MaxIt))
    print("Number of nodes in current iteration: " + str(wNorm.shape[0]))

    # Loop over the number of data
    for c in range(nData):

        # Select Input

        # Counter of cycles inside the algorithm
        nx = nx + 1
        # pick first input vector from permuted inputs
        x = CV[c, :]

        # COMPETITION AND RANKING
        # pairwise distance between normalized input value and the normalized node means
        # X = np.tile(x, (wNorm.shape[0],1))
        X = np.expand_dims(x, axis=0)
        X_state = X[:, 0: int(nDim / 2)]
        X_deriv = X[:, int(nDim / 2):nDim]
        wNorm_state = wNorm[:, 0: int(nDim / 2)]
        wNorm_deriv = wNorm[:, int(nDim / 2):nDim]
        d_state = metrics.pairwise_distances(X=X_state, Y=wNorm_state, metric='euclidean')[0][:]
        d_deriv = metrics.pairwise_distances(X=X_deriv, Y=wNorm_deriv, metric='euclidean')[0][:]

        value_state = d_state / np.sum(d_state)
        value_deriv = d_deriv / np.sum(d_deriv)

        d = value_state + value_deriv

        # Organize distances between nodes and the first data point in an ascending order
        SortOrder = np.argsort(d)

        # Closest node index to the first data point
        s1 = SortOrder[0]
        # Second closest node index to the first data point                                              
        s2 = SortOrder[1]

        # AGING
        # Increment the age of all edges emanating from s1
        t[s1, :] = t[s1, :] + 1
        t[:, s1] = t[:, s1] + 1

        # Add Error
        dist0 = np.power(d[s1], 2)
        dist1 = np.power(d[s2], 2)
        E[s1] = E[s1] + dist0

        # Utility
        # Initial mUtility is zero in first case and dist is the error of first node
        deltaUtility = dist1 - dist0
        # Difference between error of two nodes
        utility[s1] = utility[s1] + deltaUtility

        # ADAPTATION
        # Move the nearest distance node and it's neibors towards the input signal 
        # by fractions Eb and En resp.
        # 1) move nearest node
        wNorm[s1, :] = wNorm[s1, :] + epsilon_b * (x - wNorm[s1, :])

        # Take all the connections of the closest node to the data in question
        Ns1 = np.where(C[s1, :] == 1)
        # 2) move neighbors
        # for j in Ns1[0]:
        #    wNorm[j,:] = wNorm[j,:] + epsilon_n*(x-wNorm[j,:])   
        wNorm[Ns1, :] = wNorm[Ns1, :] + epsilon_n * (x - wNorm[Ns1, :])

        # Create link
        # If s1 and s2 are connected by an edge , set the age of this edge to 
        # zero , it such edge doesn't exist create it
        C[s1, s2] = 1
        C[s2, s1] = 1
        # Age of the edge
        t[s1, s2] = 0
        t[s2, s1] = 0

        # Remove old links
        # remove edges with an age larger than Amax(a threshold value)
        C[t > T] = 0

        # Number of connections of each node                                                
        nNeighbor = np.sum(C, axis=1)
        # Eliminate alone nodes from the C and t matrix and 
        # from wNorm, E and mUtility vector
        # AloneNodes = (nNeighbor==0)
        indexAloneNodes = np.where(nNeighbor == 0)
        C = np.delete(C, (indexAloneNodes), axis=0)
        C = np.delete(C, (indexAloneNodes), axis=1)
        t = np.delete(t, (indexAloneNodes), axis=0)
        t = np.delete(t, (indexAloneNodes), axis=1)

        wNorm = np.delete(wNorm, (indexAloneNodes), axis=0)

        E = np.delete(E, (indexAloneNodes), axis=0)
        utility = np.delete(utility, (indexAloneNodes), axis=0)

        # ADD NEW NODES at every L_growing

        if ((np.remainder(nx, L_growing) == 0) and (wNorm.shape[0] < N)):
            # Determine the unit q with the maximum accumulated error
            q = np.argmax(E)
            # Maximum index related to the error related to a connected node
            f = np.argmax(C[:, q] * E)

            # Total number of nodes
            r = wNorm.shape[0] + 1
            index_r = r - 1  # to index the new node
            # Insert a new unit r halfway between q and it's neibor f with 
            # the largest error variable
            newNode = (wNorm[q, :] + wNorm[f, :]) / 2

            wNorm = np.vstack([wNorm, newNode])

            # Adding one row and column to C and t
            C_temp = np.zeros((r, r))
            C_temp[0:r - 1, 0: r - 1] = C
            C = C_temp

            t_temp = np.zeros((r, r))
            t_temp[0:r - 1, 0: r - 1] = t
            t = t_temp

            # Remove old connections and introduce the presence of the
            # new created node
            C[q, f] = 0;  # eliminating connections between the two former neighbors
            C[f, q] = 0;
            C[q, index_r] = 1;  # Creating connections between old nodes and new one
            C[index_r, q] = 1;
            C[index_r, f] = 1;
            C[f, index_r] = 1;
            t[index_r, :] = 0;
            t[:, index_r] = 0;

            # Decrease the error variable of q and f by multiplying them with a constant 'alpha'
            E[q] = alpha * E[q]
            E[f] = alpha * E[f]
            # Initialize the error of the new node equal to error of the winner node
            newError = E[q]
            E = np.append(E, newError)
            # Decrease the error variable of q and f by multiplying them with a constand 'alpha'
            # mUtility[q] = alpha*mUtility[q]
            # mUtility[f] = alpha*mUtility[f]
            newUtility = 0.5 * (utility[q] + utility[f])
            utility = np.append(utility, newUtility)

        # REMOVE NODES at every L_decay

        if (np.remainder(nx, L_decay) == 0):

            # Maximum accumulated error
            max_E = np.max(E)
            # Node node_useless having minimum mUtility
            min_utility = np.min(utility)
            node_useless = np.argmin(utility)
            # Utility factor
            CONST = min_utility * k

            if (CONST < max_E):
                # Remove the connection having smaller mUtility factor
                C = np.delete(C, (node_useless), axis=0)
                C = np.delete(C, (node_useless), axis=1)
                # Remove the node having smaller mUtility factor
                wNorm = np.delete(wNorm, (node_useless), axis=0)

                # Remove the min mUtility value from the mUtility vector
                utility = np.delete(utility, (node_useless), axis=0)
                # Remove error vector correspond to the node having min mUtility
                E = np.delete(E, (node_useless), axis=0)
                # Remove aging vector correspond to the node having min mUtility
                t = np.delete(t, (node_useless), axis=0)
                t = np.delete(t, (node_useless), axis=1)

            # E = alpha*E
            # mUtility = alpha*mUtility

        # Decrease errors
        # Decrease error variables by multiplying them with a constant delta
        E = delta * E
        # Decrease the mUtility by alpha_utility constant
        utility = delta * utility

    '''
    dataColorNode = np.zeros(nData)
    for cl in range(nData):
        x = inputNormOrd[cl,:]
        X = np.expand_dims(x, axis=0)
        d = metrics.pairwise_distances(X=X, Y=wNorm, metric='euclidean')[0][:]
        minNode = np.argmin(d)
        dataColorNode[cl] = minNode
        
    wReal = np.zeros((wNorm.shape[0], wNorm.shape[1]))
    countPerCluster = np.zeros(wNorm.shape[0])
    sumPerCluster = np.zeros((wNorm.shape[0], wNorm.shape[1]))
        
    for i in range (nData):
        currentState = CV[c, :]
        currentCluster = dataColorNode[i]
        sumPerCluster[int(currentCluster), :] += currentState
        countPerCluster[int(currentCluster)] += 1
            
    for i in range(wNorm.shape[0]):
        wReal[i, :] = sumPerCluster[i, :]/countPerCluster[i]
        
    # if there were no elements belonging to a particular node, we keep it 
    # as it is, but the mUtility of it is cut to 0
    if countPerCluster[i] == 0:
        wReal[i, :] = wNorm[i, :]
        mUtility[i] = 0
            
    wNorm = wReal
    print('wReal')
    print(wNorm)
        
    if PlotFlag == True: # to ADD
        
        print("plotting...")
    '''

# Clusters of input
dataColorNode = np.zeros(nData)
for c in range(nData):
    x = inputNormOrd[c, :]
    X = np.expand_dims(x, axis=0)
    d = metrics.pairwise_distances(X=X, Y=wNorm, metric='euclidean')[0][:]
    minNode = np.argmin(d)

    dataColorNode[c] = minNode


def reduce_dimensionality(var, perplexity=10):
    dim = var.shape[-1]

    if (dim > 2):
        tsne = TSNE(perplexity=perplexity, n_components=2, init='pca', n_iter=1000)
        var_2d = tsne.fit_transform(var)
    else:
        var_2d = np.asarray(var)

    return var_2d


# PLOTTING


# Positional case
colors = {0: 'black', 1: 'grey', 2: 'blue', 3: 'cyan', 4: 'lime', 5: 'green', 6: 'yellow', 7: 'gold', 8: 'red',
          9: 'maroon'}

plt.ion()
plt.show()

f, axarr = plt.subplots(1, 1, figsize=(10, 10))
var_2d = inputNormOrd[:, 0:2]

for number, color in colors.items():
    axarr.scatter(x=var_2d[dataColorNode == number, 0], y=var_2d[dataColorNode == number, 1], color=color,
                  label=str(number))
    axarr.legend()

axarr.grid()
plt.draw()
# plt.pause(0.002)
f.suptitle("Plotting clustering", fontsize=20)

print("----%.2f----"%(time.time()-st))