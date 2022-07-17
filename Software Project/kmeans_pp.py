from cmath import inf
from time import process_time_ns
import numpy as np
import pandas as pd
import math
import sys
import random
import mykmeanssp as cd1

def kmeansPP(k , max_iter , eps , file1 , file2): 
    #extract and merge data
    data1 = pd.read_csv(file1, sep=",", header=None)
    data2 = pd.read_csv(file2, sep=",", header=None)
    data1 = pd.read_csv(file1, sep=",", header=None).set_index(0)
    data2 = pd.read_csv(file2, sep=",", header=None).set_index(0)
    data_merged = pd.concat([data1 , data2] , axis = 1)
    sdm = data_merged.sort_index()
    N = len(sdm)
    dims = len(sdm.columns)

    #k-means++
    data = sdm.to_numpy()
    sum_dists = 0
    i = 0
    np.random.seed(0)
    new_cent = np.random.randint(0 , N) 
    dists_probs = np.full((2,N) , math.inf)
    centroids = np.zeros((1,dims)) + data[new_cent]
    indices = [0 for i in range(k)]
    indices[0] = new_cent
    while (i < k-1): 
        dists_probs[0] = np.minimum(dists_probs[0] , np.sum(np.power((data - centroids[i]) , 2) , axis = 1))
        sum_dists = np.sum(dists_probs[0] , axis = 0)
        np.true_divide(dists_probs[0] , sum_dists , out = dists_probs[1])
        new_cent = np.random.choice(N , p = dists_probs[1]) 
        i += 1
        indices[i] = (int(new_cent)) #add new centroid
        centroids = np.vstack([centroids , data[new_cent]])
    sdm.to_csv('new_file' , index=False , header=False)
    return indices

#check arguments
try:
    k = int(sys.argv[1])
except:
    print("Invalid Input!")
    sys.exit(1)
if len(sys.argv) != 6: #max_iter not given
    max_iter = 300
    try:
        eps = float(sys.argv[2])
    except:
        print("Invalid Input!")
        sys.exit(1)
    file1 = sys.argv[3]
    file2 = sys.argv[4]
else: #max_iter given
    try:
        eps = float(sys.argv[3])
    except:
        print("Invalid Input!")
        sys.exit(1)
    try:
        max_iter = int(sys.argv[2])
    except:
        print("Invalid Input!")
        sys.exit(1)
    file1 = sys.argv[4]
    file2 = sys.argv[5]

centroids = kmeansPP(k , max_iter , eps , file1 , file2) #get initial k++ centroids
output = cd1.fit(int(k) , int(max_iter) , float(eps) , centroids) #get final centroids from kmeans C

#print output by requested format
centP = [str(a) for a in centroids]
print(','.join(centP))
new_res = [[0 for i in range(len(output[0]))] for j in range(k)]
for i in range(k):
    for j in range(len(output[0])):
        new_res[i][j] = round(output[i][j], 4)
        if j == len(output[0]) - 1:
            print(new_res[i][j])
        else:
            print(new_res[i][j] , end = '')
            print("," , end = '')
print("")    
