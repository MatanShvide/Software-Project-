# Software-Project-
C and Python interface
Final stage in my software project for the second semester in my second year in Tel Aviv University.

Task:

Execute the Spectral Clustering Kmeans algorithm in a C/Python interface.

1. Read datapoints from file in Python file and send to Cmodule
2. Run Jacobi algorithm in C and return k vectors to Python
3. Choose k rows as initial centroids with k++ algorithm and send indices to C
4. Find k final centroids by kmeans algorithm in C 
5. Return to python and print

*C file can run indpendently
*It is also possible to run the spk algorithm to a certain step of your choice (Weight adjacency matrix, Diagonal, Lnorm, Jacobi, and Complete spk)
