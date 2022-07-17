#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/*declarations*/
double dist(double *, double *);
void print2Darray(double **);
double toDouble(const char* , int , int );
void extract(char *);
static PyObject *kmeansC(int , int , double , PyObject *);
int col = 0;
int row = 0;

static PyObject* fit(PyObject* self , PyObject* args){
    int k;
    int max_iter;
    double eps;
    PyObject *centroids;

    if (!PyArg_ParseTuple(args , "iidO", &k, &max_iter, &eps, &centroids)){
        printf("Invalid Input!");
        return NULL;
    }
    return kmeansC(k , max_iter , eps , centroids); /*call kmeanc algo function*/ 
}

static PyMethodDef myMethods[] = {
    {"fit" , fit, METH_VARARGS, PyDoc_STR("k-means algo")},
    {NULL , NULL , 0 , NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    "Kmeans Module",
    -1,
    myMethods
};

PyMODINIT_FUNC PyInit_mykmeanssp(void){
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m){
        return NULL;
    }
    return m;
}

double toDouble(const char* s, int start, int stop) {
    int m = 1;
    int dot = 1;
    double ret = 0;
    int i;
    for (i = stop; i >= start; i--) {
        if(s[i] == '.'){
            dot = pow(10 , stop-i); /*adjust decimals*/
            continue;
        }
        if(s[i] == '-'){
            ret *= -1;
            break;
        }
        if (s[i] > 57 || s[i] < 48){ /*check validity*/
            printf("Invalid Input!");
            exit(1);
        }
        ret += (s[i] - '0') * m;
        m *= 10;
    }
    ret /= dot;
    return ret;
}

void extract(char *fileName){ /*read file and receive dimensions*/
    char *cont;
    char line[1000]; 
    int currCol = 0;
    int maxCol = 0;
    FILE *myFile1 = NULL;
    myFile1 = fopen(fileName , "r");
    if (myFile1 == NULL){
        printf("Invalid Input!");
        exit(1);
    }
    while (fgets(line , sizeof(line) , myFile1)){
        row++;
        cont = strtok(line , ",");
        if ((maxCol != currCol) && (maxCol != 0)){ /*check dims validity*/
                printf("Invalid Input!");
                exit(1);
            }
        currCol = 0;
        while(cont != NULL){
            cont = strtok(NULL , ",");
            currCol++;
        }
        if ((maxCol < currCol) && (maxCol == 0)){
            maxCol = currCol;
        }
    }
    col = maxCol;
    fclose(myFile1);
}

double dist(double *point, double *centroid){
    double res = 0;
    int i = 0;
    for (i=0; i<col; i++){
         res += pow(*(point+i) - *(centroid+i) , 2);
    }
    return res;
}
  
static PyObject *kmeansC(int k, int max_iter , double eps , PyObject *indList){ 
    /*variable declarations*/
    int i;
    int j;
    int l;
    int cnt = 0;
    int start = 0;
    int stop = 0;
    int flag;
    char *cont;
    char line[1000];
    char lineToSplit[1000];
    double **pointsArray;
    double **sums;
    double **prevAvgs;
    double **currAvgs;
    int *denominators = (int *)calloc(k , sizeof(int));
    int *clusters;
    double *pointMat; 
    double *sumsMat;
    double *prevMat;
    double *currMat;
    FILE *myFile = NULL;
    PyObject *result;
    int *centroids;
    char *readFrom;

    centroids = (int *)calloc(k , sizeof(int)); /*indices will be extracted to this array*/
    readFrom = "new_file";
    for (i = 0 ; i < k; i++){
        centroids[i] = PyLong_AsLong(PyList_GetItem(indList , i));
    }

    extract(readFrom); /*read data*/
    if ((k < 1) || (k > row)){
        printf("Invalid Input!");
        exit(1);
    }    
    clusters = (int *)calloc(row , sizeof(int));
    pointMat = calloc(row*col , sizeof(double));
    pointsArray = calloc(row , sizeof(double *));
    if ((k < 1) || (k > row)){
        printf("Invalid Input!");
        exit(1);
    }
    assert(max_iter > 1 && "Invalid Input!");
    for(i=0; i<row; i++){ /*dynamic mat alloc*/ 
        pointsArray[i] = pointMat + i*col;
    }
    i = 0;
    flag = 1;
    myFile = fopen(readFrom , "r");
    assert(myFile != NULL && "Invalid Input!");
    while (fgets(line , sizeof(lineToSplit) , myFile)){
         j = 0;
         stop = 0;
         cont = strtok(line , ",");        
         while (cont != NULL){
             if (j == col-1){
                 stop = strlen(cont)-1;
             }
             else{
                 stop = strlen(cont);
             }
             start = 0;
             pointsArray[i][j] = toDouble(cont, start, stop-1);
             j++;
             cont = strtok(NULL , ",");
         }
         i++;
    }
    fclose(myFile);
    
    /*start K-mean algo*/

    sumsMat = calloc(k*col, sizeof(double));
    assert(sumsMat != NULL && "An Error Has Occured");
    prevMat = calloc(k*col, sizeof(double));
    assert(prevMat != NULL && "An Error Has Occured");
    currMat = calloc(k*col, sizeof(double));
    assert(currMat != NULL && "An Error Has Occured");
    sums = calloc(k , sizeof(double *));
    assert(sums != NULL && "An Error Has Occured");
    prevAvgs = calloc(k , sizeof(double *));
    assert(prevAvgs != NULL && "An Error Has Occured");
    currAvgs = calloc(k , sizeof(double *));
    assert(currAvgs != NULL && "An Error Has Occured");

    for (i=0; i<k; i++){ /*2D matrix alloc*/
        sums[i] = sumsMat + i*col;
        prevAvgs[i] = prevMat + i*col;
        currAvgs[i] = currMat + i*col;
    }
    for (i=0; i<k; i++){ /*init centroids*/
        for (j=0; j<col; j++){
            currAvgs[i][j] = pointsArray[centroids[i]][j];
        }
    }
    while ((flag == 1) && (cnt < max_iter)){  
        for (i = 0; i < k; i++){ /*init denom to zeroes*/
            denominators[i] = 0;
        }
        for (i=0; i<k; i++){ /*init sums to zeroes*/
            for (j=0; j<col; j++){
                sums[i][j] = 0;
            }
        } 
        for (i=0; i<row; i++){
            double currDist = __DBL_MAX__;
            double minDist = __DBL_MAX__;
            for (j=0; j<k; j++){ /*find closest centroid*/
                currDist = dist((double *)pointsArray[i] , (double *)currAvgs[j]);
                if (currDist < minDist){
                    minDist = currDist;
                    clusters[i] = j;
                }
            }/*end of inner loop for centroids*/
            for (l=0; l<col; l++){ /*adding point's value to relevant centroid*/
                sums[clusters[i]][l] += pointsArray[i][l]; 
            }
            denominators[clusters[i]] += 1;
        }/*end of points loop*/

        for (i=0; i<k; i++){ /*update avgs*/
            for (j=0; j<col; j++){
                currAvgs[i][j] = (sums[i][j])/(denominators[i]);
            }
        }        
        flag = 0;
        for (i = 0; i < k; i++) /*check eps cond*/
        {
            double currNorm = dist((double *)currAvgs[i] , (double *)prevAvgs[i]);
            if (currNorm > pow(eps , 2)){
                flag = 1;
            }
        }
        for (i = 0; i < k; i++){ /*prevAvgs <-- currAvgs*/
            for (j=0; j<col; j++){
                prevAvgs[i][j] = currAvgs[i][j];
            }
        }
        cnt++; /*count iters*/
    }
    /*preparing final centroids for transmissin to python*/
    result = PyList_New(0);
    PyObject *dumby;
    for (i = 0; i < k; ++i){
        dumby = PyList_New(0);
        for (j = 0; j < col; ++j){
            if (PyList_Append(dumby , PyFloat_FromDouble(currAvgs[i][j])) != 0){
                printf("An Error Has Occured");
                exit(1);
            }
        }
        if (PyList_Append(result , dumby) != 0){
                printf("An Error Has Occured");
                exit(1);
            }
    }
    
    free(clusters);
    free(denominators);
    free(pointMat);        
    free(sumsMat);
    free(prevMat);
    free(currMat);
    free(currAvgs);
    free(sums);
    free(prevAvgs);
    free(pointsArray);
    free(centroids);

    return result;
}