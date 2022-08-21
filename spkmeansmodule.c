#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"

double **pyToCmat(PyObject *mat, int x, int y){ /* convert python 2D list to C 2D array */
    
    double **res = mem2D(x, y);
    int i, j;
    
    for (i = 0; i < x; i++){
        for (j = 0; j < y; j++){
            res[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(mat, i), j));
        }
    }
    
    return res;

}

PyObject *cToPyMat(double **mat, int x, int y){

    PyObject *res = PyList_New(x);
    int i, j;
    
    for (i = 0; i < x; i++){
        PyObject *row = PyList_New(y);
        for (j = 0; j < y; j++){
            PyList_SetItem(row, j, PyFloat_FromDouble(mat[i][j]));
        }
        PyList_SetItem(res, i, row);
    }
    
    return res;

}

/*
static PyObject *wamCapi(PyObject *self, PyObject *args){

    double **wamOut, **pointsMat;
    int row, col;
    PyObject *dataPoints, *resW;

    if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &dataPoints)){
        printf("Invalid Input!");
        exit(1);
    }
    if (!PyList_Check(PyList_GetItem(dataPoints, 0))){
        printf("Invalid Input!");
        exit(1);
    }
    
    row = (int)PyList_Size(dataPoints);
    col = (int)PyList_Size(PyList_GetItem(dataPoints, 0));
    pointsMat = pyToCmat(dataPoints, row, col);
    wamOut = mainFuncCapi(pointsMat, "wam", row, col);
    resW = cToPyMat(wamOut, row, col);

    free2D(pointsMat);
    free2D(wamOut);

    return resW;
}

static PyObject *ddgCapi(PyObject *self, PyObject *args){

    double **ddgOut, **pointsMat;
    int row, col;
    PyObject *dataPoints, *resD;

    if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &dataPoints)){
        printf("Invalid Input!");
        exit(1);
    }
    if (!PyList_Check(PyList_GetItem(dataPoints, 0))){
        printf("Invalid Input!");
        exit(1);
    }

    row = (int)PyList_Size(dataPoints);
    col = (int)PyList_Size(PyList_GetItem(dataPoints, 0));
    pointsMat = pyToCmat(dataPoints, row, col);
    ddgOut = mainFuncCapi(pointsMat, "ddg", row, col);
    resD = cToPyMat(ddgOut, row, col);

    free2D(pointsMat);
    free2D(ddgOut);

    return resD;

}

static PyObject *lnormCapi(PyObject *self, PyObject *args){

    double **lnormOut, **pointsMat;
    int row, col;
    PyObject *dataPoints, *resLnorm;

    if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &dataPoints)){
        printf("Invalid Input!");
        exit(1);
    }
    if (!PyList_Check(PyList_GetItem(dataPoints, 0))){
        printf("Invalid Input!");
        exit(1);
    }

    row = (int)PyList_Size(dataPoints);
    col = (int)PyList_Size(PyList_GetItem(dataPoints, 0));
    pointsMat = pyToCmat(dataPoints, row, col);
    lnormOut = mainFuncCapi(pointsMat, "lnorm", row, col);
    resLnorm = cToPyMat(lnormOut, row, col);

    free2D(pointsMat);
    free2D(lnormOut);

    return resLnorm;

}

static PyObject *jacobiCapi(PyObject *self, PyObject *args){

    double **jacobiOut, **pointsMat;
    int row, col;
    PyObject *dataPoints, *resJacobi;

    if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &dataPoints)){
        printf("Invalid Input!");
        exit(1);
    }
    if (!PyList_Check(PyList_GetItem(dataPoints, 0))){
        printf("Invalid Input!");
        exit(1);
    }

    row = (int)PyList_Size(dataPoints);
    col = (int)PyList_Size(PyList_GetItem(dataPoints, 0));
    pointsMat = pyToCmat(dataPoints, row, col);
    jacobiOut = mainFuncCapi(pointsMat, "jacobi", row, col);
    resJacobi = cToPyMat(jacobiOut, row, col);

    free2D(pointsMat);
    free2D(jacobiOut);

    return resJacobi;

}
*/

static PyObject *pyToCtoPy(PyObject *self, PyObject *args){ /* all goals except spk - one back and forth */

    double **outputC, **pointsMat;
    int row, col;
    int goal;
    PyObject *dataPoints, *resPy;

    if (!PyArg_ParseTuple(args, "Oi", &dataPoints, &goal)){ /* assign args */
        printf("Invalid Input C-api1!");
        exit(1);
    }
    if ( (!PyList_Check(dataPoints)) || (!PyList_Check(PyList_GetItem(dataPoints, 0))) ){/* check args */
        printf("Invalid Input C-api2!");
        exit(1);
    }
    
    row = (int)PyList_Size(dataPoints); /* get dimensions */
    col = (int)PyList_Size(PyList_GetItem(dataPoints, 0));
    pointsMat = pyToCmat(dataPoints, row, col); /* convert Py data points to C data points */
    /*goalC = (int)PyInt_AsLong(goal);   convert goal(Py) to goalC(C) */
    outputC = mainFuncCapi(pointsMat, goal, row, col); /* send to main func in C and receive requested output */
    if (goal == 4){ /* 4 means "jacobi" */
        resPy = cToPyMat(outputC, row + 1, row); /* +1 for eigenvalues */
    }
    else{
        resPy = cToPyMat(outputC, row, row);
    }

    free2D(pointsMat);
    free2D(outputC);

    return resPy;
}

static PyObject *spkCapi(PyObject *self, PyObject *args){ /* first back and forth - data points ---> V matrix */
    
    double **outputVmatC, **pointsMat;
    int row, col;
    PyObject *dataPoints, *resPyVmat;

    if (!PyArg_ParseTuple(args, "O", &dataPoints)){ /* assign args */
        printf("Invalid Input spk C-api1!");
        exit(1);
    }
    if ( (!PyList_Check(dataPoints)) || (!PyList_Check(PyList_GetItem(dataPoints, 0))) ){ /* check args */
        printf("Invalid Input  spk C-api2!");
        exit(1);
    }
    
    row = (int)PyList_Size(dataPoints); /* get dimensions */
    col = (int)PyList_Size(PyList_GetItem(dataPoints, 0));
    pointsMat = pyToCmat(dataPoints, row, col); /* convert Py data points to C data points */
    outputVmatC = mainFuncCapi(pointsMat, 5, row, col); /* send to main func in C and receive V matrix */
    resPyVmat = cToPyMat(outputVmatC, row + 1, row);

    /*free2D(pointsMat);
    free2D(outputVmatC);*/

    return resPyVmat; /* send V matrix to Py - next func recieves T matrix and k indices form kmeans++ */
}

static PyObject *eigenGapCapi(PyObject *self, PyObject *args){ /* returns k (eigen heuristic) + sorted V matrix */

    double **outputSortedVmatC, **inputC_V_mat;
    int row, col, k;
    PyObject *inputPy_V_mat, *resPySortedVmat;

    if (!PyArg_ParseTuple(args, "O", &inputPy_V_mat)){ /* assign args */
        printf("Invalid Input eigenGap C-api1!");
        exit(1);
    }
    if ( (!PyList_Check(inputPy_V_mat)) || (!PyList_Check(PyList_GetItem(inputPy_V_mat, 0))) ){ /* check args */
        printf("Invalid Input eigenGap C-api2!");
        exit(1);
    }

    row = (int)PyList_Size(inputPy_V_mat); /* get dimensions */
    col = (int)PyList_Size(PyList_GetItem(inputPy_V_mat, 0));
    printf("\n row = %d , col = %d" , row , col);
    inputC_V_mat = pyToCmat(inputPy_V_mat, row, col); /* convert V matrix to C */
    printf(" input to sort \n");
    print2Darray(inputC_V_mat, row, col);
    outputSortedVmatC = sortMat(inputC_V_mat); /* sort V matrix by eigen values */
    printf(" sorted jacobi \n");
    print2Darray(outputSortedVmatC , row , col);
    resPySortedVmat = cToPyMat(outputSortedVmatC, row, col); /* convert sorted V to Py */
    k = eigenGap(outputSortedVmatC); /* obtain k by hueristic */

    /*free2D(inputC_V_mat);
    free2D(outputSortedVmatC);*/

    return Py_BuildValue("(iO)" , k , resPySortedVmat); /* return k , sorted V matrix */

}

static PyObject *kmeansCapi(PyObject *self, PyObject *args){ /* second back and forth - T mat + k indices ---> k clusters */
    
    double **tMatC, **outputClustersC/*, **initCentroidsC*/;
    int *centroidsArrayC;
    int row, col, k, i;
    PyObject *tMatPy, *initCentroidsPy, *resClustersPy;
    if (!PyArg_ParseTuple(args, "OOi", &tMatPy, &initCentroidsPy, &k)){ /* assign args */
        printf("Invalid Input kmeans C-api1!");
        exit(1);
    }
    if ( (!PyList_Check(tMatPy)) || (!PyList_Check(PyList_GetItem(tMatPy, 0))) || (!PyList_Check(initCentroidsPy)) ){ /* check args */
        printf("Invalid Input  kmeans C-api2!");
        exit(1);
    }
    row = (int)PyList_Size(tMatPy);
    col = (int)PyList_Size(PyList_GetItem(tMatPy, 0));
    tMatC = pyToCmat(tMatPy, row, col); /* convert T mat from Py to C */
    /*initCentroidsC = pyToCmat(initCentroidsPy, 1, k);  convert k init centroids from Py to C */
    centroidsArrayC = (int*)calloc(k, sizeof(int));
    for (i = 0; i < k; i++){
	centroidsArrayC[i] = (int)PyFloat_AsDouble(PyList_GetItem(initCentroidsPy, i));
    }
    outputClustersC = kmeans(tMatC, centroidsArrayC , row, col); /* send T matrix to main func to calculate k clusters in kmeans algo */
    resClustersPy = cToPyMat(outputClustersC, k, col); /* receive clusters from C */
    /*free2D(tMatC);*/
    free2D(outputClustersC);
    return resClustersPy;
}

static PyMethodDef capiMethods[] = { /* an array of the module methods */ 
    {"pyToCtoPy", (PyCFunction)pyToCtoPy, METH_VARARGS, PyDoc_STR("receive input from Py and outputs C to Py")},
    {"spkCapi", (PyCFunction)spkCapi, METH_VARARGS, PyDoc_STR("first back and forth - data points ---> V matrix")}, 
    {"eigenGapCapi", (PyCFunction)eigenGapCapi, METH_VARARGS, PyDoc_STR("returns k (eigen heuristic) + sorted V matrix")},
    {"kmeansCapi", (PyCFunction)kmeansCapi, METH_VARARGS, PyDoc_STR("second back and forth - T mat + k indices ---> k clusters")},
    {NULL, NULL, 0, NULL} 
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "spkmeansmodule", "spkmeans module", -1, capiMethods
};

PyMODINIT_FUNC PyInit_spkmeansmodule(void)
{
    PyObject *m;

    m = PyModule_Create(&moduledef);
    if (!m)
    {
        return NULL;
    }
    return m;
}