#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/*declarations*/
void print2Darray(double**, int, int);
double** mem2D(int, int);
void free2D(double **);
double norm(double*, double*, int);
double toDouble(const char*, int, int);
double** extractFromFile(char*);
double** wam(double**);
double** ddg(double**);
double** leftDiagMult(double**, double**);
double** rightDiagMult(double**, double**);
double** invertRoot(double**);
double** lnorm(double**);
double** newP(int, int, double, double);
double* obtain(double, double, double); /*obtain s,c*/
int convergance(double**, double**);
double** updateA(double**, double**, int, int, double, double);
double** updateVects(double**, double**);/*P matrix multipication*/
double* find_I_J(double** aMat);
double** jacobi(double**);
double** sortMat(double**);
void swap(int*, int*);
int eigenGap(double**);
double** trans(double**, int , int);
double dist(double *, double *);
double **kmeans(double **, int *, int, int);
double **mainFuncCapi(double **, int , int , int );

/*global variables*/
int row, col;

/*functions*/
void print2Darray(double** a, int x, int y){
    int rows, columns;
 	for(rows = 0; rows < x; rows++)
  	{
        printf("row num: %d  " , rows); /* delete on submission !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
  		for(columns = 0; columns < y; columns++)
  		{
  			if (columns == y - 1){
                printf("%.4f \n" , a[rows][columns]);
            } else{
                printf("%.4f,", a[rows][columns]);
            }
		}
  	}  	
}

double** mem2D(int x, int y){
    double** mat;
    int i;
    mat = (double**)calloc(x, sizeof(double*));
    assert(mat != NULL && "An Error Has Occured");
    for(i=0; i<x; i++){ 
        mat[i] = (double*)calloc(y, sizeof(double));
        assert(mat[i] != NULL && "An Error Has Occured");
    }
    return mat;
}

void free2D(double **mat){
    free(*mat);
    free(mat);
}

double norm(double *point1, double *point2 , int x){
    double res = 0;
    int i = 0;
    for (i=0; i<x; i++){
         res += pow(*(point1+i) - *(point2+i) , 2);
    }
    res = pow(res , 0.5);
    return res;
}

double toDouble(const char* s, int start, int stop) {
    int i;
    int dot = 1;
    int m = 1;
    double ret = 0;
    stop = stop > 10 ? 10 : stop;
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
            printf("Invalid Input!D");
            exit(1);
        }
        ret += (s[i] - '0') * m;
        m *= 10;
    }
    ret /= dot;
    return ret;
}

double ** extractFromFile(char* fileName){

    int i, j;    
    char* cont;
    char line[1000]; 
    char lineToSplit[1000];
    double **pointsArray;
    int maxCol = 0;
    int currCol = 0;
    int stop = 0;
    int start = 0;
    FILE* myFile1 = NULL;
    FILE* myFile = NULL;
    myFile1 = fopen(fileName, "r");
    if (myFile1 == NULL){
        printf("Invalid Input!1");
        exit(1);
    }
    while (fgets(line , sizeof(line) , myFile1)){ /*get dims*/
        row++;
        cont = strtok(line , ",");
        if ((maxCol != currCol) && (maxCol != 0)){ /*check dims validity*/
                printf("Invalid Input!2");
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
    /*dynamic mat alloc*/ 
    pointsArray = mem2D(row, col);
    myFile = fopen(fileName, "r");
    assert(myFile != NULL && "Invalid Input!3");
    i=0;
    while (fgets(line , sizeof(lineToSplit) , myFile)){ /*initial points matrix*/
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
             pointsArray[i][j] = toDouble(cont, start, stop - 1);
             j++;
             cont = strtok(NULL , ",");
         }
         i++;
    }
    fclose(myFile);
    return pointsArray;
}

double** wam(double** dataPoints){
    double weight;
    int i,j;
    double** weightsMat;

    weightsMat = mem2D(row , row);

    for (i = 0; i < row; i++){
        for (j = i; j < row; j++){
            if (i == j) {
                weightsMat[i][j] = 0;
            }
            else {
                weight = exp(norm(dataPoints[i] , dataPoints[j], col) / -2); 
                /*printf("\n i = %d, j = %d, weight = %.10f" , i , j , weight);*/
                weightsMat[i][j] = weight;
                weightsMat[j][i] = weight;
            }
        }
    }
    return weightsMat;
}

double** ddg(double** dataPoints){
    double** diagMat;
    double** weightsMat;
    int i, j;
    double sum = 0;
    weightsMat = mem2D(row, row);
    diagMat = mem2D(row, row);
    weightsMat = wam(dataPoints);
    for(i = 0; i < row; i++){
        sum = 0;
        for (j = 0; j < row; j++){
            sum += weightsMat[i][j];
        }
        diagMat[i][i] = sum;
    }
    return diagMat;
}

double** leftDiagMult(double** diagMat, double** otherMat){ /* other[i][j] = other[i][j] * diag[i][i] */
    int i, j;
    double** res;

    res = mem2D(row, row);

    for (i = 0; i < row; i++){
        for (j = 0; j < row; j++){
            res[i][j] = otherMat[i][j] * diagMat[i][i];
        }
    }
    return res;
}

double** rightDiagMult(double** otherMat, double** diagMat){ /* other[i][j] = other[i][j] * diag[j][j] */
    int i, j;
    double** res;

    res = mem2D(row, row);

    for (i = 0; i < row; i++){
        for (j = 0; j < row; j++){
            res[i][j] = otherMat[i][j] * diagMat[j][j];
        }
    }
    return res;
}

double** invertRoot(double** diagMat){/* testing */
    double** invSqrtMat;
    int i;
    invSqrtMat = mem2D(row, row);
    for (i = 0; i < row; i++){
        invSqrtMat[i][i] = pow(diagMat[i][i] , -0.5);
    }
    return invSqrtMat;
}

double** lnorm(double** datapoints){
    double** diagMat;
    double** mat_L;
    double** weightsMat;
    double** invSqrtMat;
    double** res;
    int i, j;

    diagMat = mem2D(row, row);
    mat_L = mem2D(row, row);
    weightsMat = mem2D(row, row);
    invSqrtMat = mem2D(row, row);
    res = mem2D(row, row);

    weightsMat = wam(datapoints);
    diagMat = ddg(datapoints);

    for (i = 0; i < row; i++){ /*create L matrix*/
        for (j = 0; j< row; j++){
            mat_L[i][j] = diagMat[i][j] - weightsMat[i][j];
        }
    }
    for (i = 0; i < row; i++){
        invSqrtMat[i][i] = pow(diagMat[i][i] , -0.5);
    }
    res = rightDiagMult( leftDiagMult(invSqrtMat , mat_L) , invSqrtMat); /*triple matrix multiplication*/
    return res;
}

double* obtain(double aij, double aii, double ajj){ /* get s and c */
    double* result;
    int sign;
    double theta, t;
    if (aij == 0.0){
        theta = strtod("Inf", NULL);
        /*result[0] = 1;
        result[1] = 0;
        return result;*/
    } else {
        theta = (ajj-aii)/(2*aij);
    }
    sign = theta>=0 ? 1 : -1;
    t = sign/(fabs(theta) + pow((pow(theta, 2) + 1), 0.5));
    result = (double*)calloc(2, sizeof(double));
    assert(result != NULL && "An Error Has Occured");
    result[0] = pow(pow(t,2) + 1, -0.5);
    result[1] = t*result[0];
    /*printf("\nt is %f c is %f s is %f", t, result[0], result[1]);*/
    return result;
}

double** newP(int iInd, int jInd, double c, double s){ /* create rotation matrix P */
    double** pMat;
    int i;
    pMat = mem2D(row,row);
    for(i=0; i<row; i++){
        pMat[i][i] = 1;
    }
    pMat[iInd][iInd] = c;
    pMat[jInd][jInd] = c;
    pMat[iInd][jInd] = s;
    pMat[jInd][iInd] = -s;

    return pMat;
}

double** updateA(double** aMat, double** prevA, int iInd, int jInd, double c, double s){ /* create A' */
    int r;
    for(r=0; r<row; r++){
        if ((r!=iInd) && (r!=jInd)){
            aMat[r][iInd] = ((c*prevA[r][iInd]) - (s*prevA[r][jInd]));
            aMat[r][jInd] = ((c*prevA[r][jInd]) + (s*prevA[r][iInd]));
        } 
    }
    aMat[iInd][iInd] = ((c*c*prevA[iInd][iInd]) + (s*s*prevA[jInd][jInd]) - (2*s*c*prevA[iInd][jInd]));
    aMat[jInd][jInd] = ((s*s*prevA[iInd][iInd]) + (c*c*prevA[jInd][jInd]) + (2*s*c*prevA[iInd][jInd]));
    aMat[iInd][jInd] = (((c*c-s*s)*prevA[iInd][jInd]) + (s*c*(prevA[iInd][iInd])-(prevA[jInd][jInd])));
    return aMat;
}

double* find_I_J(double** aMat){ /* get aij, indices */
    double aij, iRes, jRes;
    double* vals;
    int i, j, iInd, jInd;
    aij = 0.0;
    iRes = 0.0; /* added !!!!!!!!!!! */
    jRes = 0.0;
    for(i = 0; i<row; i++){
        for(j=i; j<row; j++){
            if ( (fabs(aMat[i][j])>fabs(aij)) && (i!=j) ){
                aij = aMat[i][j];
                iRes = i;
                jRes = j;
            }
        }
    }
    vals = (double*)calloc(3, sizeof(double));
    assert(vals != NULL && "An Error Has Occured");
    iInd = (int)iRes;
    jInd = (int)jRes; /* for comparison */
    if ((iInd == jInd) && (row != 1)){ /* in case of I matrix */
        iRes = 0.0;
        jRes = 1.0;
    }
    vals[0]= aij;
    vals[1] = iRes;
    vals[2] = jRes;

    return vals;
}

int convergance(double** aMat, double** prevA){
    int i, j;
    double offA = 0;
    double offPrev = 0;
    double eps = 1.0 * pow(10, -5);
    for(i=0; i<row; i++){
        for(j=0; j<row; j++){
            if (i!=j){
                offA += pow(aMat[i][j], 2);
                offPrev += pow(prevA[i][j], 2);
            }
        }
    }
    if ((offPrev - offA) <= eps){
        return 1;
    }
    return 0;

}

double** updateVects(double** tempV, double** pMat){ /* matrix multiplication */
    double** res;
    int i,j,k;
    res = mem2D(row, row);
    for(i=0; i<row; i++){
        for(j=0; j<row; j++){
            for(k=0; k<row; k++){
                res[i][j] += tempV[i][k]*pMat[k][j];
            }
        }
    }
    return res;
}

double** trans(double** mat,int x ,int y){
    double** res = mem2D(x, y);
    int i, j;
    for(i = 0; i < x; i++){
        for (j = 0; j < y; j++){
            res[i][j] = mat[j][i];    
        }
    }
    return res;
}

void swap(int *a, int *b){
    int tmp;
    tmp = *a;
    *a = *b;
    *b = tmp;
}

double** sortMat(double** eigenMat){ /* returns the matrix of eigenvals/vecs sorted decreasingly */
    double **tmp1Mat, **tmp2Mat, **res;
    int* indices;
    int i, j;
    indices = (int*)calloc(row, sizeof(int));
    assert(indices!= NULL && "An Error Has Occured");
    for(i=0; i<row; i++){
        indices[i] = i;
    }
    tmp1Mat = mem2D(row , row);
    for(i=0; i<row - 1; i++){/*bubble sort indices*/
        for(j=0; j<(row - i - 1); j++){
            if (eigenMat[0][indices[j]] > eigenMat[0][indices[j+1]]){
                swap(&indices[j], &indices[j+1]);               
            } 
        }
    }
    for(i=0; i<row; i++){/*create E transpose, vectors are rows*/
        for(j = 1; j<row + 1; j++){
            tmp1Mat[i][j-1] = eigenMat[j][i];
        }
    }
    tmp2Mat = mem2D(row, row);
    for (i = 0; i < row; i++){ /* sort vectors rows by eigenvalues */
        for (j = 0; j < row; j++){
            tmp2Mat[i][j] = tmp1Mat[indices[row - i - 1]][j];
        }
    }
    tmp2Mat = trans(tmp2Mat, row, row); /* vectors are columns */
    res = mem2D(row + 1, row);
    for (i = 0; i < row; i++){
        res[0][i] = eigenMat[0][indices[row - i - 1]]; /* value for each column */
    }
    for (i = 1; i < row + 1; i++){ /* column vcector for each value */
        for (j = 0; j < row; j++){
            res[i][j] = tmp2Mat[i - 1][j];
        }
    }
    return res;
}

int eigenGap(double** sortedVmat){
    int i, k;
    double maxGap = 0;
    double *gaps = (double*)calloc(row / 2 , sizeof(double)); 
    k = 0; /* added !!!!!!!!!!!!!!!!!!!! */
    for(i = 0; i < row / 2; i++){
        gaps[i] = fabs(sortedVmat[0][i] - sortedVmat[0][i+1]);
    }
    for (i = 0; i < row / 2; i++){
        if (gaps[i] > maxGap){
            maxGap = gaps[i];
            k = i;
        }
    }
    return k + 1;
}

double** jacobi(double** aMat){
    double **vResult, **pMat, **tempV, **prevA;
    double *parameters, *c_s;
    double aij, c, s;
    int i,j, iRes, jRes;
    int flag = 0;
    int cnt = 0;
    parameters = (double*)calloc(3, sizeof(double));
    assert(parameters != NULL && "An Error Has Occured");
    c_s = (double*)calloc(2, sizeof(double));
    assert(c_s != NULL && "An Error Has Occured");
    tempV = mem2D(row,row);
    pMat = mem2D(row,row);
    prevA = mem2D(row,row);
    for(i=0; i<row; i++){/*init V*/
        tempV[i][i] = 1;
    }
    while ((flag == 0) && (cnt<100)) { /*Jacobi main loop*/
        parameters = find_I_J(aMat);
        aij = parameters[0];
        iRes = (int)parameters[1];
        jRes = (int)parameters[2];
        if (iRes == jRes){ /* single entry matrix */
            vResult = mem2D(2, 1);
            vResult[0][0] = aMat[0][0];
            vResult[1][0] = 1.0;
            return vResult;
        }
        /*printf("\n aij = %.4f, i = %d, j = %d" , aij, iRes, jRes); */
        c_s = obtain(aij, aMat[iRes][iRes], aMat[jRes][jRes]);
        c = c_s[0];
        s = c_s[1];
        pMat = newP(iRes, jRes, c, s);
        /*printf("\n\nP'%d is ", cnt);
        print2Darray(pMat,row,row);*/
        for(i=0; i<row; i++){/*save A values*/
            for(j=0; j<row; j++){
                prevA[i][j] = aMat[i][j];
            }
        }
        /*aMat = updateA(aMat, prevA, iRes, jRes, c, s);*/
        aMat = updateVects(trans(pMat, row, row), updateVects(prevA, pMat));
        /*printf("\n\nA'%d is \n", cnt);
        print2Darray(aMat, row, row); */
        tempV = updateVects(tempV, pMat); /*   multiply V * Pi   */
        flag = convergance(aMat, prevA); /* check condintion */
        cnt+=1;
    }
    printf("\n %d iterations \n", cnt);
    vResult = mem2D((row+1), row); 
    for(i=0; i<row; i++){ /* insert eigenvalues to V */
        vResult[0][i] = aMat[i][i];
    }
    for(i=1; i<=row; i++){ /* insert eigenvectors to V */
        for(j=0; j<row; j++){
            vResult[i][j] = tempV[i-1][j];
        }
    }
    return vResult;
}

double dist(double *point, double *centroid){
    double res = 0;
    int i = 0;
    for (i=0; i<col; i++){
         res += pow(*(point+i) - *(centroid+i) , 2);
    }
    return res;
}

double **kmeans(double **pointsArray, int *initCentroids, int rows, int col){

    int i, j, l, k, flag;
    int *clusters, *denominators;
    double /***pointsArray, */**sums, **prevAvgs, **currAvgs;
    int cnt = 0;
    int max_iter = 300; /* fix to 300 !!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
    i = 0;
    flag = 1;
    /* sizeof(initCentroids) / sizeof(int );*/
    k = col;
    row = rows;
    denominators = (int *)calloc(k , sizeof(int));
    clusters = (int *)calloc(row , sizeof(int));
    /*pointsArray = mem2D(row, col);*/
    sums = mem2D(k, col);
    prevAvgs = mem2D(k, col);
    currAvgs = mem2D(k, col);
    /*pointsArray = tMat; */
    for (i=0; i<k; i++){ /*init centroids*/
            for (j=0; j<col; j++){
                currAvgs[i][j] = pointsArray[initCentroids[i]][j];
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
            if (currNorm > 0){ /* changed to sqrt + change eps to 0 !!!!!!!!!!!!!!!!!!!!!!!!!!*/
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
    printf("iters = %d \n" , cnt);

    free(clusters);
    free(denominators);
    free2D(sums);
    free2D(prevAvgs);
    free2D(pointsArray);

    return currAvgs;

}

double **mainFuncCapi(double **inputMat, int goal, int x, int y){ /* input from C api, output to C api */
    double **res, **lnormOutput;
    row = x;
    col = y;
    if (goal == 1){
        res = mem2D(x, x);
        res = wam(inputMat);
        return res;
    }
    if (goal == 2){
        res = mem2D(x, x);
        res = ddg(inputMat);
        return res;
    }
    if (goal == 3){
        res = mem2D(x, x);
        res = lnorm(inputMat);
        return res;
    }
    if (goal == 4){
        res = mem2D(x + 1, x);
        res = jacobi(inputMat); /* this inputMat is a symetric Lnorm and not a regular data set */
        return res; 
    }
    if (goal == 5){
        lnormOutput = mem2D(x, x);
        res = mem2D(x + 1, x);
        lnormOutput = lnorm(inputMat);

	printf("\n %.4f\n" , lnormOutput[0][0]);
        printf(" jacobi input is \n"); /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
        print2Darray(lnormOutput, row, row);
        res = jacobi(lnormOutput); /* change from "lnormOutput" to "inputMat" in order to check jacobi*/ 
        printf(" jacobi output is \n");
        print2Darray(res, row + 1, row);

        return res;
    }
    return inputMat; /* should not get here */
}

int main(int argc, char *argv[]){

    /*int i; romove !!!!!!!!!*/
    /*int *initCentroids;  delete later !!!!!!!!!!!!! */

    double **dataPoints, **res;
    char *goal = argv[1];
    char *fileName = argv[2];

    if (argc != 3){ 
        printf("Invalid Input!Main");
        exit(1);
    }
    if ( (strcmp(goal, "wam") != 0) && (strcmp(goal, "ddg") != 0) && (strcmp(goal, "lnorm") != 0) && (strcmp(goal, "jacobi") != 0) ){  /*check goal validity*/ 
        printf("Invalid Input!Main");
        exit(1);
    }

    dataPoints = mem2D(row, row);
    dataPoints = extractFromFile(fileName);

    /*if(strcmp(goal, "kmeansC") == 0){ 
        res = mem2D(15, col);
        initCentroids = (int*)calloc(15,sizeof(int));
        for(i = 0; i< 15; i++){
            initCentroids[i] = i;
        }
        initCentroids[0] = 2732;
        initCentroids[1] = 2982;
        initCentroids[2] = 4229;
        initCentroids[3] = 4312;
        initCentroids[4] = 4232;
        initCentroids[5] = 3111;
        initCentroids[6] = 1902;
        initCentroids[7] = 1466;
        initCentroids[8] = 272;
        initCentroids[9] = 1354;
        initCentroids[10] = 2369;
        initCentroids[11] = 3996;
        initCentroids[12] = 2359;
        initCentroids[13] = 1961;
        initCentroids[14] = 4117;
        res = kmeans(dataPoints,initCentroids);
        print2Darray(res, 15, col);
        return 0;
    }*/

    if (strcmp(goal, "jacobi") == 0){
        res = mem2D(row + 1, row);
        res = jacobi(dataPoints);
        print2Darray(res, row + 1, row);
        free2D(dataPoints);
        free2D(res);
        return 0;
    }

    res = mem2D(row, row);
    if (strcmp(goal, "wam") == 0){
        res = wam(dataPoints);
    }
    if (strcmp(goal, "ddg") == 0){
        res = ddg(dataPoints);
    }
    if (strcmp(goal, "lnorm") == 0){
        res = lnorm(dataPoints);
    }
    print2Darray(res, row, row);
    free2D(dataPoints);
    free2D(res);
    return 0; 

    /*int k;
    double** pointsMat;
    double** weightsMat;
    double** diagMat;
    double** lNorm;
    double** vResult;
    double **sortedV;

    pointsMat = mem2D(row, col);
    weightsMat = mem2D(row, row);
    diagMat = mem2D(row, row);
    lNorm = mem2D(row , row);
    vResult = mem2D(row, row);
    sortedV = mem2D(row, row);
    printf("\nMain");
    pointsMat = extractFromFile(argv[1]);
    printf("\n***********************data points*************************************");
    print2Darray(pointsMat, row, col);
    weightsMat = wam(pointsMat); 
    printf("\n***********************   W    *************************************");
    print2Darray(weightsMat, row, row);
    diagMat = ddg(pointsMat); 
    printf("\n***********************    D    *************************************");
    print2Darray(diagMat, row, row);
    lNorm = lnorm(pointsMat);  
    printf("\n*****************           L Norm          ***************");
    print2Darray(lNorm, row, row);
    vResult = jacobi(lNorm);
    printf("\n***************      Jacobi - V(n+1 x n) not sorted    **********");
    print2Darray(vResult, (row+1), row);
    sortedV = sortMat(vResult);
    printf("\n*************      Vsorted        *******************");
    print2Darray(sortedV, row+1, row);
    k = eigenGap(sortedV);
    printf("\n k = %d" , k);

    return argc - argc;*/
}