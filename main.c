//Fernanda Bezerra Nascimento - 388834
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Foi utilizada a normalização euclidiana

int Columns    = 2;
int Lines      = 2;
int VectorSize = 2;


double** MatrixAlocation(int Lines, int Columns){
    double **Matrix = (double**)malloc(Lines * sizeof(double*));

    for (int i = 0; i < Lines; i++){
        Matrix[i] = (double*) malloc(Columns * sizeof(double));

    }
    return Matrix;
}

void PrintMatrix(double **Matrix){
    for (int l = 0; l < Lines; ++l) {
        for (int m = 0; m < Columns; ++m) {
            printf("%f", Matrix[l][m]);
            printf("\t");
        }
        printf("\n");
    }
}

double EuclidianNormalization(double *AutoVector) {
    double Normalization = 0;
    for (int j = 0; j < VectorSize; ++j) {
        Normalization += pow(AutoVector[j], 2);
    }
    Normalization = sqrt(Normalization);
    return Normalization;
}

double* MatrixVectorMultiplication(double **Matrix, const double *AutoVector){
    double *ResultingVector = (double*) malloc(VectorSize * sizeof(double));
    for (int i = 0; i < Lines ; ++i) {
        ResultingVector[i] = 0;
        for (int j = 0; j < Columns; ++j) {
            ResultingVector[i] = AutoVector[j] * Matrix[i][j];
        }
    }
    return ResultingVector;
}

double* CalculatingQ(double *vector) {
    double *q = (double*) malloc(VectorSize * sizeof(double));
    for (int i = 0; i < VectorSize; ++i) {
        q[i] = vector[i]/EuclidianNormalization(vector);
    }
}

double VectorTransposeMultiplication(const double *AutoVector1, const double *AutoVector2) {
    double Result = 0;
    for (int i = 0; i < VectorSize ; ++i) {
            Result += AutoVector1[i] * AutoVector2[i];

    }
    return Result;
}


double* RegularPow(double **Matrix, double *InitialVector, double Tolerance){
    double *q = CalculatingQ(InitialVector);
    double Error;
    double Lambda;
    double PreviousLambda = 0;
    double *x = MatrixVectorMultiplication(Matrix, q);

    do{
        q = CalculatingQ(x);
        x = MatrixVectorMultiplication(Matrix, x);
        Lambda = VectorTransposeMultiplication(q,x);
        Error = fabs((Lambda - PreviousLambda)/Lambda);
    }
    while (Error > Tolerance );
    return x;
}

double** LUDecompositionForInverse(double **Matrix){
    double              **L = MatrixAlocation(Lines, Columns);
    double              **U = MatrixAlocation(Lines, Columns);
    double              **Y = MatrixAlocation(Lines, Columns);
    double **IdentityMatrix = MatrixAlocation(Lines, Columns);
    double **InverseMatrix  = MatrixAlocation(Lines, Columns);

    for (int l = 0; l < Lines ; ++l) {
        for (int m = 0; m < Columns; ++m) {
            L[l][m] = 0;
            U[l][m] = 0;
            if (l == m){
               IdentityMatrix[l][m] = 1;
            }
            else{
                IdentityMatrix[l][m] = 0;
            }
        }
    }
    printf("IdentityMatrix:\n");
    PrintMatrix(IdentityMatrix);
    printf("\n");

    for(int i = 0; i < Lines; i++) {
        for(int j = 0; j < Columns; j++) {
            L[i][j] = IdentityMatrix[i][j];

            if(i <= j) {
                U[i][j] = Matrix[i][j];
                for(int k = 0; k < i; k++) {
                    U[i][j] -= L[i][k] * U[k][j];
                }
            }

            else {
                L[i][j] = Matrix[i][j];
                for(int k = 0; k < j; k++){
                    L[i][j] -= L[i][k] * U[k][j];
                    L[i][j] /= U[j][j];
                }
                U[i][j]  = 0;
            }
        }
    }
    printf("L:\n");
    PrintMatrix(L);
    printf("\n");
    printf("U:\n");
    PrintMatrix(U);
    printf("\n");

    for(int c = 0; c < Columns; c++) {
        for (int i = 0; i < Lines; ++i) {
            Y[i][c] = IdentityMatrix[i][c];
            for(int k = 0; k < i; k++) {
                Y[i][c] -= L[i][k] * Y[k][c];
            }
        }
    }
    printf("\nY:\n");
    PrintMatrix(Y);
    printf("\n");

    for(int c = 0; c < Columns; c++) {
        for (int i = Lines - 1; i >= 0; --i) {
            InverseMatrix[i][c] = Y[i][c];
            for(int k = 0; k < i; k++) {
                InverseMatrix[i][c] -= U[i][k] * InverseMatrix[k][c];
            }
            InverseMatrix[i][c] /= U[i][i];
        }
    }
    return InverseMatrix;
}

double* InversePow(double **Matrix, double *InitialVector, double Tolerance){
    double **InverseMatrix = LUDecompositionForInverse(Matrix);
    return RegularPow(InverseMatrix, InitialVector, Tolerance);
}


int main() {
    double **Matrix = MatrixAlocation(Lines, Columns);
    Matrix[0][0] = 1;
    Matrix[0][1] = 2;
    Matrix[1][0] = 3;
    Matrix[1][1] = 0;
    double **InverseMatrix = LUDecompositionForInverse(Matrix);
    PrintMatrix(InverseMatrix);

    printf("Hello, World!\n");
    system("pause");
    return 0;
}