//Fernanda Bezerra Nascimento - 388834
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Definição das structs
struct EigenValueVector{
    double   eigenValue;
    double *eigenVector;
};

struct SetOfEigenValueVector{
    double   *eigenValues;
    double **eigenVectors;
};

struct HouseHolderAnswer{
    double **tridiagonalMatrix;
    double **houseHolderMatrix;
};

struct QRMatrix{
    double **qMatrix;
    double **rMatrix;
};

//Definição do tamanho da matriz e do vetor

int columns    = 4;
int lines      = 4;
int vectorSize = 4;

//Foi utilizada a normalização euclidiana

/// Função para alocação do espaço de um vetor
/// \return vetor com o espaço alocado
double* VectorAllocation(){
    double *vector = (double*) malloc(vectorSize * sizeof(double));
    return vector;
}

/// Função para inicializar o vetor com 0
/// \param vector
void InitializeVector(double *vector) {
    for (int l = 0; l < lines ; ++l) {
        vector[l] = 0;
    }
}

/// Função para alocação do espaço para uma Matriz
/// \param lines = quantidade de linhas
/// \param columns = quantidade de colunas
/// \return a matriz com o espaço alocado
double** MatrixAllocation(int lines, int columns){
    double **matrix = (double**)malloc(lines * sizeof(double*));

    for (int i = 0; i < lines; i++){
        matrix[i] = (double*) malloc(columns * sizeof(double));
    }
    return matrix;
}

/// Imprime a Matriz
/// \param matrix
void PrintMatrix(double **matrix){
    for (int l = 0; l < lines; ++l) {
        for (int m = 0; m < columns; ++m) {
            printf("%f", matrix[l][m]);
            printf("\t");
        }
        printf("\n");
    }
    printf("\n");
}

/// Imprime Vetor
/// \param vector
void PrintVector(double *vector){
    for (int i = 0; i < lines; ++i) {
        printf("%f", vector[i]);
        printf("\n");
    }
    printf("\n");
}

/// Normalização do vetor
/// \param vector
/// \return a norma do vetor
double EuclidianNormalization(double *vector) {
    double normalization = 0;
    for (int j = 0; j < vectorSize; ++j) {
        normalization += pow(vector[j], 2);
    }
    normalization = sqrt(normalization);
    return normalization;
}

/// Faz a normalização euclidiana de uma matriz
/// \param matrix
/// \return a norma da matriz
double MatrixEuclidianNormalization(double **matrix){
    double normalization = 0;
    for (int i = 0; i < lines; ++i) {
        for (int j = 0; j < columns; ++j) {
            normalization += pow(matrix[i][j], 2);
        }
    }
    normalization = sqrt(normalization);
    return normalization;
}

/// Multiplica matriz e vetor
/// \param matrix
/// \param vector
/// \return o vetor resultante da multiplicação
double* MatrixVectorMultiplication(double **matrix, const double *vector){
    double *resultingVector = (double*) malloc(vectorSize * sizeof(double));
    for (int i = 0; i < lines ; ++i) {
        resultingVector[i] = 0;
        for (int j = 0; j < columns; ++j) {
            resultingVector[i] += vector[j] * matrix[i][j];
        }
    }
    return resultingVector;
}

/// Faz a matriz diagonal de uma dada matriz
/// \param matrix
/// \return matriz diagonal
double** MakeDiagonalMatrix(double **matrix){
    double **diagonalMatrix = MatrixAllocation(lines, columns);
    for (int l = 0; l < lines ; ++l) {
        for (int m = 0; m < columns; ++m) {
            if (l == m){
                diagonalMatrix[l][m] = matrix[l][m];
            }
            else{
                diagonalMatrix[l][m] = 0.0;
            }
        }
    }
    return diagonalMatrix;
}

/// Normaliza o vetor
/// \param vector
/// \return o vetor normalizado
double* NormalizingVector(double *vector) {
    double *normalizedVector = (double*) malloc(vectorSize * sizeof(double));
    double     normalization = EuclidianNormalization(vector);
    for (int i = 0; i < vectorSize; ++i) {
        normalizedVector[i] = vector[i]/normalization;
    }
    return normalizedVector;
}

/// Multiplicação de vetores
/// \param vector1
/// \param vector2
/// \return o resultado da multiplicacao
double VectorMultiplication(const double *vector1, const double *vector2) {
    double result = 0;
    for (int i = 0; i < vectorSize ; ++i) {
        result += vector1[i] * vector2[i];
    }
    return result;
}

/// Multiplicacao de um vetor transposto por um vetor normal
/// \param vector1
/// \param vector2
/// \return Matriz resultante
double** VectorTranposeVectorMultiplication(const double *vector1, const double *vector2){
    double **resultingMatrix  = MatrixAllocation(lines, columns);
    for (int i = 0; i < columns; ++i) {
        for (int j = 0; j < lines ; ++j) {
            resultingMatrix[j][i] = vector1[j] * vector2[i];
        }
    }
    return resultingMatrix;
}

/// Subtração entre dois vetores
/// \param vector1
/// \param vector2
/// \return Vetor resultante
double* VectorSubtraction(const double *vector1, const double *vector2){
    double *resultingVector = VectorAllocation();
    for (int i = 0; i < lines; ++i) {
        resultingVector[i] = vector1[i] - vector2[i];
    }
    return resultingVector;
}

/// Transforma uma matriz numa matriz identidade
/// \param matrix
void MakeIdentityMatrix(double **matrix) {
    for (int l = 0; l < lines ; ++l) {
        for (int m = 0; m < columns; ++m) {
            if (l == m){
                matrix[l][m] = 1.0;
            }
            else{
                matrix[l][m] = 0.0;
            }
        }
    }
}

/// Função para calcular a matriz resultante da multiplicação de uma matriz por um valor escalar
/// \param matrix
/// \param value
/// \return Matriz resultado da multiplicação
double** MatrixValueMultiplication(double **matrix, double value){
    double **resultingMatrix = MatrixAllocation(lines, columns);
    for (int l = 0; l < lines ; ++l) {
        for (int m = 0; m < columns; ++m) {
            resultingMatrix[l][m] = matrix[l][m] * value;

        }
    }
    return resultingMatrix;
}

/// Função para multiplicar matrizes
/// \param matrix1
/// \param matrix2
/// \return Matriz resultado da multiplicação
double** MatrixMultiplication(double **matrix1, double **matrix2){
    double **resultingMatrix = MatrixAllocation(lines, columns);
    for (int i = 0; i < lines ; ++i) {
        for (int j = 0; j < columns ; ++j) {
            resultingMatrix[i][j] = 0;
            for (int k = 0; k < columns; ++k) {
                resultingMatrix[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
    return resultingMatrix;
}

/// Transpõe uma matriz
/// \param matrix
/// \return a matriz transposta
double** MatrixTranspose(double **matrix){
    double **transposeMatrix = MatrixAllocation(lines, columns);
    for(int i=0; i < lines; ++i)
        for(int j=0; j < columns; ++j) {
            transposeMatrix[j][i] = matrix[i][j];
        }
    return transposeMatrix;
}

/// Função que calcula a matrix resultante da subtração de duas matrizes
/// \param matrix1
/// \param matrix2
/// \return Matriz resultado da subtração
double** MatrixSubtraction(double **matrix1, double **matrix2){
    double **resultingMatrix = MatrixAllocation(lines, columns);
    for (int l = 0; l < lines ; ++l) {
        for (int m = 0; m < columns; ++m) {
            resultingMatrix[l][m] = matrix1[l][m] - matrix2[l][m];

        }
    }
    return resultingMatrix;
}

/// Função para checar os vetores e identificar se é necessário mudar o sinal
/// \param index
/// \param vectorP
/// \param vectorPLine
/// \param vectorPNormalization
void SignalVerification(int index, const double *vectorP, double *vectorPLine, double vectorPNormalization) {
    if(vectorP[index + 1] > 0){
        vectorPLine[index + 1] = vectorPNormalization * -1.0;
    }

    else {
        vectorPLine[index + 1] = vectorPNormalization;
    }
}

/// Construção da matriz de householder
/// \param matrix
/// \param index
/// \return Matriz de householder
double** ConstructHouseHolder(double **matrix, int index){
    double **houseHolderMatrix = MatrixAllocation(lines, columns);
    double            *vectorP = VectorAllocation();
    double        *vectorPLine = VectorAllocation();

    double               *vectorN;
    double     *vectorNNormalized;
    double   vectorPNormalization;

    InitializeVector(vectorP);
    InitializeVector(vectorPLine);

    MakeIdentityMatrix(houseHolderMatrix);

    for (int i = index + 1; i < vectorSize ; ++i) {
        vectorP[i] = matrix[i][index];
    }

    vectorPNormalization = EuclidianNormalization(vectorP);

    SignalVerification(index, vectorP, vectorPLine, vectorPNormalization);

    vectorN                = VectorSubtraction(vectorP, vectorPLine);
    vectorNNormalized      = NormalizingVector(vectorN);

    houseHolderMatrix = MatrixSubtraction(houseHolderMatrix, MatrixValueMultiplication(
            VectorTranposeVectorMultiplication(vectorNNormalized, vectorNNormalized),2));

    return houseHolderMatrix;
}


/// Calcula a matriz de HouseHolder e a matriz tridiagonal
/// \param matrix
/// \return struct com a Matriz de HouseHolder e a Matriz tridiagonal
struct HouseHolderAnswer HouseHolderMethod(double **matrix){
    double **houseHolderMatrix = MatrixAllocation(lines, columns);
    double **houseHolderMatrixAux;

    struct HouseHolderAnswer houseHolderAnswer;

    MakeIdentityMatrix(houseHolderMatrix);
    for (int i = 0; i < columns - 2 ; ++i) {
        houseHolderMatrixAux = ConstructHouseHolder(matrix, i);
        matrix               = MatrixMultiplication(houseHolderMatrixAux, MatrixMultiplication(matrix,
                                                                                                houseHolderMatrixAux));
        houseHolderMatrix    = MatrixMultiplication(houseHolderMatrix, houseHolderMatrixAux);

    }

    houseHolderAnswer.houseHolderMatrix = houseHolderMatrix;
    houseHolderAnswer.tridiagonalMatrix = matrix;
    return houseHolderAnswer;
}

/// Calcula o maior autovalor pelo metodo da potencia regular
/// \param matrix        = Matriz Original
/// \param initialVector = Vetor chute
/// \param tolerance     = Tolerancia para o erro
/// \return struct com o maior autovalor e o autovetor correspondente
struct EigenValueVector RegularPow(double **matrix, double *initialVector, double tolerance){
    struct EigenValueVector answer;
    answer.eigenVector = VectorAllocation();
    double *q = NormalizingVector(initialVector);
    double error;
    double eigenValue;
    double previousEigenValue = 0;
    double *x = MatrixVectorMultiplication(matrix, q);

    do{
        q = NormalizingVector(x);
        x = MatrixVectorMultiplication(matrix, q);
        eigenValue = VectorMultiplication(q, x);
        error = fabs((eigenValue - previousEigenValue)/eigenValue);
        previousEigenValue = eigenValue;
    }
    while (error > tolerance);
    answer.eigenValue = eigenValue;
    answer.eigenVector = q;
    return answer;
}

/// Função que faz a decomposição LU e encontra a inversa da matriz
/// \param matrix = Matriz original
/// \return ponteiro para a inversa da Matriz original
double** LUDecompositionForInverse(double **matrix){
    double              **l = MatrixAllocation(lines, columns);
    double              **u = MatrixAllocation(lines, columns);
    double              **y = MatrixAllocation(lines, columns);
    double **identityMatrix = MatrixAllocation(lines, columns);
    double  **inverseMatrix = MatrixAllocation(lines, columns);

    //Decomposição LU -> Matriz = l * u
    //Povoamento da Matriz identidade e das Matrizes l e u
    MakeIdentityMatrix(identityMatrix);

    //Decomposição LU
    for(int i = 0; i < lines; i++) {
        for(int j = 0; j < columns; j++) {
            l[i][j] = identityMatrix[i][j];

            if(i <= j) {
                u[i][j] = matrix[i][j];
                for(int k = 0; k < i; k++) {
                    u[i][j] -= l[i][k] * u[k][j];
                }
            }

            else {
                l[i][j] = matrix[i][j];
                for(int k = 0; k < j; k++){
                    l[i][j] -= l[i][k] * u[k][j];
                }
                l[i][j] /= u[j][j];
                u[i][j]  = 0.0;
            }
        }
    }

    //Descobrindo o y da formula para descobrir a Matriz inversa
    // (Matriz * MatrizInversa = Identidade -> l * u * MatrizInversa = Identidade -> y = u * MatrizInversa)
    for(int c = 0; c < columns; c++) {
        for (int i = 0; i < lines; ++i) {
            y[i][c] = identityMatrix[i][c];
            for(int k = 0; k < i; k++) {
                y[i][c] -= l[i][k] * y[k][c];
            }
        }
    }

    //Descobrindo a MatrizInversa
    for(int c = 0; c < columns; c++) {
        for (int i = lines - 1; i >= 0; --i) {
            inverseMatrix[i][c] = y[i][c];
            for(int k = i + 1; k < lines; k++) {
                inverseMatrix[i][c] -= u[i][k] * inverseMatrix[k][c];
            }
            inverseMatrix[i][c] /= u[i][i];
        }
    }
    return inverseMatrix;
}

/// Funçao que calcula o maior autovalor pelo o metodo da potencia inversa
/// \param matrix        = Matriz original
/// \param initialVector = Vetor chute
/// \param tolerance     = Tolerancia para o erro
/// \return struct com o maior autovalor e o autovetor correspondente
struct EigenValueVector InversePow(double **matrix, double *initialVector, double tolerance){
    struct EigenValueVector answer;
    struct EigenValueVector regularPowAnswer;
    double **inverseMatrix = LUDecompositionForInverse(matrix);

    regularPowAnswer   = RegularPow(inverseMatrix, initialVector, tolerance);
    answer.eigenValue  = 1/regularPowAnswer.eigenValue;
    answer.eigenVector = regularPowAnswer.eigenVector;
    return answer;
}

/// Metodo da Potencia com deslocamento
/// \param matrix
/// \param initialVector
/// \param tolerance
/// \param displacement
/// \return uma struct com o autovalor com o deslocamento e o autovetor
struct EigenValueVector DisplacementPow(double **matrix, double *initialVector, double tolerance,
                                        double displacement){
    double **identityMatrix = MatrixAllocation(lines, columns);
    double **displacementMatrix;
    double **displacementIdentityMultiplication;

    struct EigenValueVector answer;
    struct EigenValueVector inversePowAnswer;

    MakeIdentityMatrix(identityMatrix);
    displacementIdentityMultiplication = MatrixValueMultiplication(identityMatrix, displacement);
    displacementMatrix                 = MatrixSubtraction(matrix, displacementIdentityMultiplication);
    inversePowAnswer                   = InversePow(displacementMatrix, initialVector, tolerance);
    answer.eigenValue                  = displacement + inversePowAnswer.eigenValue;
    answer.eigenVector                 = inversePowAnswer.eigenVector;
    return answer;
}

/// Constroi a matriz de jacobi para zerar os valores abaixo da diagonal principal no método QR
/// \param matrix
/// \param index
/// \return Matriz de jacobi
double** ConstructJacobianMatrixTranspose(double **matrix, int index){
    double **jacobianMatrix = MatrixAllocation(lines, columns);
    double      tetaTangent = matrix[index + 1][index]/matrix[index][index];
    double             teta = atan(tetaTangent);

    MakeIdentityMatrix(jacobianMatrix);

    jacobianMatrix [index]     [index]     =  cos(teta);
    jacobianMatrix [index]     [index+1]   =  sin(teta);
    jacobianMatrix [index + 1] [index]     = -sin(teta);
    jacobianMatrix [index + 1] [index + 1] =  cos(teta);

    return jacobianMatrix;
}

/// Faz a decomposição da Matriz em duas matrizes uma Q e outra R
/// \param matrix
/// \return struct com as duas matrizes(Q e R)
struct QRMatrix QRDecomposition(double **matrix){
    struct QRMatrix qrAnswer;
    double **qMatrixTranspose = MatrixAllocation(lines, columns);
    double **qMatrix;
    double **rMatrix;
    double **jacobianMatrixTranspose;

    MakeIdentityMatrix(qMatrixTranspose);
    for (int i = 0; i < lines - 1 ; ++i) {

        jacobianMatrixTranspose = ConstructJacobianMatrixTranspose(matrix, i);
        matrix                  = MatrixMultiplication(jacobianMatrixTranspose, matrix);
        qMatrixTranspose        = MatrixMultiplication(jacobianMatrixTranspose, qMatrixTranspose);

    }
    qMatrix          = MatrixTranspose(qMatrixTranspose);
    rMatrix          = matrix;
    qrAnswer.qMatrix = qMatrix;
    qrAnswer.rMatrix = rMatrix;

    return qrAnswer;
}

/// Calcula auto valor e auto vetor pelo metodo QR
/// \param matrix
/// \param tolerance
/// \return Struct com o autovalor e o autovetor
struct SetOfEigenValueVector QRMethod(double **matrix, double tolerance){
    double **x;
    double **a;
    double **aWithoutDiagonal;
    double **diagonalMatrix;
    double normalization;
    double *eigenValues = VectorAllocation();

    struct HouseHolderAnswer houseHolderAnswer;
    struct QRMatrix qrAnswer;
    struct SetOfEigenValueVector setOfEigenValueVectorAnswer;

    houseHolderAnswer = HouseHolderMethod(matrix);
    x = houseHolderAnswer.houseHolderMatrix;
    a = houseHolderAnswer.tridiagonalMatrix;

    do{
        qrAnswer = QRDecomposition(a);

        a = MatrixMultiplication(qrAnswer.rMatrix, qrAnswer.qMatrix);
        x = MatrixMultiplication(x, qrAnswer.qMatrix);
        diagonalMatrix   = MakeDiagonalMatrix(a);
        aWithoutDiagonal = MatrixSubtraction(a, diagonalMatrix);
        normalization    = MatrixEuclidianNormalization(aWithoutDiagonal);

    }

    while(normalization > tolerance);
    for (int i = 0; i < vectorSize; ++i) {
        eigenValues[i] = diagonalMatrix[i][i];
    }

    setOfEigenValueVectorAnswer.eigenVectors = x;
    setOfEigenValueVectorAnswer.eigenValues  = eigenValues;

    return setOfEigenValueVectorAnswer;
}


int main() {
    struct EigenValueVector           answerRegularPow;
    struct EigenValueVector           answerInversePow;
    struct EigenValueVector      answerDisplacementPow;
    struct SetOfEigenValueVector              answerQr;

    double       **matrix = MatrixAllocation(lines, columns);
    double      tolerance = 0.001;
    double *initialVector = (double*)malloc(vectorSize * sizeof(double*));
    double   displacement = 5.0;

    initialVector[0] = 1.0;
    initialVector[1] = 1.0;
    initialVector[2] = 1.0;
    initialVector[3] = 1.0;

    matrix[0][0] =  4.0;
    matrix[0][1] =  1.0;
    matrix[0][2] = -2.0;
    matrix[0][3] =  2.0;

    matrix[1][0] =  1.0;
    matrix[1][1] =  2.0;
    matrix[1][2] =  0.0;
    matrix[1][3] =  1.0;

    matrix[2][0] = -2.0;
    matrix[2][1] =  0.0;
    matrix[2][2] =  3.0;
    matrix[2][3] = -2.0;

    matrix[3][0] =  2.0;
    matrix[3][1] =  1.0;
    matrix[3][2] = -2.0;
    matrix[3][3] = -1.0;

    answerRegularPow = RegularPow(matrix, initialVector, tolerance);
    printf("Regular Pow eigenValue : %f\n\n", answerRegularPow.eigenValue);
    printf("Regular Pow eigenVector: \n\n");
    PrintVector(answerRegularPow.eigenVector);

    answerInversePow = InversePow(matrix, initialVector, tolerance);
    printf("Inverse Pow eigenValue: %f\n\n", answerInversePow.eigenValue);
    printf("Inverse Pow eigenVector: \n\n");
    PrintVector(answerInversePow.eigenVector);

    answerDisplacementPow = DisplacementPow(matrix, initialVector, tolerance, displacement);
    printf("Displacement Pow eigenValue: %f\n\n", answerDisplacementPow.eigenValue);
    printf("Displacement Pow eigenVector: \n\n");
    PrintVector(answerDisplacementPow.eigenVector);

    printf("HouseHolder matrix:\n");
    PrintMatrix(HouseHolderMethod(matrix).houseHolderMatrix);
    printf("Tridiagonal matrix:\n");
    PrintMatrix(HouseHolderMethod(matrix).tridiagonalMatrix);

    answerQr = QRMethod(matrix, tolerance);
    printf("QrMethod last eigenValue: \n");
    PrintVector(answerQr.eigenValues);
    printf("QrMethod eigenVectors: \n");
    PrintMatrix(answerQr.eigenVectors);

    system("pause");
    return 0;

}
