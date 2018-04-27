//Fernanda Bezerra Nascimento - 388834
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct AutoVectorAutoValue{
    double      Lambda;
    double *Autovector;
};

//Foi utilizada a normalização euclidiana

int Columns    = 3;
int Lines      = 3;
int VectorSize = 3;

double* VectorAlocation(){
    double *Vector = (double*) malloc(VectorSize * sizeof(double));
    return Vector;
}

/// Função para alocação da Matriz
/// \param Lines = quantidade de linhas
/// \param Columns = quantidade de colunas
/// \return a matriz com o espaço alocado
double** MatrixAlocation(int Lines, int Columns){
    double **Matrix = (double**)malloc(Lines * sizeof(double*));

    for (int i = 0; i < Lines; i++){
        Matrix[i] = (double*) malloc(Columns * sizeof(double));
    }
    return Matrix;
}

/// Imprime a Matriz
/// \param Matrix
void PrintMatrix(double **Matrix){
    for (int l = 0; l < Lines; ++l) {
        for (int m = 0; m < Columns; ++m) {
            printf("%f", Matrix[l][m]);
            printf("\t");
        }
        printf("\n");
    }
}

/// Normalização do vetor
/// \param Vector
/// \return a norma do vetor
double EuclidianNormalization(double *Vector) {
    double Normalization = 0;
    for (int j = 0; j < VectorSize; ++j) {
        Normalization += pow(Vector[j], 2);
    }
    Normalization = sqrt(Normalization);
    return Normalization;
}

/// Multiplica matriz e vetor
/// \param Matrix
/// \param Vector
/// \return o vetor resultante da multiplicação
double* MatrixVectorMultiplication(double **Matrix, const double *Vector){
    double *ResultingVector = (double*) malloc(VectorSize * sizeof(double));
    for (int i = 0; i < Lines ; ++i) {
        ResultingVector[i] = 0;
        for (int j = 0; j < Columns; ++j) {
            ResultingVector[i] += Vector[j] * Matrix[i][j];
        }
    }
    return ResultingVector;
}

/// Normaliza o vetor
/// \param Vector
/// \return o vetor normalizado
double* NormalizingVector(double *Vector) {
    double *NormalizedVector = (double*) malloc(VectorSize * sizeof(double));
    double     Normalization = EuclidianNormalization(Vector);
    for (int i = 0; i < VectorSize; ++i) {
        NormalizedVector[i] = Vector[i]/Normalization;
    }
    return NormalizedVector;
}

/// Multiplicação de vetores
/// \param Vector1
/// \param Vector2
/// \return o resultado da multiplicacao
double VectorMultiplication(const double *Vector1, const double *Vector2) {
    double Result = 0;
    for (int i = 0; i < VectorSize ; ++i) {
        Result += Vector1[i] * Vector2[i];
    }
    return Result;
}

/// Calcula o maior autovalor pelo metodo da potencia regular
/// \param Matrix        = Matriz Original
/// \param InitialVector = Vetor chute
/// \param Tolerance     = Tolerancia para o erro
/// \return struct com o maior autovalor e o autovetor correspondente
struct AutoVectorAutoValue RegularPow(double **Matrix, double *InitialVector, double Tolerance){
    struct AutoVectorAutoValue Answer;
    Answer.Autovector = VectorAlocation();
    double *q = NormalizingVector(InitialVector);
    double Error;
    double Lambda;
    double PreviousLambda = 0;
    double *x = MatrixVectorMultiplication(Matrix, q);

    do{
        q = NormalizingVector(x);
        x = MatrixVectorMultiplication(Matrix, q);
        Lambda = VectorMultiplication(q, x);
        Error = fabs((Lambda - PreviousLambda)/Lambda);
        PreviousLambda = Lambda;
    }
    while (Error > Tolerance);
    Answer.Lambda = Lambda;
    Answer.Autovector = x;
    return Answer;
}

/// Função que faz a decomposição LU e encontra a inversa da matriz
/// \param Matrix = Matriz original
/// \return ponteiro para a inversa da Matriz original
double** LUDecompositionForInverse(double **Matrix){
    double              **L = MatrixAlocation(Lines, Columns);
    double              **U = MatrixAlocation(Lines, Columns);
    double              **Y = MatrixAlocation(Lines, Columns);
    double **IdentityMatrix = MatrixAlocation(Lines, Columns);
    double  **InverseMatrix = MatrixAlocation(Lines, Columns);

    //Decomposição LU -> Matriz = L * U
    //Povoamento da Matriz identidade e das Matrizes L e U
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

    //Decomposição LU
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

    //Descobrindo o Y da formula para descobrir a Matriz inversa
    // (Matriz * MatrizInversa = Identidade -> L * U * MatrizInversa = Identidade -> Y = U * MatrizInversa)
    for(int c = 0; c < Columns; c++) {
        for (int i = 0; i < Lines; ++i) {
            Y[i][c] = IdentityMatrix[i][c];
            for(int k = 0; k < i; k++) {
                Y[i][c] -= L[i][k] * Y[k][c];
            }
        }
    }

    //Descobrindo a MatrizInversa
    for(int c = 0; c < Columns; c++) {
        for (int i = Lines - 1; i >= 0; --i) {
            InverseMatrix[i][c] = Y[i][c];
            for(int k = i + 1; k < Lines; k++) {
                InverseMatrix[i][c] -= U[i][k] * InverseMatrix[k][c];
            }
            InverseMatrix[i][c] /= U[i][i];
        }
    }
    return InverseMatrix;
}

/// Funçao que calcula o maior autovalor pelo o metodo da potencia inversa
/// \param Matrix        = Matriz original
/// \param InitialVector = Vetor chute
/// \param Tolerance     = Tolerancia para o erro
/// \return struct com o maior autovalor e o autovetor correspondente
struct AutoVectorAutoValue InversePow(double **Matrix, double *InitialVector, double Tolerance){
    struct AutoVectorAutoValue Answer;
    struct AutoVectorAutoValue RegularPowAnswer;
    double **InverseMatrix = LUDecompositionForInverse(Matrix);

    RegularPowAnswer  = RegularPow(InverseMatrix, InitialVector, Tolerance);
    Answer.Lambda     = 1/RegularPowAnswer.Lambda;
    Answer.Autovector = RegularPowAnswer.Autovector;
    return Answer;
}

/// Função para calcular a matriz resultante da multiplicação de uma matriz por um valor escalar
/// \param Matrix
/// \param Displacement
/// \return Matriz resultado da multiplicação
double** MatrixValueMultiplication(double **Matrix, double Displacement){
    double **ResultingMatrix = MatrixAlocation(Lines, Columns);
    for (int l = 0; l < Lines ; ++l) {
        for (int m = 0; m < Columns; ++m) {
            ResultingMatrix[l][m] = Matrix[l][m] * Displacement;

        }
    }
    return ResultingMatrix;
}

/// Função que calcula a matrix resultante da subtração de duas matrizes
/// \param Matrix1
/// \param Matrix2
/// \return Matriz resultado da subtração
double** MatrixSubtraction(double **Matrix1, double **Matrix2){
    double **ResultingMatrix = MatrixAlocation(Lines, Columns);
    for (int l = 0; l < Lines ; ++l) {
        for (int m = 0; m < Columns; ++m) {
            ResultingMatrix[l][m] = Matrix1[l][m] - Matrix2[l][m];

        }
    }
    return ResultingMatrix;
}

/// Metodo da Potencia com deslocamento
/// \param Matrix
/// \param InitialVector
/// \param Tolerance
/// \param Displacement
/// \return uma struct com o autovalor com o deslocamento e o autovetor
struct AutoVectorAutoValue DisplacementPow(double **Matrix, double *InitialVector, double Tolerance,
                         double Displacement){
    double **IdentityMatrix = MatrixAlocation(Lines, Columns);
    struct AutoVectorAutoValue Answer;
    struct AutoVectorAutoValue InversePowAnswer;
    double **DisplacementMatrix;
    double **DisplacementIdentityMultiplication;

    for (int l = 0; l < Lines ; ++l) {
        for (int m = 0; m < Columns; ++m) {
            if (l == m){
                IdentityMatrix[l][m] = 1;
            }
            else{
                IdentityMatrix[l][m] = 0;
            }
        }
    }
    DisplacementIdentityMultiplication = MatrixValueMultiplication(IdentityMatrix, Displacement);
    DisplacementMatrix                 = MatrixSubtraction(Matrix, DisplacementIdentityMultiplication);
    InversePowAnswer                   = RegularPow(DisplacementMatrix, InitialVector, Tolerance);
    Answer.Lambda                      = Displacement + InversePowAnswer.Lambda;
    Answer.Autovector                  = InversePowAnswer.Autovector;
    return Answer;
}

int main() {
    struct AutoVectorAutoValue AnswerRegularPow;
    struct AutoVectorAutoValue AnswerInversePow;

    double       **Matrix = MatrixAlocation(Lines, Columns);
    double      Tolerance = 0.0000001;
    double *InitialVector = (double*)malloc(VectorSize * sizeof(double*));

    InitialVector[0] = 1;
    InitialVector[1] = 1;
    InitialVector[2] = 1;
    Matrix[0][0] = 4;
    Matrix[0][1] = -1;
    Matrix[0][2] = 1;
    Matrix[1][0] = -1;
    Matrix[1][1] = 3;
    Matrix[1][2] = -2;
    Matrix[2][0] = 1;
    Matrix[2][1] = -2;
    Matrix[2][2] = 3;

    AnswerRegularPow = RegularPow(Matrix, InitialVector, Tolerance);
    printf("%f\n", AnswerRegularPow.Lambda);
    AnswerInversePow = InversePow(Matrix, InitialVector, Tolerance);
    printf("%f\n", AnswerInversePow.Lambda);

    system("pause");
    return 0;
}