//Fernanda Bezerra Nascimento - 388834
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Definição das stucts para os metodos
typedef struct node node;

typedef struct element element;

struct node{
    int incidence;
    double coordinate;
};

struct element{
    node node[2];
    double **kMatrix;
    double *fVector;
    double *gVector;
};


/// Função para alocação do espaço de um vetor
/// \param vectorSize
/// \return vetor com o espaço alocado
double* VectorAllocation(int vectorSize){
    double *vector = (double*) malloc(vectorSize * sizeof(double));
    return vector;
}

/// Função para inicializar o vetor com 0
/// \param vector, vectorSize
void InitializeVector(double *vector,int vectorSize) {
    for (int l = 0; l < vectorSize ; ++l) {
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

/// Função para inicializar a matriz com 0
/// \param matriz, lines, columns
void InitializeMatrix(double **matrix, int lines, int columns){
    for (int i = 0; i < lines ; ++i) {
        for (int j = 0; j < columns; ++j) {
            matrix[i][j] = 0;
        }
    }
}
/// Imprime a Matriz
/// \param matrix
void PrintMatrix(double **matrix, int lines, int columns){
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
void PrintVector(double *vector, int lines, char vectorName){
    for (int i = 0; i < lines; ++i) {
        printf("%c%d: ",vectorName,i);
        printf("%f", vector[i]);
        printf("\n");
    }
    printf("\n");
}

/// Transforma uma matriz numa matriz identidade
/// \param matrix
void MakeIdentityMatrix(double **matrix, int lines, int columns) {
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

/// Função que faz a decomposição LU e encontra a inversa da matriz
/// \param matrix = Matriz original
/// \return ponteiro para a inversa da Matriz original
double* LUMethodForLinearSystems(double **matrix, double const *vector, int lines, int columns){
    double              **l = MatrixAllocation(lines, columns);
    double              **u = MatrixAllocation(lines, columns);
    double               *y = VectorAllocation(lines);
    double         *results = VectorAllocation(lines);
    double **identityMatrix = MatrixAllocation(lines, columns);
    double  **inverseMatrix = MatrixAllocation(lines, columns);

    //Decomposição LU -> Matriz = l * u
    //Povoamento da Matriz identidade e das Matrizes l e u
    MakeIdentityMatrix(identityMatrix, lines, columns);

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
    for (int i = 0; i < lines; ++i) {
        y[i] = vector[i];
        for(int k = 0; k < i; k++) {
            y[i] -= l[i][k] * y[k];
        }
    }

    for(int i=lines-1; i>=0; i--) {
        results[i]= y[i];
        for(int j=i+1; j< lines; j++) {
            results[i]-=u[i][j]*results[j];
        }
        results[i]/=u[i][i];
    }
    return results;
    }

/// Calcula a matriz A do metodo de elementos finitos
/// \param elements
/// \param matrixSize
/// \param numberOfPartitions
/// \return matriz A calculada
double** MatrixAFiniteElements(element *elements, int matrixSize, int numberOfPartitions){
    double **matrixA = MatrixAllocation(matrixSize, matrixSize);
    InitializeMatrix(matrixA, matrixSize, matrixSize);
    int numberOfNodes = 2; //Cada elemento sempre é ligado em 2 nós
    int aI, aJ;
    for(int position = 0; position < numberOfPartitions; position++){
        for(int i = 0; i < numberOfNodes; i++){
            for(int j = 0; j < numberOfNodes; j++){
                aI = elements[position].node[i].incidence;
                aJ = elements[position].node[j].incidence;
                if(aI != -1 && aJ != -1){
                    matrixA[aI][aJ] += elements[position].kMatrix[i][j];
                }
            }
        }
    }
    //matrixA[matrixSize-1][matrixSize-1] += 2.0;
    return matrixA;
}

/// Calcula o vetor B do método de elementos finitos (com as contribuições de f e g)
/// \param elements
/// \param numberOfPartitions
/// \param vectorSize
/// \return o vetor B
double* VectorBFiniteElements(element *elements, int numberOfPartitions, int vectorSize){
    double *vectorB = VectorAllocation(vectorSize);
    InitializeVector(vectorB, vectorSize);
    int numberOfNodes = 2;
    int bI;
    for(int position = 0; position < numberOfPartitions; position++){
        for(int i = 0; i < numberOfNodes; i++){
            bI = elements[position].node[i].incidence;
            if(bI != -1){
                vectorB[bI] += (elements[position].fVector[i] - elements[position].gVector[i]);

            }
        }
    }

    return vectorB;
}

/// Calcula o PVC pelo metodo dos elementos finitos
/// \param matrixASize = Tamanho da matriz A
/// \param numberOfPartitions = numero de partições da discretização do domínio
/// \param vectorBSize = tamanho do vetor B
/// \param elements = vetor de structs de elementos que são as partições
void FiniteElementsMethod(int matrixASize, int numberOfPartitions, int vectorBSize, element *elements) {
    double** matrixA =  MatrixAFiniteElements(elements, matrixASize, numberOfPartitions);
    double*  vectorB = VectorBFiniteElements(elements, numberOfPartitions,vectorBSize);
    printf("              ++++++++ METODO DOS ELEMENTOS FINITOS ++++++++\n");
    printf("Matriz A:\n");
    PrintMatrix(matrixA, matrixASize,matrixASize);
    printf("Vetor B:\n");
    PrintVector(vectorB,vectorBSize, 'B');
    printf("Resultado do Sistema A * y = B:\n");
    PrintVector(LUMethodForLinearSystems(matrixA, vectorB, matrixASize, matrixASize),matrixASize, 'y');

}

int main() {

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////                             MONTAGEM DOS ELEMENTOS PARA A EXECUÇÃO DO MÉTODO                            /////
    /////                                        DE ELEMENTOS FINITOS                                             /////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int        matrixASize = 5;
    int numberOfPartitions = 6;
    int        vectorBSize = 5;
    int      sizeOfKMatrix = 2;
    int     sizeOfFGVector = 2;

    double lenghOfInterval = 0.3;
    double          deltaX = lenghOfInterval / numberOfPartitions;

    element      *elements = malloc(numberOfPartitions * sizeof(element));

    elements[0].node[0].incidence  = -1;
    elements[0].node[1].incidence  = 0;
    elements[0].node[0].coordinate = 0.2;
    elements[0].node[1].coordinate = 0.25;

    elements[1].node[0].incidence  = 0;
    elements[1].node[1].incidence  = 1;
    elements[1].node[0].coordinate = 0.25;
    elements[1].node[1].coordinate = 0.30;

    elements[2].node[0].incidence  = 1;
    elements[2].node[1].incidence  = 2;
    elements[2].node[0].coordinate = 0.30;
    elements[2].node[1].coordinate = 0.35;

    elements[3].node[0].incidence  = 2;
    elements[3].node[1].incidence  = 3;
    elements[3].node[0].coordinate = 0.35;
    elements[3].node[1].coordinate = 0.40;

    elements[4].node[0].incidence  = 3;
    elements[4].node[1].incidence  = 4;
    elements[4].node[0].coordinate = 0.40;
    elements[4].node[1].coordinate = 0.45;

    elements[5].node[0].incidence  = 4;
    elements[5].node[1].incidence  = -1;
    elements[5].node[0].coordinate = 0.45;
    elements[5].node[1].coordinate = 0.50;

    double le;
    double coordinateJ;
    double coordinateI;
    double sumOfCoordinates;

    for(int i = 0; i < numberOfPartitions; i++){
        elements[i].kMatrix = MatrixAllocation(sizeOfKMatrix,sizeOfKMatrix);
        elements[i].gVector = VectorAllocation(sizeOfFGVector);
        elements[i].fVector = VectorAllocation(sizeOfKMatrix);

        coordinateI      = elements[i].node[0].coordinate;
        coordinateJ      = elements[i].node[1].coordinate;
        sumOfCoordinates = coordinateI + coordinateJ;
        le = elements[i].node[1].coordinate - elements[i].node[0].coordinate;

        elements[i].kMatrix[0][0] = (- 19 * le * log(sumOfCoordinates + le) + 19 * le * log(sumOfCoordinates - le)
                                     - 19 * (sumOfCoordinates) * log(sumOfCoordinates + le)
                                     + 19 * (sumOfCoordinates) * log(sumOfCoordinates - le))
                                    / (2 * pow(le,2));

        elements[i].kMatrix[0][1] = (19 *(sumOfCoordinates + le) * (log(sumOfCoordinates + le) -
                                    log(sumOfCoordinates - le))) / (2 * pow(le,2));

        elements[i].kMatrix[1][0] = (19 *(sumOfCoordinates - le) * (log(sumOfCoordinates + le) -
                                    log(sumOfCoordinates - le))) / (2 * pow(le,2));

        elements[i].kMatrix[1][1] =  (19 * le*log(sumOfCoordinates + le) - 19 * le * log( sumOfCoordinates - le)
                                      - 19 * sumOfCoordinates * log(sumOfCoordinates +le) +
                                      19 * sumOfCoordinates*log(sumOfCoordinates - le)) / (2 * pow(le, 2));

        elements[i].gVector[0]    = 0;
        elements[i].gVector[1]    = 0;


        elements[i].fVector[0]    = -15 *(le/2);
        elements[i].fVector[1]    = -15 *(le/2);

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////                                          TÉRMINO DA MONTAGEM                                            /////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Execução do método dos elementos finitos
    FiniteElementsMethod(matrixASize, numberOfPartitions, vectorBSize, elements);

    system("pause");
    return 0;
}
