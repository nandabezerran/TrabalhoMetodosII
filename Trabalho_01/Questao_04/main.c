//FERNANDA BEZERRA NASCIMENTO - QUARTA QUESTÃO - MATRÍCULA: 388834 - MÉTODOS NUMÉRICOS II
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

const double ncPesos[8][5] = {
                {1, 1,    0,    0,  0}, //Newton Cotes pesos fechados
                {1, 4,    1,    0,  0},
                {1, 3,    3,    1,  0},
                {7, 32,  12,   32,  7},
                {1, 1,    0,    0,  0}, //Newton Cotes pesos abertos
                {2, -1,   2,    0,  0},
                {11, 1,   1,   11,  0},
                {11, -14, 26, -14, 11}
};

const double glRaizes[4][5] = {
        {- 0.57735, 0.57735,     0,      0,    0},
        {0, -0.77459, 0.774596,          0,    0},
        {-0.3399, 0.33998, -0.8611, 0.86113,   0},
        {0, -0.53846, 0.53846, -0.90617, 0.90617}
};

const double ncMultiplicador[8] = {1.0/2.0, 1.0/3.0, 3.0/8.0, 2.0/45.0, 3.0/2.0, 4.0/3.0, 5.0/24.0, 6.0/20.0};

const double glPesos[4][5] = {
        {1.0, 1.0,        0,             0,         0},
        {0.88888, 0.55555, 0.55555,      0,         0},
        {0.65214, 0.65214, 0.34785, 0.34785,        0},
        {0.56888, 0.478626, 0.47862, 0.23692, 0.23692}

};

/// Função para calcular a função inserida
/// \param x
/// \param y
/// \return
double fun(double x, double y) {
    return x*y + (pow(x, 2.0) * pow(y, 3.0));
}

/// Função para efetuar a mudança de variável, no caso calculariámos X(alfa)
/// \param x = o nosso alfa
/// \param a = o nosso limite inferior do intervalo
/// \param b = o nosso limite superior do intervalo
/// \return x(alfa)
double funXGl(const double x, double a, double b){
    return (((a + b)/2) + (x * ((a - b)/2)));
}

/// Calcular a integral dupla
/// \param a
/// \param b
/// \param c
/// \param d
/// \param metodoX
/// \param metodoY
/// \return
double calculoIntegralDupla(double a, double b,
                            double c, double d, int metodoX, int metodoY) {
    int dxGrau=0, dyGrau=0;

    //Pegar o grau correto de acordo com o método entrado
    dxGrau = ((metodoX - 1)%4) + 1;
    dyGrau = ((metodoY - 1)%4) + 1;

    int dxPointsNum = dxGrau + 1;
    double pointsDx[dxPointsNum];
    double dxH;

    //Fazendo a altura de acordo com NC Fechado e Aberto e Gauss Legendre
    dxH = (b - a) / (dxGrau + ((metodoX > 4) ? 2.0 : 0));
    for (int i = 0; i < dxPointsNum; ++i) {
        if (metodoX <= 4) {
            pointsDx[i] = ((dxH) * i) + a;
        }
        else if (metodoX <= 8) {
            pointsDx[i] = ((dxH) * i) + a + dxH;
        }
        else{
            pointsDx[i] = funXGl(glRaizes[dxGrau - 1][i], a, b);
        }
    }

    int dyPointsNum = dyGrau + 1;
    double pointsDy[dyPointsNum];
    double dyH;

    //Fazendo a altura de acordo com NC Fechado e Aberto e Gauss Legendre
    if (metodoY > 4){
        dyH = (d - c) / (dyGrau + 2.0);
    }
    else if (metodoY <= 8){
        dyH = (d - c) / dyGrau;
    }

    for (int j = 0; j < dyPointsNum; ++j) {
        if (metodoY <= 4) {
            pointsDy[j] = (dyH * j) + c;
        }
        else if (metodoY <= 8) {
            pointsDy[j] = ((dyH) * j) + c + dyH;
        }
        else{
            pointsDy[j] = funXGl(glRaizes[dyGrau - 1][j], c, d);
        }
    }

    double result = 0.0;
    double weightX, weightY;
    for (int k = 0; k < dxPointsNum; ++k) {
        for (int l = 0; l < dyPointsNum; ++l) {
            if (metodoX <= 8) {
                weightX = ncPesos[metodoX - 1][k];
            }
            if (metodoX > 8) {
                weightX = glPesos[dxGrau - 1][k];
            }
            if (metodoY <= 8) {
                weightY = ncPesos[metodoY - 1][l];
            }
            if (metodoY > 8) {
                weightY = glPesos[dyGrau - 1][l];
            }

            double weight = weightX * weightY;
            result += (weight) * (fun(pointsDx[k], pointsDy[l]));
        }
    }
    double multDx;
    double multDy;

    if (metodoX > 8 ){
        multDx = ((b-a)/2);
    }
    else {
        multDx = dxH * ncMultiplicador[metodoX - 1];
    }

    if(metodoY > 8){
        multDy = ((d-c)/2);
    }
    else{
        multDy = dyH * ncMultiplicador[metodoY - 1];
    }
    result *= multDx * multDy;

    return result;
}

/// Calcular a convergencia da integral dupla
/// \param botY
/// \param topY
/// \param botX
/// \param topX
/// \param toler
/// \param metodoX
/// \param metodoY
/// \return
double convergenciaIntegralDupla(double botY, double topY, double botX, double topX, double toler, int metodoX,
                                 int metodoY) {
    double partNum = 1.0;
    double error;
    double prevResult;
    double result = 0.0;
    do {

        prevResult = result;
        result = 0.0;
        double partXSize = (topX - botX) / partNum;
        double partYSize = (topY - botY) / partNum;
        for (int i = 0; i < partNum; ++i) {
            for (int j = 0; j < partNum; ++j) {

                double a = botX + (i * partXSize);
                double b = a + partXSize;
                double c = botY + (j * partYSize);
                double d = c + partYSize;

                result += calculoIntegralDupla(a, b, c, d, metodoX, metodoY);
            }
        }
        partNum *= 2.0;
        error = fabs(result - prevResult) / result;
    } while (error > toler || partNum <= .0);

    return result;
}

int main() {

    double botX;
    double topX;
    double botY;
    double topY;
    double toler;

    printf("Entre com o intervalo da integral de Y\n");
    printf("g1(x):>>> ");
    //scanf("%f",&botY);
    botY = 0;
    printf("g2(x):>>> ");
    scanf("%f",&topY);

    printf("Entre com o intervalo da integral de X\n");
    printf("a :>>> ");
    //scanf("%f",&botX);
    botX = 0;
    printf("b :>>> ");
    scanf("%f",&topX);
    topX = M_PI/2;


    printf("Entre com a tolerância da integral de X\n");
    scanf("%f",&toler);
    int metodoX, metodoY;

    printf("\n--------Métodos:---------\n");
    printf("1 -  Newton Cotes Fechada Primeiro Grau\n2 -  Newton Cotes Fechada Segundo Grau\n3 -  Newton Cotes Fechada Terceiro Grau\n"
                   "4 -  Newton Cotes Fechada Quarto Grau\n5 -  Newton Cotes Aberta Primeiro Grau\n6 -  Newton Cotes Aberta Segundo Grau\n"
                   "7 -  Newton Cotes Aberta Terceiro Grau\n8 -  Newton Cotes Aberta Quarto Grau\n9 -  Gauss Legendre Grau 2\n10 - Gauss Legendre Grau 3\n"
                   "11 - Gauss Legendre Grau 4\n12 - Gauss Legendre Grau 5\n>>> ");
    printf("Metodo para integração em X:\n>>> ");
    scanf("%d", &metodoX);
    printf("Metodo para integração em Y:\n>>> ");
    scanf("%d", &metodoY);

    double resultado = convergenciaIntegralDupla(botY, topY, botX, topX, toler, metodoX, metodoY);
    printf("Método para dx:%d Método para dy:%d\nIntegral:\n>>> %f\n\n", metodoX, metodoY, resultado);

    return 0;
    system("pause");
}