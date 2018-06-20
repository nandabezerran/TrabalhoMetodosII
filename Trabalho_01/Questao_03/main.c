//FERNANDA BEZERRA NASCIMENTO - TERCEIRA QUESTÃO - MATRÍCULA: 388834 - MÉTODOS NUMÉRICOS II
#include <stdio.h>
#include <tgmath.h>
#include <stdlib.h>

/// Função para calcular a função inserida
/// \param x
/// \return resultado da função
double fun (double x){
    //Para uma função diferente para integral, edite a linha abaixo
    return (double)cos(x);
    //return (double)(1.0/pow(x, 2.0/3.0));
}

/// Função para a mudança de variável para exponecial simples
/// \param a = o nosso limite inferior do intervalo
/// \param b = o nosso limite superior do intervalo
/// \param alfa
/// \return x(alfa)
double x_ExpSimples(double a, double b, double alfa){
    return (double)(((a+b)/2.0) + (((b - a)/2.0) * tanh(alfa)));
}

/// Função para a mudança de variável para exponecial dupla
/// \param a = o nosso limite inferior do intervalo
/// \param b = o nosso limite superior do intervalo
/// \param alfa
/// \return x(alfa)
double x_ExpDupla(double a, double b, double alfa){
    return (double)(((a+b)/2.0) + (((b - a)/2.0) * tanh((M_PI/2.0) * sinh(alfa))));
}

/// Cálculo da exponencial simples
/// \param a = o nosso limite inferior do intervalo
/// \param b = o nosso limite superior do intervalo
/// \param alfa
/// \return resultado da exponencial
double g_ExpSimples(double a, double b, double alfa){
    double expSimples = x_ExpSimples(a, b, alfa);
    double fun1 = fun(expSimples);
    return (double) (fun1 * ((b - a) / pow(cosh(alfa), 2.0)));
}

/// Calculo de u da midança da de variavel na derivada
/// \param alfa
/// \return
double u_Alfa(double alfa){
    return (double)((M_PI/2) * sinh(alfa));
}

/// Calculo da derivada
/// \param a
/// \param b
/// \param alfa
/// \return
double d_xAlfa(double a, double b, double alfa){
    return (double)(((M_PI * (b - a))/4.0) * (cosh(alfa)/pow(cosh(u_Alfa(alfa)),2.0)));
}

/// Cálculo da exponencial dupla
/// \param a = o nosso limite inferior do intervalo
/// \param b = o nosso limite superior do intervalo
/// \param alfa
/// \return resultado da exponencial
double g_ExpDupla(double a, double b, double alfa){
    return fun(x_ExpDupla(a, b, alfa)) * d_xAlfa(a, b, alfa);
}

/// Função para escolha do método da exponencial
/// \param a = o nosso limite inferior do intervalo
/// \param b = o nosso limite superior do intervalo
/// \param alfa
/// \param aux
/// \return resultado da exponencial
double Exponencial(double a, double b, double alfa, int aux){
    if (aux == 1){
        return g_ExpSimples(a, b, alfa);
    }
    else if (aux == 2){
        return g_ExpDupla(a, b, alfa);
    }
}
/// Função para efetuar a mudança de variável, no caso calculariámos X(alfa)
/// \param x = o nosso alfa
/// \param a = o nosso limite inferior do intervalo
/// \param b = o nosso limite superior do intervalo
/// \return x(alfa)
double fun_x(double x, double a, double b){
    return (((a + b)/2.0) + (x * ((a - b)/2.0)));
}

/// Calula a distancia entre cada ponto no nosso intervalo
/// \param a = o nosso limite inferior do intervalo
/// \param b = o nosso limite superior do intervalo
/// \param num_part = numero de partições
/// \return delta_x que é a distância entre cada ponto
double delta_x (double a, double b, int num_part){
    return (double)fabs((b-a)/num_part);
}

/// Essa função utilizamos em todos os graus do gauss legendre
/// \param a          = o nosso limite inferior do intervalo
/// \param b          = o nosso limite superior do intervalo
/// \param c1         = o nosso limite inferior do intervalo
/// \param c2         = o nosso limite superior do intervalo
/// \param toler      = tolerancia (caso utilizarmos a estratégia 1)
/// \param num_part_u = numero de partições fixas (caso utilizarmos a estratégia 2)
/// \param estrat     = estratégia que será utilizada
/// \param w          = vetor com os pesos
/// \param x          = vetor com as raizes
/// \param grau       = grau do gauss legendre
/// \return resultado da integral
double intr_gauss(double a, double b, double c1, double c2, double toler, const double *w, const double *x, double grau, int num_part,
                  int n, int aux) {
    double xi, intr = 0, intr_aux, xf, sum, erro_r, deltaX;
        do {
            intr_aux = intr;
            intr = 0;
            deltaX = delta_x(c1, c2, num_part);
            for (n = 0; n < num_part; n++) {
                xi = c1 + (n * deltaX);
                xf = xi + deltaX;
                sum = 0;
                for (int i = 0; i < grau; i++) {
                    sum += (w[i] * (Exponencial(a, b, fun_x(x[i], xi, xf), aux)));
                }
                intr += ((xf - xi) / 2.0) * sum;
            }

            num_part *= 2;
            erro_r = (double)fabs((intr - intr_aux) / intr);

        } while (erro_r > toler);


    return intr;
}

double trapezoidal_fechada(double a, double b, double c1, double c2, double toler,int aux){
    double xi, intr = 0, xf, intr_aux, h;
    int n = 0;
    double erro_r;
    double deltaX;
    int num_part = 1;
    do {
        intr_aux = intr;
        intr = 0;
        deltaX = delta_x(c1, c2, num_part);
        h = deltaX ;
        for (n = 0; n < num_part; n++) {
            xi = c1 + (n * deltaX);
            xf = xi + deltaX;
            intr += ((h/2)* (Exponencial(a, b, xi, aux) + Exponencial(a, b, xf, aux)));

        }
        num_part *= 2;
        erro_r = (double)fabs((intr - intr_aux) / intr);


    } while (erro_r > toler);
    return intr;
}

double trapezoidal_aberta(double a, double b, double c1, double c2, double toler, int aux) {
    double xi, intr = 0, intr_aux, xm_1, xm_2, h;
    int n = 0;
    double erro_r;
    double deltaX;
    int num_part = 1;
    do {
        intr_aux = intr;
        intr = 0;
        deltaX = delta_x(c1, c2, num_part);
        h = deltaX / 3;
        for (n = 0; n < num_part; n++) {
            xi = c1 + (n * deltaX);
            xm_1 = xi + h;
            xm_2 = xm_1 + h;
            intr += ((3 * h) / 2) * (Exponencial(a, b, xm_1, aux) + Exponencial(a, b, xm_2, aux));
        }
        num_part *= 2;
        erro_r = (double)fabs((intr - intr_aux) / intr);
    } while (erro_r > toler);

    return intr;
}

double NewtonCotes_2_fechada (double a, double b, double c1, double c2, double toler, int aux){
    double xi, intr = 0, xf, intr_aux, xm, h;
    int num_part = 1;
    int n = 0;
    double erro_r;
    double deltaX;
    do {
        intr_aux = intr;
        intr = 0;
        deltaX = delta_x(c1, c2, num_part);
        h = deltaX / 2;
        for (n = 0; n < num_part; n++) {
            xi = c1 + (n * deltaX);
            xf = xi + deltaX;
            xm = xi + h;
            intr += (h / 3) * ((Exponencial(a, b, xi, aux)) + 4 * (Exponencial(a, b, xm, aux)) + (Exponencial(a, b, xf, aux)));
        }
        num_part *= 2;
        erro_r = (double)fabs((intr - intr_aux) / intr);

    } while (erro_r > toler);

    return intr;

}

double NewtonCotes_3_fechada (double a, double b, double c1, double c2, double toler, int aux){
    double xi, intr = 0, xf, intr_aux, xm_1, xm_2, h;
    int num_part = 1;
    int n = 0;
    double erro_r;
    double deltaX;
    do {
        intr_aux = intr;
        intr = 0;
        deltaX = delta_x(c1, c2, num_part);
        h = deltaX / 3;
        for (n = 0; n < num_part; n++) {
            xi = c1 + (n * deltaX);
            xf = xi + deltaX;
            xm_1 = xi + h;
            xm_2 = xm_1 + h;
            intr += ((3 * h) / 8) *
                    ((Exponencial(a, b, xi, aux)) + 3 * (Exponencial(a, b, xm_1, aux)) + 3 * (Exponencial(a, b, xm_2, aux)) + (Exponencial(a, b, xf, aux)));
            }
        num_part *= 2;
        erro_r = (double)fabs((intr - intr_aux) / intr);

    } while (erro_r > toler);

    return intr;

}

double NewtonCotes_4_fechada (double a, double b, double c1, double c2, double toler, int aux){
    double xi, intr = 0, xf, intr_aux, xm_1, xm_2, xm_3, h;
    int num_part = 1;
    int n = 0;
    double erro_r;
    double deltaX;
    do {
        intr_aux = intr;
        intr = 0;
        deltaX = delta_x(c1, c2, num_part);
        h = deltaX / 4;
        for (n = 0; n < num_part; n++) {
            xi = c1 + (n * deltaX);
            xf = xi + deltaX;
            xm_1 = xi + h;
            xm_2 = xm_1 + h;
            xm_3 = xm_2 + h;
            intr += ((2 * h) / 45) *
                    (7 * (Exponencial(a, b, xi, aux)) + 32 * (Exponencial(a, b, xm_1, aux)) + 12 * (Exponencial(a, b, xm_2, aux)) +
                         32 * (Exponencial(a, b, xm_3, aux)) + 7 * (Exponencial(a, b, xf, aux)));
        }
        num_part *= 2;
        erro_r = (double)fabs((intr - intr_aux) / intr);
    } while (erro_r > toler);

    return intr;

}

double NewtonCotes_2_aberta (double a, double b, double c1, double c2, double toler, int aux){
    double xi, intr = 0, intr_aux, xm_1, xm_2, xm_3, h;
    int num_part = 1;
    int n = 0;
    double erro_r;
    double deltaX;
    do {
        intr_aux = intr;
        intr = 0;
        deltaX = delta_x(c1, c2, num_part);
        h = deltaX / 4;
        for (n = 0; n < num_part; n++) {
            xi = c1 + (n * deltaX);
            xm_1 = xi + h;
            xm_2 = xm_1 + h;
            xm_3 = xm_2 + h;
            intr += ((4 * h) / 3) * (2 * (Exponencial(a, b, xm_1, aux)) - (Exponencial(a, b, xm_2, aux)) + 2 * (Exponencial(a, b, xm_3, aux)));
        }
        num_part *= 2;
        erro_r = (double)fabs((intr - intr_aux) / intr);
    } while (erro_r > toler);


    return intr;

}

double NewtonCotes_3_aberta (double a, double b, double c1, double c2, double toler, int aux){
    double xi, intr = 0, intr_aux, xm_1, xm_2, xm_3, xm_4, h, xf;
    int num_part = 1;
    int n = 0;
    double erro_r;
    double deltaX;
    do {
        intr_aux = intr;
        intr = 0;
        deltaX = delta_x(c1, c2, num_part);
        h = deltaX / 5;
        for (n = 0; n < num_part; n++) {
            xi = c1 + (n * deltaX);
            xm_1 = xi + h;
            xm_2 = xm_1 + h;
            xm_3 = xm_2 + h;
            xm_4 = xm_3 + h;
            intr += ((5 * h) / 24) * (11 * (Exponencial(a, b, xm_1, aux)) + (Exponencial(a, b, xm_2, aux)) + (Exponencial(a, b, xm_3, aux)) +
                                      11 * (Exponencial(a, b, xm_4, aux)));
        }
        num_part *= 2;
        erro_r = (double)fabs((intr - intr_aux) / intr);

    } while (erro_r > toler);

    return intr;

}

double NewtonCotes_4_aberta (double a, double b, double c1, double c2, double toler, int aux){
    double xi, intr = 0, intr_aux, xm_1, xm_2, xm_3, xm_4, xm_5, h, xf;
    int num_part = 1;
    int n = 0;
    double erro_r;
    double deltaX;
    do {
        intr_aux = intr;
        intr = 0;
        deltaX = delta_x(c1, c2, num_part);
        h = deltaX / 6;
        for (n = 0; n < num_part; n++) {
            xi = c1 + (n * deltaX);
            xm_1 = xi + h;
            xm_2 = xm_1 + h;
            xm_3 = xm_2 + h;
            xm_4 = xm_3 + h;
            xm_5 = xm_4 + h;
            xf = xm_5 + xf;
            intr += ((6 * h) / 20) *
                        (11 * (Exponencial(a, b, xm_1, aux)) - (14 * Exponencial(a, b, xm_2, aux)) + (26 * Exponencial(a, b, xm_3, aux)) -
                         (14 * Exponencial(a, b, xm_4, aux)) + (11 * Exponencial(a, b, xm_5, aux)));
        }
        num_part *= 2;
        erro_r = (double)fabs((intr - intr_aux) / intr);

    } while (erro_r > toler);

    return intr;

}

double GaussLegendre_2(double a, double b, double c1, double c2, double toler, int aux){
    double intr = 0, w[2], x[2], grau;
    int num_part = 1;
    int n = 0;
    w[0] = 1;
    w[1] = 1;
    x[0] = - 0.57735;
    x[1] = 0.57735;
    grau = 2;
    intr = intr_gauss(a, b, c1, c2, toler, w, x, grau, num_part, n, aux);
    return intr;
}


double GaussLegende_3(double a, double b, double c1, double c2, double toler, int aux){
    double intr = 0, w[3], x[3], grau;
    int num_part = 1;
    int n = 0;
    w[0] = 0.88888;
    w[1] = 0.55555;
    w[2] = 0.55555;
    x[0] = 0;
    x[1] = -0.77459;
    x[2] = 0.774596;
    grau = 3;
    intr = intr_gauss(a, b, c1, c2, toler, w, x, grau, num_part, n, aux);
    return intr;
}

double GaussLegende_4(double a, double b, double c1, double c2, double toler, int aux){
    double intr = 0, w[4], x[4], grau;
    int num_part = 1;
    int n = 0;
    w[0] = 0.65214;
    w[1] = 0.65214;
    w[2] = 0.34785;
    w[3] = 0.34785;
    x[0] = -0.3399;
    x[1] = 0.33998;
    x[2] = -0.8611;
    x[3] = 0.86113;
    grau = 4;
    intr = intr_gauss(a, b, c1, c2, toler, w, x, grau, num_part, n, aux);
    return intr;
}

double GaussLegendre_5(double a, double b, double c1, double c2, double toler, int aux){
    double intr = 0, w[5], x[5], grau;
    int num_part = 1;
    int n = 0;
    w[0] = 0.56888;
    w[1] = 0.47862;
    w[2] = 0.47862;
    w[3] = 0.23692;
    w[4] = 0.23692;
    x[0] = 0;
    x[1] = -0.53846;
    x[2] = 0.53846;
    x[3] = -0.90617;
    x[4] = 0.90617;

    //Utilizando a primeira estratégia (tolerância)
    grau = 5;
    intr = intr_gauss(a, b, c1, c2, toler, w, x, grau, num_part, n, aux);
    return intr;
}

/// Função que recebe o número do método digitado pelo usuário e retorna o método desejado
/// \param metodo     = metodo desejado
/// \param a          = o nosso limite inferior do intervalo
/// \param b          = o nosso limite superior do intervalo
/// \param toler      = tolerancia (caso utilizarmos a estratégia 1)
/// \param num_part_u = numero de partições fixas (caso utilizarmos a estratégia 2)
/// \param estrat     = estratégia que será utilizada
/// \return resultado da integral
double integral (int metodo, double a, double b, double c1, double c2, double toler, int aux){
    if (metodo == 1){
        return trapezoidal_fechada (a, b, c1, c2, toler, aux);
    }
    if (metodo == 2){
        return NewtonCotes_2_fechada (a, b, c1, c2, toler, aux);
    }
    if (metodo == 3){
        return NewtonCotes_3_fechada (a, b, c1, c2, toler, aux);
    }
    if (metodo == 4){
        return NewtonCotes_4_fechada (a, b, c1, c2, toler, aux);
    }
    if (metodo == 5){
        return trapezoidal_aberta (a, b, c1, c2, toler, aux);
    }
    if (metodo == 6){
        return NewtonCotes_2_aberta (a, b, c1, c2, toler, aux);
    }
    if (metodo == 7){
        return NewtonCotes_3_aberta (a, b, c1, c2, toler, aux);
    }
    if (metodo == 8){
        return NewtonCotes_4_aberta (a, b, c1, c2, toler, aux);
    }
    if (metodo == 9){
        return GaussLegendre_2(a, b, c1, c2, toler, aux);
    }
    if (metodo == 10){
        return GaussLegende_3 (a, b, c1, c2, toler, aux);
    }
    if (metodo == 11){
        return GaussLegende_4 (a, b, c1, c2, toler, aux);
    }
    if (metodo == 12){
        return GaussLegendre_5(a, b, c1, c2, toler, aux);
    }
}

/// Função que concentra a lógica do calculo da exponencial
/// \param a
/// \param b
/// \param toler_1
/// \param toler_2
/// \param metodo
/// \param aux
/// \return
double Questao03(double a, double b,double toler_1, double toler_2, int metodo, int aux){
    double c1, c2;
    double intr, intr_antr, erro_r;
    c1 = -1;
    c2 = 1;
    intr = 0;
    do {
        intr_antr = intr;
        intr = integral(metodo, a, b, c1, c2, toler_2, aux);
        c1 -= 1;
        c2 +=  1 ;
        erro_r = (double)fabs((intr - intr_antr) / intr);

    }while (erro_r > toler_1);

    return intr;
}

int main() {
    double toler_1, toler_2, a, b;
    int metodo;
    int aux;
    printf("\n--------Métodos:---------\n");
    printf("1 -  Newton Cotes Fechada Primeiro Grau\n2 -  Newton Cotes Fechada Segundo Grau\n3 -  Newton Cotes Fechada Terceiro Grau\n"
                   "4 -  Newton Cotes Fechada Quarto Grau\n5 -  Newton Cotes Aberta Primeiro Grau\n6 -  Newton Cotes Aberta Segundo Grau\n"
                   "7 -  Newton Cotes Aberta Terceiro Grau\n8 -  Newton Cotes Aberta Quarto Grau\n9 -  Gauss Legendre Grau 2\n10 - Gauss Legendre Grau 3\n"
                   "11 - Gauss Legendre Grau 4\n12 - Gauss Legendre Grau 5\n>>> ");
    //scanf("%d", &metodo);
    metodo = 1;
    printf("Entre com a tolerancia para os cortes:\n>>> ");
    //scanf("%lf", &toler_1);
    toler_1 = 0.001;
    printf("Entre com a tolerancia para a integral:\n>>> ");
    //scanf("%lf", &toler_2);
    toler_2 = 0.001;
    printf("Entre com a :\n>>> ");
    //scanf("%lf", &a);
    a = 0;
    printf("Entre com b :\n>>> ");
    //scanf("%lf", &b);
    b = 1;
    printf("Deseja utilizar qual estrategia?\n1 - Exponencial Simples\n2 - Exponencial Dupla :\n>>> ");
    //scanf("%d", &aux);
    aux = 1;
    printf("\n--------Resultado:---------\n");
    printf("%lf", Questao03(a, b, toler_1, toler_2, metodo, aux));
    system("pause");
}