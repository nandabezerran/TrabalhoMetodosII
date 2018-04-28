//FERNANDA BEZERRA NASCIMENTO - PRIMEIRA QUESTÃO - MATRÍCULA: 388834 - MÉTODOS NUMÉRICOS II
#include <stdio.h>
#include <math.h>
#include <tgmath.h>

double fun (double x){
    //Para uma função diferente para integral, edite a linha abaixo
    return 1/ (1 + pow(sin(x),2));
}

double fun_x(double x, double a, double b){
    return (((a + b)/2) + (x * ((a - b)/2)));
}
double delta_x (double a, double b, int num_part){
    return fabs((b-a)/num_part);
}

double intr_gauss(double a, double b, double toler, int num_part_u, int estrat, double xi, double intr,
                  double intr_aux, double xf, const double *w, const double *x, double sum, double grau, int num_part,
                  int n, double erro_r, double deltaX) {
    if (estrat == 1) {
        do {
            intr_aux = intr;
            intr = 0;
            deltaX = delta_x(a, b, num_part);
            for (n = 0; n < num_part; n++) {
                xi = a + (n * deltaX);
                xf = xi + deltaX;
                sum = 0;
                for (int i = 0; i < grau; i++) {
                    sum += (w[i] * (fun(fun_x(x[i], xi, xf))));
                }
                intr += ((xf - xi) / 2) * sum;
            }

            num_part *= 2;
            erro_r = fabs((intr - intr_aux) / intr);

        } while (erro_r > toler);
    }
    else if (estrat == 2){
        intr = 0;
        deltaX = delta_x(a, b, num_part_u);
        for (n = 0; n < num_part_u; n++) {
            xi = a + (n * deltaX);
            xf = xi + deltaX;
            sum = 0;
            for (int i = 0; i < grau; i++) {
                sum += (w[i] * (fun(fun_x(x[i], xi, xf))));
            }
            intr += ((xf - xi) / 2) * sum;
        }
    }
    return intr;
}

double trapezoidal_fechada(double a, double b, double toler, int num_part_u, int estrat){
    double xi, intr = 0, xf, intr_aux, h;
    int n = 0;
    double erro_r;
    double deltaX;
    if (estrat == 1) {
        int num_part = 1;
        do {
            intr_aux = intr;
            intr = 0;
            deltaX = delta_x(a, b, num_part);
            h = deltaX / 2;
            for (n = 0; n < num_part; n++) {
                xi = a + (n * deltaX);
                xf = xi + deltaX;
                intr += (h * (fun(xi) + fun(xf)));
            }
            num_part *= 2;
            erro_r = fabs((intr - intr_aux) / intr);

        } while (erro_r > toler);
    }
    else if (estrat == 2){
        intr = 0;
        deltaX = delta_x(a, b, num_part_u);
        h = deltaX / 2;
        for (n = 0; n < num_part_u; n++) {
            xi = a + (n * deltaX);
            xf = xi + deltaX;
            intr += (h * (fun(xi) + fun(xf)));
        }
    }
    return intr;
}

double trapezoidal_aberta(double a, double b, double toler, int num_part_u, int estrat) {
    double xi, intr = 0, xf, intr_aux, xm_1, xm_2, h;
    int n = 0;
    double erro_r;
    double deltaX;
    if (estrat == 1) {
        int num_part = 1;
        do {
            intr_aux = intr;
            intr = 0;
            deltaX = delta_x(a, b, num_part);
            h = deltaX / 3;
            for (n = 0; n < num_part; n++) {
                xi = a + (n * deltaX);
                xm_1 = xi + h;
                xm_2 = xm_1 + h;
                intr += ((3 * h) / 2) * (fun(xm_1) + fun(xm_2));
            }
            num_part *= 2;
            erro_r = fabs((intr - intr_aux) / intr);

        } while (erro_r > toler);
    }
    else if (estrat == 2){
        intr = 0;
        deltaX = delta_x(a, b, num_part_u);
        h = deltaX / 3;
        for (n = 0; n < num_part_u; n++) {
            xi = a + (n * deltaX);
            xm_1 = xi + h;
            xm_2 = xm_1 + h;
            intr += ((3 * h) / 2) * (fun(xm_1) + fun(xm_2));
        }
    }
    return intr;
}

double NewtonCotes_2_fechada (double a, double b, double toler, int num_part_u, int estrat){
    double xi, intr = 0, xf, intr_aux, xm, h;
    int num_part = 1;
    int n = 0;
    double erro_r;
    double deltaX;
    if (estrat == 1) {
        do {
            intr_aux = intr;
            intr = 0;
            deltaX = delta_x(a, b, num_part);
            h = deltaX / 2;
            for (n = 0; n < num_part; n++) {
                xi = a + (n * deltaX);
                xf = xi + deltaX;
                xm = xi + h;
                intr += (h / 3) * ((fun(xi)) + 4 * (fun(xm)) + (fun(xf)));
            }
            num_part *= 2;
            erro_r = fabs((intr - intr_aux) / intr);

        } while (erro_r > toler);
    }
    else if (estrat == 2){
        intr = 0;
        deltaX = delta_x(a, b, num_part_u);
        h = deltaX / 2;
        for (n = 0; n < num_part_u; n++) {
            xi = a + (n * deltaX);
            xf = xi + deltaX;
            xm = xi + h;
            intr += (h / 3) * ((fun(xi)) + 4 * (fun(xm)) + (fun(xf)));
        }

    }
    return intr;

}

double NewtonCotes_3_fechada (double a, double b, double toler, int num_part_u, int estrat){
    double xi, intr = 0, xf, intr_aux, xm_1, xm_2, h;
    int num_part = 1;
    int n = 0;
    double erro_r;
    double deltaX;
    if (estrat == 1) {
        do {
            intr_aux = intr;
            intr = 0;
            deltaX = delta_x(a, b, num_part);
            h = deltaX / 3;
            for (n = 0; n < num_part; n++) {
                xi = a + (n * deltaX);
                xf = xi + deltaX;
                xm_1 = xi + h;
                xm_2 = xm_1 + h;
                intr += ((3 * h) / 8) *
                        ((fun(xi)) + 3 * (fun(xm_1)) + 3 * (fun(xm_2)) + (fun(xf)));
            }
            num_part *= 2;
            erro_r = fabs((intr - intr_aux) / intr);

        } while (erro_r > toler);
    }
    else if (estrat == 2){
        intr = 0;
        deltaX = delta_x(a, b, num_part_u);
        h = deltaX / 3;
        for (n = 0; n < num_part_u; n++) {
            xi = a + (n * deltaX);
            xf = xi + deltaX;
            xm_1 = xi + h;
            xm_2 = xm_1 + h;
            intr += ((3 * h) / 8) * ((fun(xi)) + 3 * (fun(xm_1)) + 3 * (fun(xm_2)) + (fun(xf)));
        }

    }
    return intr;

}

double NewtonCotes_4_fechada (double a, double b, double toler, int num_part_u, int estrat){
    double xi, intr = 0, xf, intr_aux, xm_1, xm_2, xm_3, h;
    int num_part = 1;
    int n = 0;
    double erro_r;
    double deltaX;
    if (estrat == 1) {
        do {
            intr_aux = intr;
            intr = 0;
            deltaX = delta_x(a, b, num_part);
            h = deltaX / 4;
            for (n = 0; n < num_part; n++) {
                xi = a + (n * deltaX);
                xf = xi + deltaX;
                xm_1 = xi + h;
                xm_2 = xm_1 + h;
                xm_3 = xm_2 + h;
                intr += ((2 * h) / 45) *
                        (7 * (fun(xi)) + 32 * (fun(xm_1)) + 12 * (fun(xm_2)) +
                         32 * (fun(xm_3)) + 7 * (fun(xf)));
            }
            num_part *= 2;
            erro_r = fabs((intr - intr_aux) / intr);

        } while (erro_r > toler);
    }
    else if(estrat == 2){
        intr = 0;
        deltaX = delta_x(a, b, num_part_u);
        h = deltaX / 4;
        for (n = 0; n < num_part_u; n++) {
            xi = a + (n * deltaX);
            xf = xi + deltaX;
            xm_1 = xi + h;
            xm_2 = xm_1 + h;
            xm_3 = xm_2 + h;
            intr += ((2 * h) / 45) *
                    (7 * (fun(xi)) + 32 * (fun(xm_1)) + 12 * (fun(xm_2)) +
                     32 * (fun(xm_3)) + 7 * (fun(xf)));
        }
    }
    return intr;

}

double NewtonCotes_2_aberta (double a, double b, double toler, int num_part_u, int estrat){
    double xi, intr = 0, intr_aux, xm_1, xm_2, xm_3, h;
    int num_part = 1;
    int n = 0;
    double erro_r;
    double deltaX;
    if (estrat == 1) {
        do {
            intr_aux = intr;
            intr = 0;
            deltaX = delta_x(a, b, num_part);
            h = deltaX / 4;
            for (n = 0; n < num_part; n++) {
                xi = a + (n * deltaX);
                xm_1 = xi + h;
                xm_2 = xm_1 + h;
                xm_3 = xm_2 + h;
                intr += ((4 * h) / 3) * (2 * (fun(xm_1)) - (fun(xm_2)) + 2 * (fun(xm_3)));
            }
            num_part *= 2;
            erro_r = fabs((intr - intr_aux) / intr);

        } while (erro_r > toler);
    }
    else if (estrat == 2){
        intr = 0;
        deltaX = delta_x(a, b, num_part_u);
        h = deltaX / 4;
        for (n = 0; n < num_part_u; n++) {
            xi = a + (n * deltaX);
            xm_1 = xi + h;
            xm_2 = xm_1 + h;
            xm_3 = xm_2 + h;
            intr += ((4 * h) / 3) * (2 * (fun(xm_1)) - (fun(xm_2)) + 2 * (fun(xm_3)));
        }
    }
    return intr;

}

double NewtonCotes_3_aberta (double a, double b, double toler, int num_part_u, int estrat){
    double xi, intr = 0, intr_aux, xm_1, xm_2, xm_3, xm_4, h;
    int num_part = 1;
    int n = 0;
    double erro_r;
    double deltaX;
    if (estrat == 1) {
        do {
            intr_aux = intr;
            intr = 0;
            deltaX = delta_x(a, b, num_part);
            h = deltaX / 5;
            for (n = 0; n < num_part; n++) {
                xi = a + (n * deltaX);
                xm_1 = xi + h;
                xm_2 = xm_1 + h;
                xm_3 = xm_2 + h;
                xm_4 = xm_3 + h;
                intr += ((5 * h) / 24) * (11 * (fun(xm_1)) + (fun(xm_2)) + (fun(xm_3)) +
                                          11 * (fun(xm_4)));
            }
            num_part *= 2;
            erro_r = fabs((intr - intr_aux) / intr);

        } while (erro_r > toler);
    }
    else if (estrat == 2){
        intr = 0;
        deltaX = delta_x(a, b, num_part_u);
        h = deltaX / 5;
        for (n = 0; n < num_part_u; n++) {
            xi = a + (n * deltaX);
            xm_1 = xi + h;
            xm_2 = xm_1 + h;
            xm_3 = xm_2 + h;
            xm_4 = xm_3 + h;
            intr += ((5 * h) / 24) * (11 * (fun(xm_1)) + (fun(xm_2)) + (fun(xm_3)) +
                                      11 * (fun(xm_4)));
        }
    }
    return intr;

}

double NewtonCotes_4_aberta (double a, double b, double toler, int num_part_u, int estrat){
    double xi, intr = 0, intr_aux, xm_1, xm_2, xm_3, xm_4, xm_5, h;
    int num_part = 1;
    int n = 0;
    double erro_r;
    double deltaX;
    if (estrat == 1) {
        do {
            intr_aux = intr;
            intr = 0;
            deltaX = delta_x(a, b, num_part);
            h = deltaX / 6;
            for (n = 0; n < num_part; n++) {
                xi = a + (n * deltaX);
                xm_1 = xi + h;
                xm_2 = xm_1 + h;
                xm_3 = xm_2 + h;
                xm_4 = xm_3 + h;
                xm_5 = xm_4 + h;
                intr += ((6 * h) / 20) *
                        (11 * (fun(xm_1)) - (14 * fun(xm_2)) + (26 * fun(xm_3)) -
                         (14 * fun(xm_4)) + (11 * fun(xm_5)));
            }
            num_part *= 2;
            erro_r = fabs((intr - intr_aux) / intr);

        } while (erro_r > toler);
    }
    else if (estrat == 2){
        intr = 0;
        deltaX = delta_x(a, b, num_part_u);
        h = deltaX / 6;
        for (n = 0; n < num_part_u; n++) {
            xi = a + (n * deltaX);
            xm_1 = xi + h;
            xm_2 = xm_1 + h;
            xm_3 = xm_2 + h;
            xm_4 = xm_3 + h;
            xm_5 = xm_4 + h;
            intr += ((6 * h) / 20) *
                    (11 * (fun(xm_1)) - (14 * fun(xm_2)) + (26 * fun(xm_3)) -
                     (14 * fun(xm_4)) + (11 * fun(xm_5)));
        }
    }
    return intr;

}

double GaussLegendre_2(double a, double b, double toler, int num_part_u, int estrat){
    double xi, intr = 0, intr_aux, xf, w[2], x[2], sum, grau;
    int num_part = 1;
    int n = 0;
    double erro_r;
    double deltaX;
    w[0] = 1;
    w[1] = 1;
    x[0] = - 0.57735;
    x[1] = 0.57735;
    grau = 2;
    intr = intr_gauss(a, b, toler, num_part_u, estrat, xi, intr, intr_aux, xf, w, x, sum, grau, num_part, n,
                      erro_r, deltaX);
    return intr;
}


double GaussLegende_3(double a, double b, double toler, int num_part_u, int estrat){
    double xi, intr = 0, intr_aux, xf, w[3], x[3], sum, grau;
    int num_part = 1;
    int n = 0;
    double erro_r;
    double deltaX;
    w[0] = 0.88888;
    w[1] = 0.55555;
    w[2] = 0.55555;
    x[0] = 0;
    x[1] = -0.77459;
    x[2] = 0.774596;
    grau = 3;
    intr = intr_gauss(a, b, toler, num_part_u, estrat, xi, intr, intr_aux, xf, w, x, sum, grau, num_part, n,
                      erro_r, deltaX);
    return intr;
}

double GaussLegende_4(double a, double b, double toler, int num_part_u, int estrat){
    double xi, intr = 0, intr_aux, xf, w[4], x[4], sum, grau;
    int num_part = 1;
    int n = 0;
    double erro_r;
    double deltaX;
    w[0] = 0.65214;
    w[1] = 0.65214;
    w[2] = 0.34785;
    w[3] = 0.34785;
    x[0] = -0.3399;
    x[1] = 0.33998;
    x[2] = -0.8611;
    x[3] = 0.86113;
    grau = 4;
    intr = intr_gauss(a, b, toler, num_part_u, estrat, xi, intr, intr_aux, xf, w, x, sum, grau, num_part, n,
                      erro_r, deltaX);
    return intr;
}

double GaussLegendre_5(double a, double b, double toler, int num_part_u, int estrat){
    double xi, intr = 0, intr_aux, xf, w[5], x[5], sum, grau;
    int num_part = 1;
    int n = 0;
    double erro_r;
    double deltaX;
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
    intr = intr_gauss(a, b, toler, num_part_u, estrat, xi, intr, intr_aux, xf, w, x, sum, grau, num_part, n,
                      erro_r, deltaX);
    return intr;
}

//Função que recebe o número do método digitado pelo usuário e retorna o método desejado
double integral (int metodo, double a, double b, double toler, int num_part_u, int estrat){
    if (metodo == 1){
        return trapezoidal_fechada (a, b, toler, num_part_u, estrat);
    }
    if (metodo == 2){
        return NewtonCotes_2_fechada (a, b, toler, num_part_u, estrat);
    }
    if (metodo == 3){
        return NewtonCotes_3_fechada (a, b, toler, num_part_u, estrat);
    }
    if (metodo == 4){
        return NewtonCotes_4_fechada (a, b, toler, num_part_u, estrat);
    }
    if (metodo == 5){
        return trapezoidal_aberta (a, b, toler, num_part_u, estrat);
    }
    if (metodo == 6){
        return NewtonCotes_2_aberta (a, b, toler, num_part_u, estrat);
    }
    if (metodo == 7){
        return NewtonCotes_3_aberta (a, b, toler, num_part_u, estrat);
    }
    if (metodo == 8){
        return NewtonCotes_4_aberta (a, b, toler, num_part_u, estrat);
    }
    if (metodo == 9){
        return GaussLegendre_2(a, b, toler, num_part_u, estrat);
    }
    if (metodo == 10){
        return GaussLegende_3 (a, b, toler, num_part_u, estrat);
    }
    if (metodo == 11){
        return GaussLegende_4 (a, b, toler, num_part_u, estrat);
    }
    if (metodo == 12){
        return GaussLegendre_5(a, b, toler, num_part_u, estrat);
    }
}

void main(int argc, char **argv){
    double a_b[6];
    int metodo;
    int a, b, estrat;
    int valid = 0;

    //Povoamento do vetor auxiliar(a_b) para o usuario escolher o intervalo
    a_b[0] = -1; a_b[1] = M_PI/8;
    a_b[2] = M_PI/4; a_b[3] = 1;//(3*M_PI)/2;
    a_b[4] = M_PI/2; a_b[5] = M_PI;

    //Menu de intervalos
    while (valid != 1){
        printf("\nEscolha o intervalo:\n");
        printf("----------Para a:-----------\n>>> 0 - 0	>>> 3 - 3π/2\n>>> 1 - π/8	>>> 4 - π/2\n>>> 2 - π/4	>>> 5 - π\n>>> ");
        //scanf("%d", &a);
        a = 0;
        printf("\n----------Para b:-----------\n>>> 0 - 0	>>> 3 - 3π/2\n>>> 1 - π/8	>>> 4 - π/2\n>>> 2 - π/4	>>> 5 - π\n>>> ");
        //scanf("%d", &b);
        b = (3*M_PI)/2;
        if (b >= a){
            valid = 1;
        }
        else {
            printf("Intervalo inválido, 'a' deve ser maior que 'b'\n");
        }
    }


    //Menu de estratégias
    printf("Escolha a estratégia:\n");
    printf("---------Estrategia:---------\n>>> 1 - Fornecer a tolerancia para a precisao desejada, epsilon \n>>> 2 - Fornecer o numero de partiçoes, N\n>>> ");
    scanf("%d", &estrat);


    //Escolha da tolerância ou número de partições
    double toler = 0;
    int num_part_u = 0;
    if (estrat == 1){
        printf("Entre com a tolerancia:\n>>> ");
        scanf("%lf", &toler);
    }

    else if (estrat == 2){
        printf("Entre com o numero de particoes:\n>>> ");
        scanf("%d", &num_part_u);
    }

    //Escolha do método
    printf("\n--------Métodos:---------\n");
    printf("1 -  Newton Cotes Fechada Primeiro Grau\n2 -  Newton Cotes Fechada Segundo Grau\n3 -  Newton Cotes Fechada Terceiro Grau\n"
                   "4 -  Newton Cotes Fechada Quarto Grau\n5 -  Newton Cotes Aberta Primeiro Grau\n6 -  Newton Cotes Aberta Segundo Grau\n"
                   "7 -  Newton Cotes Aberta Terceiro Grau\n8 -  Newton Cotes Aberta Quarto Grau\n9 -  Gauss Legendre Grau 2\n10 - Gauss Legendre Grau 3\n"
                   "11 - Gauss Legendre Grau 4\n12 - Gauss Legendre Grau 5\n>>> ");
    scanf("%d", &metodo);


    //Exibição do resultado
    printf("\n--------Resultado:---------\n");
    printf("%lf", integral(metodo ,a_b[a], a_b[b], toler, num_part_u, estrat));
}
