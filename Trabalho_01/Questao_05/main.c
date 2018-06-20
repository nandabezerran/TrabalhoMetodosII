//FERNANDA BEZERRA NASCIMENTO - Integral Dupla Para Resolver O Problema da AP1 - MATRÍCULA: 388834 - MÉTODOS NUMÉRICOS II
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double fun (double x, double y){
    return cos(x*y);
}

/// Função para efetuar a mudança de variável, no caso calculariámos X(alfa)
/// \param x = o nosso alfa
/// \param a = o nosso limite inferior do intervalo
/// \param b = o nosso limite superior do intervalo
/// \return x(alfa)
double fun_x(double x, double a, double b){
    return (((a + b)/2) + (x * ((a - b)/2)));
}

double fun_y(double x_alfa, double y){
    return ((x_alfa/2) + (y * (x_alfa/2)));
}

double jacobiano(double x){
    return (M_PI/4) * x/2;
}

int main() {
    double intr = 0, w[3], x[3];
    double a = 0, b = M_PI/2;
    w[0] = 0.88888;
    w[1] = 0.55555;
    w[2] = 0.55555;
    x[0] = 0;
    x[1] = -0.77459;
    x[2] = 0.774596;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3 ; ++j) {
            intr+=w[i]*w[j]*fun(fun_x(x[i],a,b), fun_y(fun_x(x[i],a,b),x[j]))*jacobiano(fun_x(x[i],a,b));
        }

    }
    printf("%f", intr);
    system("pause");
    return 0;

}