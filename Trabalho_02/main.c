//Fernanda Bezerra Nascimento - 388834
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////                             IMPLEMENTAÇÃO DO PROBLEMA DA AP1 DERIVAÇÃO                                  /////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Nesse problema temos 5 pontos f(x0), f(x0 + deltaX), f(x0 + 2*deltaX), f(x0 + 3*deltaX) e f(x0 + 4*deltaX)
    /// Temos como deltaX = 1/18
    /// Faremos a mudança de variável de x para alfa, tendo alfa base = 2
    /// Na reta x temos os pontos x-2, x-1, x0, x1, x2 e esses pontos são respecivamente representados na reta alfa como
    /// 0, 1, 2, 3, 4
    /// Fazemos a função g(alfa) = somatorio de 0 a 3 ( a k) delta^kf0
    /// Tiramos a derivada dessa função e então substituimos os valores

int main(){

    double deltaX = 1.0/18.0;
    double x0 = 1.0;
    double y  = 0.5;

    double  functionPoint0 = cos((x0 - 2*deltaX)* y),
            functionPoint1 = cos((x0 - deltaX)* y),
            functionPoint2 = cos((x0)*y),
            functionPoint3 = cos((x0 + deltaX)*y),
            functionPoint4 = cos((x0 + 2*deltaX)*y);
    double  gFunctionDerivative;

    int alfa = 2;
    gFunctionDerivative = (functionPoint1 - functionPoint0) + 0.5*(functionPoint2 - 2 * functionPoint1 + functionPoint0)
                           *(2 * alfa - 1) + (1/6)*(functionPoint3 - 3 * functionPoint2
                                                    + 3 * functionPoint1 - functionPoint0)
                                             *(3*pow(alfa,2) - 6 * alfa + 2) + (1/24)*(functionPoint4
                                              - 4 * functionPoint3+ 6 * functionPoint2 - 4 * functionPoint1
                                              + functionPoint0) * (4 * pow(alfa, 3)
                                                                                                                                               - 18 * pow(alfa, 2)
                                                                                                                                               + 22 * alfa - 6 );
    printf("Resultado da derivada de cos(xy) em relacao a x: %f\n", gFunctionDerivative);
    system("pause");
    return 0;

}