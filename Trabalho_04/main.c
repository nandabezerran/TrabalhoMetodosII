//Fernanda Bezerra Nascimento - 388834
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/// Define a função que iremos utilizar para o método
/// \param x
/// \param y
/// \return o resultado da função
double Function(double x, double y){
    return (x + y + x * y);
}

double ForwardEuler(double initialPosition, double initialTime, double wantedTime, double step){
    double     nextPosition;
    double   actualPosition;
    double             time;

    nextPosition     = -0;
    actualPosition   = initialPosition;
    time             = initialTime;

    while (wantedTime > time){
        nextPosition = actualPosition + (time * Function(actualPosition,time));
        time += step;
        actualPosition = nextPosition;

    }

    return nextPosition;

}

int main() {
    double initialPosition = 1;
    double      intialTime = 0;
    double            step = 0.025;
    double      wantedTime = 0.1;

    printf("Forward Euler result: %f\n",ForwardEuler(initialPosition,intialTime,wantedTime,step));

    system("pause");
    return 0;
}