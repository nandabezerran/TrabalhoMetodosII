//Fernanda Bezerra Nascimento - 388834
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef struct node node;

typedef struct element element;
/// Define a função que iremos utilizar para o método
/// \param x
/// \param y
/// \return o resultado da função
double Function(double x, double y){
    return (1 - x * y);
}

///Definição de Structs para o método de eliminação finita
struct node{
    double coordinate;
    int    id;
};

struct element{
    node    nodes[2];
    double **kMatrix;
    double  *fvector;
    double  *gvector;
};


/// Calcula equações diferenciais ordinarias pelo metodo de Forward Euler recebendo um valor inicial
/// \param initialPosition
/// \param initialTime
/// \param wantedTime
/// \param step
/// \return resultado da equação
double ForwardEuler(double initialPosition, double initialTime, double wantedTime, double step){
    double     nextPosition;
    double   actualPosition;
    double             time;

    nextPosition     = -0;
    actualPosition   = initialPosition;
    time             = initialTime;

    while (wantedTime > time){
        nextPosition   = actualPosition + (step * Function(actualPosition,time));
        time           = time + step;
        actualPosition = nextPosition;

    }

    return nextPosition;

}

/// Calcula a solução de uma equação diferencial ordinaria pelos metodos de runge-kutta
/// \param initialPosition
/// \param initialTime
/// \param wantedTime
/// \param step
/// \param order
/// \return solução da equação
double RangeKutta(double initialPosition, double initialTime, double wantedTime, double step, int order){
    double                   nextPosition;
    double        aproximatedNextPosition;
    double                 actualPosition;
    double                           time;
    double                       nextTime;

    nextPosition   = -0;
    actualPosition = initialPosition;
    time           = initialTime;     // Utilizado a partir do Runge-Kutta de ordem 2
    nextTime       = time + step;     // Utilizado a partir do Runge-Kutta de ordem 2

    if(order == 2){
        while (wantedTime > time){
            aproximatedNextPosition = actualPosition + step*(Function(actualPosition, time));
            nextPosition            = actualPosition + ((step/2) * (Function(actualPosition, time) +
                                                                 Function(aproximatedNextPosition, nextTime)));
            time      = time + step;
            nextTime  = time + step;

            actualPosition = nextPosition;

        }
        return nextPosition;
    }

    if(order == 3){
        double    aproximatedNextHalfPosition;
        double                   nextHalfTime;

        nextHalfTime   = time + (step/2);
        while (wantedTime > time){
            aproximatedNextHalfPosition = actualPosition + ((step/2) * (Function(actualPosition, time)));
            aproximatedNextPosition     = actualPosition + (step * Function(actualPosition, time));
            nextPosition                = actualPosition + ((step/6) * ((Function(actualPosition, time)) +
                                                            4*(Function(aproximatedNextHalfPosition, nextHalfTime))+
                                                            Function(aproximatedNextPosition, nextTime)));
            time          = time + step;
            nextHalfTime  = time + step/2;
            nextTime      = time + step;

            actualPosition = nextPosition;

        }
        return nextPosition;
    }

    if(order == 4){
        double      aproximatedNextThirdPosition;
        double   aproximatedNextTwoThirdPosition;
        double                     nextThirdTime;
        double                  nextTwoThirdTime;
        nextThirdTime     = time + (step/3);
        nextTwoThirdTime  = time + 2*(step/3);

        while (wantedTime > time){
            aproximatedNextThirdPosition    = actualPosition + ((step/3) * (Function(actualPosition, time)));
            aproximatedNextTwoThirdPosition = actualPosition  + ((2*(step/3)) * Function(actualPosition, time));
            aproximatedNextPosition         = actualPosition + (step * (Function(actualPosition, time)));
            nextPosition                    = actualPosition + ((step/8) * ((Function(actualPosition, time)) +
                                                     3*(Function(aproximatedNextThirdPosition, nextThirdTime))+
                                                    (3*Function(aproximatedNextTwoThirdPosition, nextTwoThirdTime))+
                                                     Function(aproximatedNextPosition, nextTime)));

            time              = time + step;
            nextThirdTime     = time + (step/3);
            nextTwoThirdTime  = time + 2*(step/3);
            nextTime          = time + step;

            actualPosition = nextPosition;

        }
        return nextPosition;
    }

}

double PredictorCorrectorMethod(double initialPosition, double initialTime, double wantedTime,
                                   double step, int order){
    double        nextPosition;
    double      actualPosition;
    double   predictedPosition;
    double   correctedPosition;
    double            nextTime;
    double                time;

    actualPosition    = initialPosition;
    time              = initialTime;
    nextTime          = time + step;
    correctedPosition = -0;

    if (order == 2) {
        while(wantedTime > time) {
            //Fase 0: Inicialização
            nextPosition = ForwardEuler(actualPosition, time, nextTime, step);

            //Fase 1: Predição
            predictedPosition = nextPosition + step * ((-1 / 2) * Function(actualPosition, time) +
                                                       (3 / 2) * Function(nextPosition, nextTime));

            //Fase 2: Correção
            correctedPosition = nextPosition + (step / 2) * (Function(actualPosition, nextTime) +
                                                             Function(predictedPosition, nextTime + step));

            time           = time + step;
            nextTime       = time + step;
            actualPosition = correctedPosition;
        }
        return correctedPosition;
    }

    if (order == 3) {
        double nextNextPosition;
        double nextNextTime = nextTime + step;
        while(wantedTime > time) {
            //Fase 0: Inicialização
            nextPosition = ForwardEuler(actualPosition, time, nextTime, step);
            nextNextPosition = ForwardEuler(nextPosition, nextTime, nextNextTime, step);

            //Fase 1: Predição
            predictedPosition = nextPosition + ((step/12) * (5* Function(actualPosition, time) -
                                                                         16 * Function(nextPosition, nextTime) +
                                                                         23 * Function(nextNextPosition,nextNextTime)));

            //Fase 2: Correção
            correctedPosition = nextPosition + ((step / 12) * (- Function(nextPosition, nextTime) +
                                                             8 * Function(nextNextPosition, nextNextTime) +
                                                             5 * Function(predictedPosition, nextNextTime + step)));

            time           = time + step;
            nextTime       = time + step;
            nextNextTime   = nextTime + step;
            actualPosition = correctedPosition;
        }
        return correctedPosition;
    }
}

int main() {
    double initialPosition = 1;
    double      intialTime = 0;
    double            step = 0.1;
    double      wantedTime = 0.1;

    int           rangeKuttaOrder = 2;
    int   predictorCorrectorOrder = 3;

    //printf("Forward Euler result: %f\n",ForwardEuler(initialPosition,intialTime,wantedTime,step));

    //printf("%d Order Runge-Kutta result: %f\n",rangeKuttaOrder,
          // RangeKutta(initialPosition,intialTime,wantedTime,step,rangeKuttaOrder));

    //printf("%d Order Predictor-Corrector result: %f\n",predictorCorrectorOrder,PredictorCorrectorMethod
                                                //(initialPosition,intialTime, wantedTime,step, predictorCorrectorOrder));


    system("pause");
    return 0;
}