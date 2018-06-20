//Fernanda Bezerra Nascimento - 388834
#include <stdio.h>
#include <stdlib.h>

/// Define a função que iremos utilizar para o método
/// \param x
/// \param y
/// \return o resultado da função
double Function(double x, double y){
    return (1 - x * y);
}

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

/// Calcula a solução de uma equação diferencial ordinaria pelos metodos do preditor corretor
/// \param initialPosition
/// \param initialTime
/// \param wantedTime
/// \param step
/// \param order
/// \return solução da equação
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
            predictedPosition = nextPosition + step * ((-1.0 / 2) * Function(actualPosition, time) +
                                                       (3.0 / 2) * Function(nextPosition, nextTime));

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
            predictedPosition = nextPosition + ((step/12) * (5.0* Function(actualPosition, time) -
                                                                         16.0 * Function(nextPosition, nextTime) +
                                                                         23.0 * Function(nextNextPosition,nextNextTime)));

            //Fase 2: Correção
            correctedPosition = nextPosition + ((step / 12) * (- Function(nextPosition, nextTime) +
                                                             8.0 * Function(nextNextPosition, nextNextTime) +
                                                             5.0 * Function(predictedPosition, nextNextTime + step)));

            time           = time + step;
            nextTime       = time + step;
            nextNextTime   = nextTime + step;
            actualPosition = correctedPosition;
        }
        return correctedPosition;
    }

    if (order == 4){
        double nextNextPosition;
        double nextNextNextPosition;
        double nextNextTime = nextTime + step;
        double nextNextNextTime = nextNextTime + step;
        while (wantedTime > time){
            //Fase 0: Inicialização
            nextPosition = ForwardEuler(actualPosition, time, nextTime, step);
            nextNextPosition = ForwardEuler(nextPosition, nextTime, nextNextTime, step);
            nextNextNextPosition = ForwardEuler(nextNextPosition, nextNextTime, nextNextNextTime, step);

            //Fase 1: Predição
            predictedPosition = nextPosition + (step/24) * ((-9.0 * Function(actualPosition, nextTime))
                                                              + (37.0 * Function(nextPosition, nextTime))
                                                              + (-59.0 * Function(nextNextPosition,nextNextTime))
                                                          + (55.0 * Function(nextNextNextPosition,nextNextNextPosition)));

            //Fase 2: Correção
            correctedPosition = nextPosition + (step/24) * (Function(nextPosition, nextTime) -5.0 * Function(nextNextPosition, nextNextTime)
                                                        +19.0 * Function(nextNextNextPosition, nextNextNextTime)
                                                        + 9.0 * Function(predictedPosition, nextNextNextTime + step));
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
    double      wantedTime = 2.4;

    int           rangeKuttaOrder;
    int   predictorCorrectorOrder;

    printf("Forward Euler result: %f\n",ForwardEuler(initialPosition,intialTime,wantedTime,step));

    for (rangeKuttaOrder = 2; rangeKuttaOrder <= 4; ++rangeKuttaOrder) {
        printf("%d Order Runge-Kutta result: %f\n",rangeKuttaOrder,
               RangeKutta(initialPosition,intialTime,wantedTime,step,rangeKuttaOrder));
    }


    for (predictorCorrectorOrder = 2; predictorCorrectorOrder < 5; ++predictorCorrectorOrder) {
        printf("%d Order Predictor-Corrector result: %f\n",predictorCorrectorOrder,PredictorCorrectorMethod
                (initialPosition,intialTime, wantedTime,step, predictorCorrectorOrder));
    }


    system("pause");
    return 0;
}