//FERNANDA BEZERRA NASCIMENTO - SEGUNDA QUESTÃO - MATRÍCULA: 388834 - MÉTODOS NUMÉRICOS II
#include <stdio.h>
#include <math.h>
#include <tgmath.h>

double fun (double x){
    //Para uma função diferente para integral, edite a linha abaixo
    return cos(x);
}

double GaussHermit_2(){
    double intr = 0, w[2], x[2];

    w[0] = 0.88622;
    w[1] = 0.88622;
    x[0] = -0.70710;
    x[1] = 0.70710;
    for (int i = 1; i <= 2; i++){
        intr += w[i] * (fun(x[i]));
    }
    return intr;
}

double GaussHermit_3(){
    double intr = 0, w[3], x[3];

    w[0] = 0.29540;
    w[1] = 1.18163;
    w[2] = 0.29540;
    x[0] = -1.22474;
    x[1] = 0;
    x[2] = 1.22474;
    for (int i = 1; i <= 3; i++){
        intr += w[i] * (fun(x[i]));
    }
    return intr;
}

double GaussHermit_4(){
    double intr = 0, w[2], x[2];

    w[0] = 0.08131;
    w[1] = 0.80491;
    w[2] = 0.80491;
    w[3] = 0.08131;
    x[0] = -1.65068;
    x[1] = -0.52464;
    x[2] = 0.524647;
    x[3] = 1.650680;
    for (int i = 1; i <= 5; i++){
        intr += w[i] * (fun(x[i]));
    }
    return intr;
}

double GaussLaguerre_2(){
    double intr = 0, w[2], x[2];

    w[0] = 0.85355;
    w[1] = 0.14644;
    x[0] = 0.58578;
    x[1] = 3.41421;
    for (int i = 1; i <= 2; i++){
        intr += w[i] * (fun(x[i]));
    }
    return intr;
}

double GaussLaguerre_3(){
    double intr = 0, w[3], x[3];

    w[0] = 0.71109;
    w[1] = 0.27851;
    w[2] = 0.01038;
    x[0] = 0.41577;
    x[1] = 2.29428;
    x[2] = 6.28994;
    for (int i = 1; i <= 3; i++){
        intr += w[i] * (fun(x[i]));
    }
    return intr;
}

double GaussLaguerre_4(){
    double intr = 0, w[2], x[2];

    w[0] = 0.60315;
    w[1] = 0.35741;
    w[2] = 0.03888;
    w[3] = 5.39294;
    x[0] = 0.32254;
    x[1] = 1.74576;
    x[2] = 4.53662;
    x[3] = 9.39507;
    for (int i = 1; i <= 5; i++){
        intr += w[i] * (fun(x[i]));
    }
    return intr;
}

double fun_xk(int k, int n){
    return cos(((k - 1/2)* M_PI)/n);
}

double GaussChebyshev(int n){
    double intr = 0, w = M_PI/2;

    for (int i = 1; i <= n; i++){
        intr += w * (fun(fun_xk(i,n)));
    }
    return intr;
}

double integral (int metodo){
    if (metodo == 1){
        return GaussHermit_2 ();
    }
    if (metodo == 2){
        return GaussHermit_3 ();
    }
    if (metodo == 3){
        return GaussHermit_4 ();
    }
    if (metodo == 4){
        return GaussLaguerre_2();
    }
    if (metodo == 5){
        return GaussLaguerre_3();
    }
    if (metodo == 6){
        return GaussLaguerre_4();
    }

}
void main(int argc, char **argv){
    double a_b[6];
    int metodo;

    //Povoamento do vetor auxiliar(a_b) para o usuario escolher o intervalo
    a_b[0] = 0; a_b[1] = M_PI/8;
    a_b[2] = M_PI/4; a_b[3] = (3*M_PI)/3;
    a_b[4] = M_PI/2; a_b[5] = M_PI;



    //Escolha do método
    printf("\n--------Métodos:---------\n");
    printf("1 -  Gauss Hermit Grau 2\n2 -  Gauss Hermit Grau 3\n3 -  Gauss Hermit Grau 4\n"
                   "4 -  Gauss Laguerre Grau 2\n5 -  Gauss Laguerre Grau 3\n6 -  Gauss Laguerre Grau 4\n"
                   "7 -  Gauss Chebyshev \n>>> ");
    scanf("%d", &metodo);
    if (metodo == 7){
        int n;
        printf("Escolha o n:\n>>> ");
        scanf("%d", n);
        printf("\n--------Resultado:---------\n");
        printf("%lf",GaussChebyshev (n));

    //Exibição do resultado
    }
    else{
        printf("\n--------Resultado:---------\n");
        printf("%lf", integral(metodo));
    }
}