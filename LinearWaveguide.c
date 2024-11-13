//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//%%%%%%%%%%%%%%%%%%%% Circular Waveguide %%%%%%%%%%%%%%%%%%//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//%%%%%%%%%%%% Elaborado por: Igor Beder Burti R. %%%%%%%%%%//
//%%%%%%%%%%%%%% Última modificação: 06/06/2024 %%%%%%%%%%%%//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#define I _Complex_I

int main(){

    FILE *fidelidade20, *fidelidade3, *fidelidade4, *fidelidade5, *fidelidade6;
    FILE *fidelidadeN00N121920, *fidelidadeN00N45, *fidelidadeN00N56, *fidelidadeN00N35;
    FILE *squeeze2p, *squeeze2q;
    FILE *numeromedio1, *numeromedio20;
    int i, l, k, n;

    n = 20;

    double _Complex A[n+1][n+1];
    double J[n+1], J0, F20, F3, F4, F5, F6, FN56, FN34, FN45, FN35, s2, q2;
    double N1, N20;
    double _Complex k1A[n+1], k2A[n+1], k3A[n+1], k4A[n+1];
    double _Complex ic, x1;
    double t, tmax, dt, r, phi;

    //fidelidade20 = fopen("fidelidade20.dat", "w");
    //fidelidade3 = fopen("fidelidade3.dat", "w");
    //fidelidade4 = fopen("fidelidade4.dat", "w");
    //fidelidade5 = fopen("fidelidade5.dat", "w");
    //fidelidade6 = fopen("fidelidade6.dat", "w");
    //fidelidadeN00N121920 = fopen("fidelidadeN00N121920.dat", "w");
    //fidelidadeN00N45 = fopen("fidelidadeN00N45.dat", "w");
    //fidelidadeN00N56 = fopen("fidelidadeN00N56.dat", "w");
    //fidelidadeN00N35 = fopen("fidelidadeN00N35.dat", "w");
    //squeeze2p = fopen("squeeze3p.dat", "w");
    //squeeze2q = fopen("squeeze3q.dat", "w");
    numeromedio1 = fopen("NumeroMedio1_N20_N00N.dat", "w");
    //numeromedio20 = fopen("NumeroMedio20.dat", "w");

    ic = -1.0*I + 0.0;
    J0 = 1.0;

    for (i = 1; i <= n; i++){

        J[i] = J0*sqrt(i*(n-i));
        //J[i] = cos(M_PI + i)*cos(M_PI + i);
        //J[i] = 1.0;
    }

    for (i = 1; i <= n; i++){ // Condição inicial 

        for(k = 1; k <= n; k++){
            
            if(i == k){
            
                A[i][k] = 1.0;

            }

            else  A[i][k] = 0.0;

        }

    }

    dt = (2* M_PI)/999.0;
    dt = 0.001;
    tmax = 2* M_PI;
    //tmax = 100.0;
    for (t = dt; t <= tmax; t = t + dt){

        for (l = 1; l <= n; l++){
        
        /* Range-Kutta
        3 eqs.:
        1. i*dA_{1}/dt = J_1*A_{2} 
        2. i*dA_{N}/dt = J_{N-1}*A_{N-1})
        3. i*dA_{k}/dt = J_{k}*A_{k+1} + J_{k-1}*A_{k-1})
        */

        /* k1 */

        k1A[1] = ic*dt*J[1]*A[2][l]; // Lembrando que defini ic = -i
        k1A[n] = ic*dt*J[n-1]*A[n-1][l];

        for (k = 2; k <= n-1; k++){

        k1A[k] = ic*dt*(J[k-1]*A[k-1][l] + J[k]*A[k+1][l]);
        
        }

        /* k2 */         // k2 => f = f + h*k/2
        
        k2A[1] = ic*dt*J[1]*(A[2][l] + 0.5*k1A[2]);
        k2A[n] = ic*dt*J[n-1]*(A[n-1][l] + 0.5*k1A[n-1]);
        
        for(k = 2; k <= n-1; k++){

        k2A[k] = ic*dt*(J[k-1]*(A[k-1][l] + 0.5*k1A[k-1]) + J[k]*(A[k+1][l] + 0.5*k1A[k+1]));

        }

        /* k3 */         // k3 = k2
        
        k3A[1] = ic*dt*J[1]*(A[2][l] + 0.5*k2A[2]);
        k3A[n] = ic*dt*J[n-1]*(A[n-1][l] + 0.5*k2A[n-1]);        
        
        for(k = 2; k <= n-1; k++){

        k3A[k] = ic*dt*(J[k-1]*(A[k-1][l] + 0.5*k2A[k-1]) + J[k]*(A[k+1][l] + 0.5*k2A[k+1]));

        }

        /* k4 */         // k4 = k3*2
        
        k4A[1] = ic*dt*J[1]*(A[2][l] + k3A[2]);
        k4A[n] = ic*dt*J[n-1]*(A[n-1][l] + k3A[n-1]);
        
        for(k = 2; k <= n-1; k++){

        k4A[k] = ic*dt*(J[k-1]*(A[k-1][l] + k3A[k-1]) + J[k]*(A[k+1][l] + k3A[k+1]));

        }

        for (k = 1; k <= n; k++){ // Solução dos A_k

            A[k][l] = A[k][l] + (k1A[k]/6.0) + (k2A[k]/3.0) + (k3A[k]/3.0) + (k4A[k]/6.0);

        }

        }
        /* Calculo do número médio de fótons no waveguide j: N_j = n*|A_{ij}|^2*/

        // %%%%%%% Estado Fock %%%%%%% //

        //N1 = pow(cabs(A[1][1]), 2.0); //
        //N20 = 2.0*pow(cabs(A[1][20]), 2.0);
        //x1 = A[1][20]*A[1][20];
        //F20 = pow(cabs(x1), 2.0);
        //F5 = pow(cabs(A[1][5]), 2.0);
        //F6 = pow(cabs(A[1][6]), 2.0);

        // %%%%%%%%%%%%%%%%%%%%%%%%%%% //

        // %%%%%%% Estado N00N %%%%%%% //

        x1 = A[1][19]*A[1][19] + A[1][20]*A[1][20] + A[2][19]*A[2][19] + A[2][20]*A[2][20];
        FN34 = 0.25*cabs(x1)*cabs(x1);
        N1 = pow(cabs(A[1][1] + A[2][1]), 2.0);
        //x1 = A[4][1]*A[4][1] + A[4][2]*A[4][2] + A[5][1]*A[5][1] + A[5][2]*A[5][2];
        //FN45 =  0.25*cabs(x1)*cabs(x1);

        //x1 = A[3][1]*A[3][1] + A[3][2]*A[3][2] + A[5][1]*A[5][1] + A[5][2]*A[5][2];
        //FN35 =  0.25*cabs(x1)*cabs(x1);

        //x1 = A[5][1]*A[5][1] + A[5][2]*A[5][2] + A[6][1]*A[6][1] + A[6][2]*A[6][2];
        //FN56 =  0.25*cabs(x1)*cabs(x1);

        // %%%%%%%%%%%%%%%%%%%%%%%%%%% //

        // %%%%%%% Estado Squeeze %%%%%%% //
        //r = 0.7;
        //phi = M_PI;

        //x1 = A[2][1]*A[2][1]*exp(I*phi);
        //s2 = cabs(A[2][1])*cabs(A[2][1])*sinh(r)*sinh(r) + 0.25*sinh(2*r)*(x1 + conj(x1));
        //q2 = cabs(A[2][1])*cabs(A[2][1])*sinh(r)*sinh(r) - 0.25*sinh(2*r)*(x1 + conj(x1));


    //fprintf(fidelidade20, "%.8g %.8g\n", t, F20);
    //fprintf(fidelidade3, "%.8g %.8g\n", t, F3);
    //fprintf(fidelidade4, "%.8g %.8g\n", t, F4);
    //fprintf(fidelidade5, "%.8g %.8g\n", t, F5);
    //fprintf(fidelidade6, "%.8g %.8g\n", t, F6);
    //fprintf(fidelidadeN00N121920, "%.8g %.8g\n", t, FN34);
    //fprintf(fidelidadeN00N45, "%.8g %.8g\n", t, FN45);
    //fprintf(fidelidadeN00N35, "%.8g %.8g\n", t, FN35);
    //fprintf(fidelidadeN00N56, "%.8g %.8g\n", t, FN56);
    //fprintf(squeeze2p, "%.8g %.8g\n", t, s2);
    //fprintf(squeeze2q, "%.8g %.8g\n", t, q2);
    fprintf(numeromedio1, "%.8g %.8g\n", t, N1);
    //fprintf(numeromedio20, "%.8g %.8g\n", t, N20);

    }

    return 0;
}
