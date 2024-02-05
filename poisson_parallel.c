#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "omp.h"

#define MAX 1e+10 //Numero máximo de iterações
#define TOl 1e-13 //Tolerância do erro
#define h  1e-4 //Tamanho do passo x
#define k  1e-4 //Tamanho do passo y
#define l 4400
#define c 4400

void condcontorno(double **u, double **aux) //Inicializando as condições de contorno
{
	int i, y1, y2, r, s, v, o;
	double y, t;

	for(i=c-(c*0.125); i<=c+(c*0.125); i++) //condição dos "is" marca o tamanho dos passos na malha
	{
	    r = i+c*0.5;

	    t = i*(8*M_PI/l);
	    y = 0.125*l*sin(t); //equação da onda
	    y1= (int) y + l*0.5;
	    y2= (int) -1*y + l*0.5;

	    u[y1][r] = -1;
	    u[y2][r] = -1;
	    aux[y1][r] = -1;
		aux[y2][r] = -1;
	}

	for(i=0; i<(c*0.1); i++)
	{
	    s = i+c*0.1;
	    v = (9*c*0.1)-i;
	    o = c*0.1;

	    u[o][s] = 1;
	    u[s][o] = 1;
	    u[c-(o)][v] = 1;
	    u[v][c-(o)] = 1;

	    aux[o][s] = 1;
	    aux[s][o] = 1;
	    aux[c-o][v] = 1;
	    aux[v][c-o] = 1;
	}
}
double mtdffin(double **w, double **mt_aux)
{
	int i, j;
	double  df, err = 0, r;

#pragma omp parallel shared(mt_aux, df, w) private(i, j)
{

	for(i=1;i<l-1;i++)
	{
	    #pragma omp for
		for(j=1;j<c-1;j++)
		{
			df = w[i][j]; //Salvando o valor atual da matriz no ponto [i][j]
			if(mt_aux[i][j] == 0)
			{
				w[i][j] = (-h*h + w[i+1][j] + w[i-1][j] + w[i][j-1] + w[i][j+1])*0.25; //Formula encontrada para calcular os pontos pelo metodo das diferenças finitas.

				df = fabs(df - w[i][j]); //calculando o erro
			}
		}
	}
}
	if(df > err) //Verificando se o erro é menor que o anterior.
	{
		err = df;
	}

	return err;
}
void salvasol(double **w) //Salvando a solução no arquivo.
{
	FILE *fp;
	int i, j;

	fp = fopen("matriz_solucao.dat", "w");

	for(i=0;i<l;i++)
	{
		for(j=0;j<c;j++)
		{
			fprintf(fp, "%lf\t", w[i][j]);
		}

		fprintf(fp,"\n");
	}

 fclose(fp);
}
int main()
{
	double **M, **Maux, er = 1;
	int i, j, cont;

	M = (double**)malloc(l*sizeof(double*)); //Alocando a matriz principal.
	for(i=0;i<l;i++)
	{
		M[i] = (double*)malloc(c*sizeof(double*));
	}


	Maux = (double**)malloc(l*sizeof(double*)); //Alocando a matriz auxiliar.
	for(i=0;i<l;i++)
	{
		Maux[i] = (double*)malloc(c*sizeof(double*));
	}

	for(i=0;i<l;i++) //iniciandos as matrizes com zeros
	{
		for(j=0;j<c;j++)
		{
			M[i][j] = 0;
			Maux[i][j] = 0;
		}
	}

	condcontorno(M, Maux); //Aplicando as condições de contorno as matrizes.

	for(cont = 1; er>TOl && cont<MAX; cont++) //Aplicando o metodo, verificando se o erro está dentro da tolerancia e se o numero de iterações não foi exedido.
	{
		er = mtdffin(M, Maux);

		//printf("\nPasso: %d\n %.22lf\n", cont, er);
	}

	salvasol(M); //Salvando a solução

 return 0;
}
