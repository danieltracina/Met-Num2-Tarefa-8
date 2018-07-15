#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double x, double y) { // A função que aparece do lado direito da igualdade da EDP
	return(4.0);
}

double g(double x, double y) { // A função de contorno
	return(pow(x-y,2));
}

double u(double x, double y) { // Solução exata
	return(pow(x-y,2));
}

int main() {
	int n = 2, m = 4, N = 20000, l;
	double a = 0, b = 1, c = 0, d = 2, h, k, **matriz, *x, *y;
	double lambda, sigma, z, NORM, TOL = 0.0001;
	FILE *fp = fopen("solucao.txt", "w");
	
	//PASSOS DO ALGORITMO 12.1
	//Passo 1
	h = (b-a)/n;
	k = (d-c)/m;
	
	//Passo 2
	x = (double*)malloc((n+1)*sizeof(double));
	for(int i=0; i<n+1; i++)
		x[i] = a +i*h;
	
	//Passo 3
	y = (double*)malloc((m+1)*sizeof(double));
	for(int j=0; j<m+1; j++)
		y[j] = c +j*k;
	
	//Passo 4
	matriz = (double**)malloc((m+1)*sizeof(double*));
	for(int i=0; i<m+1; i++)
		matriz[i] = (double*)malloc((n+1)*sizeof(double));
	
	for(int i=0; i<m+1; i++)
		for(int j=0; j<n+1; j++)
			matriz[i][j] = 0;
	
	//Passo 5
	lambda = pow(h,2)/pow(k,2);
	sigma = 2*(1+lambda);
	l = 1;
	
	//Passo 6
	while(l<=N) {
		//Passo 7
		z = (-pow(h,2)*f(x[1],y[m-1])+g(a,y[m-1])+lambda*g(x[1],d)+lambda*matriz[1][m-2]+matriz[2][m-1])/sigma;
		
		NORM = fabs(z-matriz[1][m-1]);
		matriz[1][m-1] = z;
		
		//Passo 8
		for(int i=2; i<=n-2; i++) {
			z = (-pow(h,2)*f(x[i],y[m-1])+lambda*g(x[i],d)+matriz[i-1][m-1]+matriz[i+1][m-1]+lambda*matriz[i][m-2])/sigma;
			
			if(fabs(matriz[i][m-1]-z) > NORM)
				NORM = fabs(matriz[i][m-1]-z);
			matriz[i][m-1] = z;
		}
		
		//Passo 9
		z = (-pow(h,2)*f(x[n-1],y[m-1])+g(b,y[m-1])+lambda*g(x[n-1],d)+matriz[n-2][m-1]+lambda*matriz[n-1][m-2])/sigma;
		
		if(fabs(matriz[n-1][m-1]-z) > NORM)
			NORM = fabs(matriz[n-1][m-1]-z);
		matriz[n-1][m-1] = z;
		
		//Passo 10
		for(int j=m-2; j>=2; j--) {
			//Passo 11
			z = (-pow(h,2)*f(x[1],y[j])+g(a,y[j])+lambda*matriz[1][j+1]+lambda*matriz[1][j-1]+matriz[2][j])/sigma;
			if(fabs(matriz[1][j]-z) > NORM)
				NORM = fabs(matriz[1][j]-z);
			matriz[1][j] = z;
			
			//Passo 12
			for(int i=2; i<=n-2; i++) {
				z = (-pow(h,2)*f(x[i],y[j])+matriz[i-1][j]+lambda*matriz[i][j+1]+matriz[i+1][j]+lambda*matriz[i][j-1])/sigma; 
				if(fabs(matriz[i][j]-z) > NORM)
					NORM = fabs(matriz[i][j]-z);
				matriz[i][j] = z;
			}
			
			//Passo 13
			z = (-pow(h,2)*f(x[n-1],y[j])+g(b,y[j])+matriz[n-2][j]+lambda*matriz[n-1][j+1]+lambda*matriz[n-1][j-1])/sigma;
			if(fabs(matriz[n-1][j]-z) > NORM)
				NORM = fabs(matriz[n-1][j]-z);
			matriz[n-1][j] = z;
		}	
			//Passo 14
			z = (-pow(h,2)*f(x[1],y[1])+g(a,y[1])+lambda*g(x[1],c)+lambda*matriz[1][2]+matriz[2][1])/sigma;
			if(fabs(matriz[1][1]-z) > NORM)
				NORM = fabs(matriz[1][1]-z);
			matriz[1][1] = z;
			
			//Passo 15
			for(int i=2; i<=n-2; i++) {
				z = (-pow(h,2)*f(x[i],y[1])+lambda*g(x[i],c)+matriz[i-1][1]+lambda*matriz[i][2]+matriz[i+1][1])/sigma;
				if(fabs(matriz[i][1]-z) > NORM)
					NORM = fabs(matriz[i][1] -z);
				matriz[i][1] = z;
			}
			
			//Passo 16
			z = (-pow(h,2)*f(x[n-1],y[1])+g(b,y[1])+lambda*g(x[n-1],c)+matriz[n-2][1]+lambda*matriz[n-1][2])/sigma;
			if(fabs(matriz[n-1][1]-z) > NORM)
				NORM = fabs(matriz[n-1][1]-z);
			matriz[n-1][1] = z;
			
			//Passo 17
			if(NORM <= TOL) {
				//Passo 18 (Aqui ao invés de imprimir na tela, imprimirei os resultados no arquivo)
				fprintf(fp, "i\tj\tx[i]\t\ty[j]\t\tw[i][j]\t\tu(xi,yj)\n");
				for(int i=1; i<=n-1; i++)
					for(int j=1; j<=m-1; j++)
						fprintf(fp, "%d\t%d\t%lf\t%lf\t%lf\t%lf\n", i, j, x[i], y[j], matriz[i][j], u(x[i],y[j]));	
				//Passo 19
				printf("Procedimento bem sucedido;\n");
				break;
			}
			
			//Passo 20
			l++;
		}
		
	//Passo 21 (Aqui além de imprimir a mensagem de erro, irei imprimir a matriz encontrada depois da execução do algoritmo)
	if(l > N) {
		printf("\nProcedimento mal sucedido, número de iterações excedido.\n");
		fprintf(fp, "i\tj\tx[i]\t\ty[j]\t\tw[i][j]\t\tu(xi,yj)\n");
			for(int i=1; i<=n-1; i++)
				for(int j=1; j<=m-1; j++)
					fprintf(fp, "%d\t%d\t%lf\t%lf\t%lf\t%lf\n", i, j, x[i], y[j], matriz[i][j], u(x[i],y[j]));	
								
	}
	
	fclose(fp);
	return(0);
}	