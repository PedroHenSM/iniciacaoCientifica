/**
 * CRPSO + APM + variantes APM
 * Problemas clássicos de engenharia
 * Treliças
 * Primeiro objetivo: minimizar o peso
 * Restrições: tensão e deslocamento
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "test/foo.h"
#include "eureka/F101Truss10Bar.h"
#include "eureka/F103Truss25Bar.h"
#include "eureka/F105Truss60Bar.h"
#include "eureka/F107Truss72Bar.h"
// #include "eureka/F109Truss942Bar.h"
// #include "eureka/EurekaOptimaException.h"
// #include "eureka/TrussBarStructureStaticProblem.h"
// #include "eureka/TrussBarStructureStaticSimulator.h"
#include "eureka/Problem.h"


// #define USE_WELDED_MOD1_7
#define USE_STRING
#define N_VAR 3 // STRING: 3, REDUCER: 7, WELDED: 4, USE_WELDED_MOD1_7: 4, PRESSURE: 4, CANTILEVER: 10, T10C: 10, T10D: 10, T25D: 25, T52D: 52, T60C: 60, T60D: 60, T72C: 72, T942C: 942
#define N_CON 4 // STRING: 4, REDUCER: 11, WELDED: 5, USE_WELDED_MOD1_7: 7, PRESSURE: 4, CANTILEVER: 11, T10C: 18, T10D:18, T25D: 43, T52D: 32, T60C: 198, T72C: 240, T942C:
#define MAX_PAR 50
#define MAX_TIME 1500 // STRING: 400 (20k), REDUCER: 600 (30k) , WELDED_MOD: 400 (20K), PRESSURE: 600 (30k) (DeMelo)
// #define MAX_TIME 640 // S:720 R:720 W:640/6400 P:1600 C:700 T10C:5600 T10D:3600/1800 T25D:400 T52D:350 T60C:240 T72C:700 T942C:300(acho) 
#define MAX_RUN 30
#define INIT 999999999

#define D 3
#define Dkl 3

using namespace problem;

typedef struct best_particle_s {
	double * position;
	double fitness;
	double fitnessAPM;
	float * v;
} best_particle_t;

typedef struct particle_s {
	double * position;
	double * velocity;
	double fitness;
	double fitnessAPM;
	best_particle_t pBest;
	float * v;
} particle_t;

typedef struct swarm_s {
	particle_t particles[MAX_PAR];
	best_particle_t gBest;
} swarm_t;

//=======================================================================================================================================

void reduced_matrix (int numNodes, int n, double *A, double *elasticity, double *length, int **nodesCoord, double **C, double **kr, int numConst, int *in, int **vector, double **c, double *fr);
void LU (int numNodes, int numConst, double ** kr, double * fr, double * u, int * in);
void tension_stress (int n, int numNodes, double * u, double ** c, double * elasticity, double * A, double * length, int ** vector, double * F, double * sumConst, double * area);
void calculations (int n, int ** nodesCoord, int ** vector, int numNodes, double * nodesC, int numConst, double * nodesL, int * in, double * fr, double ** C, double * length);

/**
 * Truss function
 */	
void truss (double * A, double * sumConst, double * u, int cont) { 
	int i, j, n, m = D * 2;
	int * index, * in, ** nodesCoord, ** vector, numConst, numNodes, nodesConst, nodesLoad; // nodesConst = nós restritos; numConst = # de restricoes; n = num elementos (barras); nodesLoad = nós carregados
	char * input;
    double * area, * elasticity, * length, ** kr, * nodesC, * nodesL, * fr, ** C, ** coefXY, * F;
	FILE * inputTruss;

	#ifdef USE_T10C
	inputTruss = fopen("ins/input3D_t10.txt","r");
	if(inputTruss == NULL){	
		printf("\nFile open error...\n");
		exit(1);
	}
	#endif
	
	#ifdef USE_T10D
	inputTruss = fopen("ins/input3D_t10.txt","r");
	if(inputTruss == NULL){	
		printf("\nFile open error...\n");
		exit(1);
	}
	#endif

	#ifdef USE_T25C
	inputTruss = fopen("ins/input3D_t25.txt","r");
	if(inputTruss == NULL){	
		printf("\nFile open error...\n");
		exit(1);
	}
	#endif

	#ifdef USE_T25D
	inputTruss = fopen("ins/input3D_t25.txt","r");
	if(inputTruss == NULL){	
		printf("\nFile open error...\n");
		exit(1);
	}
	#endif

	#ifdef USE_T52C
	inputTruss = fopen("ins/input3D_t52.txt","r");
	if(inputTruss == NULL){	
		printf("\nFile open error...\n");
		exit(1);
	}
	#endif

	#ifdef USE_T52D
	inputTruss = fopen("ins/input3D_t52.txt","r");
	if(inputTruss == NULL){	
		printf("\nFile open error...\n");
		exit(1);
	}
	#endif

	#ifdef USE_T60C
	if (cont == 0)
		input = "ins/input3D_t60_0.txt";
	else if (cont == 1)
		input = "ins/input3D_t60_1.txt";
	else 
		input = "ins/input3D_t60_2.txt";
	inputTruss = fopen(input,"r");
	if(inputTruss == NULL){	
		printf("\nFile open error...\n");
		exit(1);
	}
	#endif


	#ifdef USE_T72C
	if (cont == 0)
		input = "ins/input3D_t72_0.txt";
	else 
		input = "ins/input3D_t72_1.txt";
	inputTruss = fopen(input,"r");
	if(inputTruss == NULL){	
		printf("\nFile open error...\n");
		exit(1);
	}
	#endif

	#ifdef USE_T942C
	inputTruss = fopen("ins/input3D_t942.txt","r");
	if(inputTruss == NULL){	
		printf("\nFile open error...\n");
		exit(1);
	}
	#endif

	fscanf(inputTruss,"%d %d", &n, &numNodes);

	//Allocates memory
	index = (int*) malloc(n * sizeof(int)); // armazena os indices dos elementos
	area = (double*) malloc(n * sizeof(double)); // areas das secoes de cada elemento	
	elasticity = (double*) malloc(n * sizeof(double)); // Elasticidade de cada elemento
	length = (double*) malloc(n * sizeof(double)); // Comprimento de cada elemento
	nodesCoord = (int**) malloc(n * sizeof(int*)); // Armazena os numNodess extremos dos elementos 
	for(i = 0; i < n; i++) {
		nodesCoord[i] = (int*) malloc(2 * sizeof(int));
	}
	C = (double**) malloc(numNodes * sizeof(double*));
	for (i = 0; i < numNodes; i++) {
		C[i] = (double*) malloc(D * sizeof(double));
	}
	nodesC = (double*) malloc(Dkl * numNodes * sizeof(double)); // Nos Restritos 
	nodesL = (double*) malloc(Dkl * numNodes * sizeof(double)); // Nos carregados
	in = (int*) malloc(Dkl * numNodes * sizeof(int));	
	coefXY = (double**) malloc(n * sizeof(double*)); // coeficientes cx, cy, ...
	for (i = 0;i < n; i++) {
		coefXY[i] = (double*) malloc(D * sizeof(double));
	}
	F = (double*) malloc(n * sizeof(double)); 
	vector = (int**) malloc(n * sizeof(int*)); // Monta o vetor de orientacao global
	for(i = 0; i < n; i++) {
    	vector[i] = (int*) malloc(m * sizeof(int));
	}

	//Reading
	numConst = 0;
	int i3;
	for (i = 0; i < numNodes; i++) {
		fscanf(inputTruss,"%d ", &i3);
		for(j = 0; j < Dkl; j++) {
			fscanf(inputTruss," %lf", &nodesC[i3*Dkl-(Dkl-j)]);
			if(nodesC[i3*Dkl-(Dkl-j)] == 1) 
				numConst++; 
		}
		for(j = 0; j < D; j++) {
			fscanf(inputTruss,"%lf ", &C[i3-1][j]);
			
		}
	}
		
	for (i = 0; i < D * numNodes; i++)
		nodesL[i] = 0;

	// Constraints and loads
	int aux, cont2, yy, xx;
		fscanf(inputTruss," %d", &nodesLoad);
		for (i = 0; i < nodesLoad; i++) {
			fscanf(inputTruss, "%d", &i3);
			fscanf(inputTruss, "%d", &aux);
			fscanf(inputTruss, "%lf", &nodesL[i3*Dkl-(Dkl-(aux-1))]);
		
		}

	for (i = 0; i < n; i++) {
		fscanf(inputTruss, "%d %lf %lf", &index[i], &elasticity[i], &area[i]);
	}

	for (i = 0; i < n; i++) {
		fscanf(inputTruss, "%d %d %d %d %d ", &index[i], &nodesCoord[i][0], &nodesCoord[i][1], &yy, &xx);
	}

	// Allocates memory
	kr = (double**) malloc((numNodes * Dkl - numConst + 1) * sizeof(double*));
	for(i = 0; i < (numNodes*Dkl-numConst+1); i++) {
		kr[i] =  (double*)malloc((numNodes*Dkl-numConst+1) * sizeof(double));    			    
	}

	fr =(double*) malloc((numNodes*Dkl-numConst+1) * sizeof(double));   
	
	fclose(inputTruss);

	calculations(n, nodesCoord, vector, numNodes, nodesC, numConst, nodesL, in, fr, C, length);
	reduced_matrix(numNodes, n, A, elasticity, length, nodesCoord, C, kr, numConst, in, vector, coefXY, fr);
	LU  (numNodes, numConst, kr, fr, u, in);
	tension_stress(n, numNodes, u, coefXY, elasticity, A, length, vector, F, sumConst, area);
	
	// Deallocates memory
	free (index);
	free (area);
	free (elasticity);
	free (length);
	free (nodesC);
	free (nodesL);
	free (in);
	free (F);
	free (fr);	
	for (i = 0; i < (numNodes*Dkl-numConst+1); i++)
		free (kr[i]);
	free (kr);
	for (i = 0; i < n; i++) {
		free (vector[i]);
		free (coefXY[i]);
		free (nodesCoord[i]);
	}
	free (vector);
	free (coefXY);
	free (nodesCoord);
	for (i = 0; i < numNodes; i++)
		free (C[i]);
	free (C);
}

/**
 * Reduced_matrix function
 */	
void reduced_matrix (int numNodes, int n, double * A, double * elasticity, double * length, int ** nodesCoord, double ** C, double ** kr, int numConst, int * in, int ** vector, double ** coefXY, double * fr) {
	int i,j;
	double k;

	// inicializacao da matriz kr
	for (i = 0; i < (numNodes*D-numConst); i++)
		for (j = 0; j < (numNodes*D-numConst); j++) {
		    kr[i][j] = 0.;
		}
	
	for(i = 0; i < n; i++) {
		for(j = 0; j < D; j++) {
			coefXY[i][j] = (C[nodesCoord[i][0]-1][j] - C[nodesCoord[i][1]-1][j]) / length[i];
		}

		k = A[i] * elasticity[i] / length[i];

		kr[in[vector[i][0]]][in[vector[i][0]]] += coefXY[i][0]*coefXY[i][0]*k;
		kr[in[vector[i][0]]][in[vector[i][1]]] += coefXY[i][1]*coefXY[i][0]*k;
		kr[in[vector[i][0]]][in[vector[i][2]]] += coefXY[i][2]*coefXY[i][0]*k;
		kr[in[vector[i][0]]][in[vector[i][3]]] +=-coefXY[i][0]*coefXY[i][0]*k;
		kr[in[vector[i][0]]][in[vector[i][4]]] +=-coefXY[i][1]*coefXY[i][0]*k;
		kr[in[vector[i][0]]][in[vector[i][5]]] +=-coefXY[i][2]*coefXY[i][0]*k;

		kr[in[vector[i][1]]][in[vector[i][0]]] += coefXY[i][0]*coefXY[i][1]*k;
		kr[in[vector[i][1]]][in[vector[i][1]]] += coefXY[i][1]*coefXY[i][1]*k;
		kr[in[vector[i][1]]][in[vector[i][2]]] += coefXY[i][2]*coefXY[i][1]*k;
		kr[in[vector[i][1]]][in[vector[i][3]]] +=-coefXY[i][0]*coefXY[i][1]*k;
		kr[in[vector[i][1]]][in[vector[i][4]]] +=-coefXY[i][1]*coefXY[i][1]*k;
		kr[in[vector[i][1]]][in[vector[i][5]]] +=-coefXY[i][2]*coefXY[i][1]*k;

		kr[in[vector[i][2]]][in[vector[i][0]]] += coefXY[i][0]*coefXY[i][2]*k;
		kr[in[vector[i][2]]][in[vector[i][1]]] += coefXY[i][1]*coefXY[i][2]*k;
		kr[in[vector[i][2]]][in[vector[i][2]]] += coefXY[i][2]*coefXY[i][2]*k;
		kr[in[vector[i][2]]][in[vector[i][3]]] +=-coefXY[i][0]*coefXY[i][2]*k;
		kr[in[vector[i][2]]][in[vector[i][4]]] +=-coefXY[i][1]*coefXY[i][2]*k;
		kr[in[vector[i][2]]][in[vector[i][5]]] +=-coefXY[i][2]*coefXY[i][2]*k;

		kr[in[vector[i][3]]][in[vector[i][0]]] +=-coefXY[i][0]*coefXY[i][0]*k;
		kr[in[vector[i][3]]][in[vector[i][1]]] +=-coefXY[i][1]*coefXY[i][0]*k;
		kr[in[vector[i][3]]][in[vector[i][2]]] +=-coefXY[i][2]*coefXY[i][0]*k;
		kr[in[vector[i][3]]][in[vector[i][3]]] += coefXY[i][0]*coefXY[i][0]*k;
		kr[in[vector[i][3]]][in[vector[i][4]]] += coefXY[i][1]*coefXY[i][0]*k;
		kr[in[vector[i][3]]][in[vector[i][5]]] += coefXY[i][2]*coefXY[i][0]*k;
	
		kr[in[vector[i][4]]][in[vector[i][0]]] +=-coefXY[i][0]*coefXY[i][1]*k;
		kr[in[vector[i][4]]][in[vector[i][1]]] +=-coefXY[i][1]*coefXY[i][1]*k;
		kr[in[vector[i][4]]][in[vector[i][2]]] +=-coefXY[i][2]*coefXY[i][1]*k;
		kr[in[vector[i][4]]][in[vector[i][3]]] += coefXY[i][0]*coefXY[i][1]*k;
		kr[in[vector[i][4]]][in[vector[i][4]]] += coefXY[i][1]*coefXY[i][1]*k;
		kr[in[vector[i][4]]][in[vector[i][5]]] += coefXY[i][2]*coefXY[i][1]*k;

		kr[in[vector[i][5]]][in[vector[i][0]]] +=-coefXY[i][0]*coefXY[i][2]*k;
		kr[in[vector[i][5]]][in[vector[i][1]]] +=-coefXY[i][1]*coefXY[i][2]*k;
		kr[in[vector[i][5]]][in[vector[i][2]]] +=-coefXY[i][2]*coefXY[i][2]*k;
		kr[in[vector[i][5]]][in[vector[i][3]]] += coefXY[i][0]*coefXY[i][2]*k;
		kr[in[vector[i][5]]][in[vector[i][4]]] += coefXY[i][1]*coefXY[i][2]*k;
		kr[in[vector[i][5]]][in[vector[i][5]]] += coefXY[i][2]*coefXY[i][2]*k;
	}    	
}

/**
 * LU function
 */	
void LU (int numNodes, int numConst, double ** kr, double * fr, double * u, int * in) {
	
	int n, i, j, k, l;
	n = (numNodes*D-numConst);
	double A[n][n+1], x[n], termo, m;
	
	// Definir A
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			A[i][j] = kr[i][j];
		}
	}	

	// acoplamento fr
	for (j = 0; j < n; j++) {
		A[j][n] = fr[j];
	}

	// Implementando Método de Gauss
	for (k = 0; k < n-1; k++) {
		for (i = k+1; i < n; i++) {
			// Multiplicadores
			m = -1 * (A[i][k] / A[k][k]);
			for (j = 0; j < n+1; j++) {
				A[i][j] = (A[k][j] * m) + A[i][j];
			}
		}
	}

	// Resolvendo o sistema
	for (i = 0; i < n; i++){
		termo = 0;	
		l = n - i;
		for (j=l; j<n;j++)
			termo = termo + (x[j] * A[n-i-1][j]);
		x[n-i-1] = (A[n-1-i][n] - termo) / A[n-i-1][n-i-1];
	}

	//impressao de U
	int c = 0;
	for (i = 0; i < (numNodes*D); i++) {
		u[i] = 0;
		if(in[i] < (numNodes*D-numConst)){
			u[i] = x[c];
			c++;
		}
	}
}

/**
 * Tension_stress function
 */	
void tension_stress(int n, int numNodes, double * u, double ** coefXY, double * elasticity, double * A, double * length, int ** vector, double * F, double * sumConst, double * area){

	int i,j;

	for (i = 0; i < n; i++) {
		sumConst[i] = 0;
		for(j = 0; j < D; j++)
			sumConst[i] += coefXY[i][j] * u[vector[i][j]] - coefXY[i][j] * u[vector[i][j+D]];
		sumConst[i] *= elasticity[i] / length[i];
		F[i] = sumConst[i] * A[i];
	}

			
	/*printf ("\nNode\tX-DISPLACEMENT\tY-DISPLACEMENT\tZ-DISPLACEMENT\n");
	for (i=0; i<(numNodes); i++) {
		printf("\n %d\t", i);
		for(j=0; j<D; j++)
			printf ("%lf\t",u[D*i+j]);
	}

	printf ("\n\nElement\t\tForce\t\tStress\nNumber\n");
	for (i=0; i<n; i++) 
		printf (" %d\t\t%e\t%e\n",i,F[i],sumConst[i]);	
	


 	 printf ("\nElemento_2\tcx\t\tcy\t\tcz\n");
	 for (i=0; i<n; i++) {
	 printf ("\n%d\t\t",i);
	 for (j = 0; j < D; j++) 
	 printf ("%lf\t",coefXY[i][j]);	
			    }

	printf("\nVetor area\n");
	for (i=0; i<n; i++) 
		printf("%f\n",area[i]);
	printf("\nVetor elasticity\n");
	for (i=0; i<n; i++)
		printf("%f\n",elasticity[i]);
	
	//printf ("\n\n");*/
}

/**
 * Calculations function
 */	
void calculations (int n, int ** nodesCoord, int ** vector, int numNodes, double * nodesC, int numConst, double * nodesL, int * in, double * fr, double ** C, double * length) {

	int i,j;	
	// Calculo dos comprimentos
	double c1; 
	for (i = 0; i < n; i++) {
		c1 = 0;
		for(j = 0; j < D; j++) {
			c1 += pow((C[nodesCoord[i][0]-1][j] - C[nodesCoord[i][1]-1][j]), 2);
		}
		length[i] = pow(c1,0.5);
		//printf("length[%d] = %lf;\n", i, length[i]);
	}
	//exit(0);
	
	for (i = 0; i < n; i++) {
		for(j = 0; j < D; j++) {
			vector[i][j] = nodesCoord[i][0] * D - (D- j);
			vector[i][j+D] = nodesCoord[i][1] * D - (D - j);
		}
	}

	int cn = 0;   
	for (i = 0; i < D*numNodes; i++) {
		if(nodesC[i] == 1) {
			cn++;
			in[i] = -17;
		} else {
			in[i] = i - cn;
		}
	}

	for (i = 0; i < D*numNodes; i++)
		if(in[i] == -17) {
			in[i] = D * numNodes - cn;
		}

	// Vector fr
	for (i = 0; i < D*numNodes; i++) 
		if(in[i] < (numNodes*D-numConst)) {
			fr[in[i]] = nodesL[i];
		}
}

//=============================================================================================================================================

/**
 * Allocate function - allocates the swarm that will be used.
 */	
swarm_t * allocate () {
	int i;
	swarm_t * swarm = (swarm_t *) malloc (sizeof(swarm_t));	
	swarm->gBest.position = (double *) malloc (N_VAR * sizeof(double));
	swarm->gBest.v = (float* )malloc (N_CON * sizeof(double));

	for (i = 0; i < MAX_PAR; i++) {
		swarm->particles[i].position = (double *) malloc (N_VAR * sizeof(double));
		swarm->particles[i].velocity = (double *) malloc (N_VAR * sizeof(double));
		swarm->particles[i].v = (float *) malloc (N_CON * sizeof(double));
		swarm->particles[i].pBest.position = (double *) malloc (N_VAR * sizeof(double));
		swarm->particles[i].pBest.v = (float *) malloc (N_CON * sizeof(double));
	}

	return swarm;
}

/**
 * Boundary function - boundary condition.
 */
void boundary (double lowerBound[], double upperBound[]) {
	int m;

	#ifdef USE_STRING
		lowerBound[0] = 2;
		upperBound[0] = 15; 
		lowerBound[1] = 0.25;
		upperBound[1] = 1.3;
		lowerBound[2] = 0.05;
		upperBound[2] = 2; 
	#endif

	#ifdef USE_REDUCER
		// lowerBound[0] = 2.6;   
		// upperBound[0] = 3.6;
		// lowerBound[1] = 0.7;   
		// upperBound[1] = 0.8;
		// lowerBound[2] = 17;   
		// upperBound[2] = 28;
		// lowerBound[3] = 7.3;   
		// upperBound[3] = 8.3;
		// lowerBound[4] = 7.8;   
		// upperBound[4] = 8.3;
		// lowerBound[5] = 2.9;   
		// upperBound[5] = 3.9;	
		// lowerBound[6] = 0;   
		// upperBound[6] = 20;		

		lowerBound[0] = 2.6;   
		upperBound[0] = 3.6;
		lowerBound[1] = 0.7;   
		upperBound[1] = 0.8;
		lowerBound[2] = 17;   
		upperBound[2] = 28;
		lowerBound[3] = 7.3;   
		upperBound[3] = 8.3;
		lowerBound[4] = 7.3;   
		upperBound[4] = 8.3;
		lowerBound[5] = 2.9;   
		upperBound[5] = 3.9;	
		lowerBound[6] = 5;   
		upperBound[6] = 5.5;		
	#endif

	#ifdef USE_WELDED
		lowerBound[0] = 0.125;
		upperBound[0] = 10;

		for (m = 1; m < N_VAR; m++) {
			lowerBound[m] = 0.1;  
			upperBound[m] = 10; 
		} 
	#endif

	#ifdef USE_WELDED_MOD1_7
		lowerBound[0] = 0.1;
		upperBound[0] = 2;

		lowerBound[1] = 0.1;
		upperBound[1] = 10;

		lowerBound[2] = 0.1;
		upperBound[2] = 10;

		lowerBound[3] = 0.1;
		upperBound[3] = 2;
	#endif

	#ifdef USE_PRESSURE
		for (m = 0; m < 2; m++) {
			lowerBound[m] = 1;  
			upperBound[m] = 80; 
		}
		for (m = 2; m < N_VAR; m++) {
			lowerBound[m] = 10;  
			upperBound[m] = 200; 
		}
	#endif

	#ifdef USE_CANTILEVER
		lowerBound[0] = 1;      
		upperBound[0] = 64;
		lowerBound[1] = 1;      
		upperBound[1] = 64;
		lowerBound[2] = 0;      
		upperBound[2] = 3;
		lowerBound[3] = 0;      
		upperBound[3] = 3;
		lowerBound[4] = 0;      
		upperBound[4] = 3;
		lowerBound[5] = 0;      
		upperBound[5] = 3;
		lowerBound[6] = 1;      
		upperBound[6] = 64; 
		lowerBound[7] = 1;      
		upperBound[7] = 64; 
		lowerBound[8] = 1;      
		upperBound[8] = 64;
		lowerBound[9] = 1;      
		upperBound[9] = 64; 
	#endif

	#ifdef USE_T10C
		for (m = 0; m < N_VAR; m++) {
			lowerBound[m] = 0.1;
			upperBound[m] = 33.50;
		} 
	#endif

	#ifdef USE_T10D
		for (m = 0; m < N_VAR; m++) {
			lowerBound[m] = 0;
			upperBound[m] = 31;
		} 
	#endif

	#ifdef USE_T25C
		for (m = 0; m < N_VAR; m++) {
			lowerBound[m] = 0.1;
			upperBound[m] = 3.4;
		} 
	#endif

	#ifdef USE_T25D
		for (m = 0; m < N_VAR; m++) {
			lowerBound[m] = 0;
			upperBound[m] = 29;
		} 
	#endif

	#ifdef USE_T52C
		for (m = 0; m < N_VAR; m++) {
			lowerBound[m] = 11.102791;
			upperBound[m] = 3350.831008;
		} 
	#endif

	#ifdef USE_T52D
		for (m = 0; m < N_VAR; m++) {
			lowerBound[m] = 0;
			upperBound[m] = 63;
		} 
	#endif

	#ifdef USE_T60C
		for (m = 0; m < N_VAR; m++) {
			lowerBound[m] = 0.5;
			upperBound[m] = 5;
		} 
	#endif

	#ifdef USE_T72C
		for (m = 0; m < N_VAR; m++) {
			lowerBound[m] = 0.1;
			upperBound[m] = 5;
		} 
	#endif

	#ifdef USE_T942C
		for (m = 0; m < N_VAR; m++) {
			lowerBound[m] = 1;
			upperBound[m] = 200;
		} 
	#endif
}

/**
 * Initialize function - initialize the position and velocity randomly.
 */
void initialize (swarm_t * swarm, double lowerBound[], double upperBound[]) {
	int i, j, m;
	float rand1, rand2;
	
	swarm->gBest.fitness = INIT;
	swarm->gBest.fitnessAPM = INIT;
	for (i = 0; i < N_CON; i++)	{
		swarm->gBest.v[i] = INIT;
	}

	for (i = 0; i < MAX_PAR; i++) {
		swarm->particles[i].fitness = INIT;
		swarm->particles[i].pBest.fitness = INIT;
		swarm->particles[i].fitnessAPM = INIT;
		swarm->particles[i].pBest.fitnessAPM = INIT;
		for (m = 0; m < N_CON; m++)	{
			swarm->particles[i].pBest.v[m] = INIT;
		}
		for (j = 0; j < N_VAR; j++) {
			// TODO PEDRO apagar depois 
			rand1 = (rand() % 10000);
			rand1 /= 10000;
			swarm->particles[i].position[j] = lowerBound[j] + rand1 * (upperBound[j] - lowerBound[j]);  

			// float bestT60C[25] = {1.14806, 2.0408, 0.501062, 1.85304, 1.75507, 0.585083, 1.9266, 1.8511, 1.00357, 1.88154, 1.86348, 0.5, 2.02915, 1.24715, 1.09826, 0.637386, 0.730659, 0.972064, 1.12724, 1.14845, 1.05482, 1.05667, 0.654224, 1.02509, 1.25693};
			// swarm->particles[i].position[j] = bestT60C[j];
			//printf(" %d = Position: %.12lf  \tVelocity: ", j, swarm->particles[i].position[j]);
		
			rand2 =  (rand() % 10000);
			rand2 /= 10000;
			swarm->particles[i].velocity[j] = lowerBound[j] + rand2 * (upperBound[j] - lowerBound[j]);  
			//printf("%.12lf\n", swarm->particles[i].velocity[j]);
		}
		//printf("\n");
	}
}


/**
 * Constraint function
 */
void constraint (swarm_t * swarm, double ** violationAcum, double * violation, double * sumConst, int * numViolation, double valuesArray[MAX_PAR][N_CON+1]) {
	int i, j, m, cont1, cont2, cont3;
	double x[N_VAR], g[N_CON];
	// printf("\ndentro de constraint");
	// exit(66);

	for (i = 0; i < MAX_PAR; i++) {
		for (j = 0; j < N_VAR; j++) {
			x[j] = swarm->particles[i].position[j];
			//printf("x(con) = %lf ", x[j]);
		}

		#ifdef USE_STRING
			/*x[0] = 11.852177;
			x[1] = 0.34747463;
			x[2] = 0.051301897;*/
			g[0] = 1 - x[1] * x[1] * x[1] * x[0] / (71785 * pow(x[2], 4));
			g[1] = (4 * x[1] * x[1] - x[2] * x[1]) / (12566 * (x[1] * pow(x[2], 3) - pow(x[2], 4))) + 1 / (5108 * x[2] * x[2]) - 1;
			g[2] = 1 - 140.45 * x[2] / (x[1] * x[1] * x[0]);
			g[3] = (x[1] + x[2]) / 1.5 - 1;
		#endif

		#ifdef USE_REDUCER
			/*x[0] = 3.5;
			x[1] = 0.7;
			x[2] = 17;
			x[3] = 7.3000035;
			x[4] = 7.7115326;
			x[5] = 3.350216;
			x[6] = 5.286655;*/
			g[0] = 27 * pow(x[0], -1) * pow(x[1], -2) * pow(x[2], -1) - 1;
			g[1] = 397.5 * pow(x[0], -1) * pow(x[1], -2) * pow(x[2], -2) - 1;
			g[2] = 1.93 * pow(x[1], -1) * pow(x[2], -1) * pow(x[3], 3) * pow(x[5], -4) - 1;
			g[3] = 1.93 * pow(x[1], -1) * pow(x[2], -1) * pow(x[4], 3) * pow(x[6], -4) - 1;
			g[4] = ((pow(pow((745 * x[3]) / (x[1] * x[2]), 2) + 16900000, 0.5)) / (0.1 * pow(x[5], 3))) - 1100;
			g[5] = ((pow(pow((745 * x[4]) / (x[1] * x[2]), 2) + 157500000, 0.5)) / (0.1 * pow(x[6], 3))) - 850;
			g[6] = x[1] * x[2] - 40;
			g[7] = -(x[0] / x[1]) + 5;
			g[8] = (x[0] / x[1]) - 12;
			g[9] = ((1.5 * x[5] + 1.9) * pow(x[3], -1)) - 1;
			g[10] = ((1.1 * x[6] + 1.9) * pow(x[4], -1)) - 1;
		#endif

		#ifdef USE_WELDED
			/*x[0] = 0.2442949; 
			x[1] = 6.2116738;
			x[2] = 8.3015486;
			x[3] = 0.2443003;*/
			double tau, sigma, Pc, delta, tau1, tau2, alfa, numerador, denominador;
			alfa = sqrt (0.25 * (x[1] * x[1] + pow(x[0] + x[2], 2)));
			tau1 = 6000 / (sqrt(2) * x[0] * x[1]);
			numerador = 6000 * (14 + 0.5 * x[1]) * alfa;
			denominador = 2 * (0.707 * x[0] * x[1] * (x[1]*x[1]/12 + 0.25 * (pow(x[0] + x[2], 2))));
			tau2 = numerador / denominador;
			tau = sqrt((tau1 * tau1) + (tau2 * tau2) + ((x[1] * tau1 * tau2) / alfa));
			sigma = 504000 / (x[2] * x[2] * x[3]);
			Pc = 64746.022 * (1 - 0.0282346 * x[2]) * x[2] * x[3] * x[3] * x[3];
			delta = 2.1952 / (x[2] * x[2] * x[2] * x[3]);
			g[0] = tau - 13600;
			g[1] = sigma - 30000;
			g[2] = x[0] - x[3];
			g[3] = - Pc + 6000;
			g[4] = delta - 0.25 ;
			//g[0] = 0;
		#endif

		#ifdef USE_WELDED_MOD1_7
		  double p = 6000;
			double l = 14;
			double e = 30 * (pow(10, 6));
			double G = 12 * (pow(10, 6));
			double tauMax = 13600;
			double sigmaMax = 30000;
			double deltaMax = 0.25;

			double m = p * (l + (x[1] / 2));
			double r = sqrt( (pow(x[1], 2) / 4) + pow(((x[0] + x[2]) / 2), 2) );
			double j = 2 * ( sqrt(2) * x[0] * x[1] * ((pow(x[1], 2) / 12) + pow( (x[0] + x[2] / 2), 2)));

			double tauLinha = p / (sqrt(2) * x[0] * x[1]);
			double tauDuasLinhas = (m * r) / j;
			double sigmaX = (6 * p * l) / ( x[3] * pow(x[2], 2));
			double deltaX = (4 * p * pow(l, 3)) / (e * pow(x[2], 3) * x[3]);
			double pcX = ((4.013 * e * sqrt((pow(x[2], 2) * pow(x[3], 6)) / 36 )) / pow(l,2)) * (1 - (x[2] / (2 * l)) * sqrt(e / (4*G)));
			double tauX = sqrt( pow(tauLinha, 2) + 2 * tauLinha * tauDuasLinhas * (x[1] / (2*r)) + pow(tauDuasLinhas, 2));


			g[0] = tauX - tauMax; // # <= 0
			g[1] = sigmaX - sigmaMax; // # <= 0
			g[2] = x[0] - x[3]; // # <= 0
			g[3] = 0.10471 * pow(x[0], 2) + 0.04811 * x[2] * x[3] * (14.0 + x[1]) - 5; // # <= 0
			g[4] = 0.125 - x[0]; // # <= 0
			g[5] = deltaX - deltaMax; // # <= 0
			g[6] = p - pcX; // # <= 0
		#endif

		#ifdef USE_PRESSURE
			int aux;
			aux = round(x[0]);
			x[0] = aux * 0.0625;
			aux = round(x[1]);
			x[1] = aux * 0.0625;
			/*x[0] = 0.8125;
			x[1] = 0.4375;
			x[2] = 42.086994;
			x[3] = 176.779128;*/
			double pi = 3.141592654;
			g[0] = -x[0] + 0.0193 * x[2];
			g[1] = -x[1] + 0.00954 * x[2];
			g[2] = -(pi * x[2] * x[2] * x[3]) - ((4. / 3.) * pi * x[2] * x[2] * x[2]) + 1296000;
			g[3] = x[3] - 240;
	
		#endif

		#ifdef USE_CANTILEVER
			double B[] = {2.4, 2.6, 2.8, 3.1}; 
			double H[] = {45.0, 50.0, 55.0, 60.0};
			int aux, cont;
			x[0] = (int) x[0];
			x[1] = (int) x[1];
			aux = round(x[2]);
			x[2] = B[aux];
			aux = round(x[4]);
			x[4] = B[aux];
			aux = round(x[3]);
			x[3] = H[aux];
			aux = round(x[5]);
			x[5] = H[aux];

			/*x[0] = 3; 
			x[1] = 60;
			x[2] = 3.1;
			x[3] = 60;
			x[4] = 2.6;
			x[5] = 50;
			x[6] = 2.2094;
			x[7] = 44.0428;
			x[8] = 2.0944;
			x[9] = 31.9867;*/

			g[0] = 6.*25*1000000/(x[0]*x[1]*x[1]) - 14000;
			g[1] = 6.*20*1000000/(x[2]*x[3]*x[3]) - 14000;
			g[2] = 6.*15*1000000/(x[4]*x[5]*x[5]) - 14000;
			g[3] = 6.*10*1000000/(x[6]*x[7]*x[7]) - 14000;
			g[4] = 6.*5*1000000/(x[8]*x[9]*x[9]) - 14000;
			for (cont = 0; cont < 5; cont++ ) {
				g[cont + 5] = x[cont * 2 + 1] / x[cont * 2] - 20;
			}
			double aj1, aj2, aj3, aj4, aj5;	
			aj1 = x[0] * pow(x[1], 3) / 12.;
			aj2 = x[2] * pow(x[3], 3) / 12.;
			aj3 = x[4] * pow(x[5], 3) / 12.;
			aj4 = x[6] * pow(x[7], 3) / 12.;
			aj5 = x[8] * pow(x[9], 3) / 12.;
			g[10] = (0.0025 / 3.) * ((1000000 / aj5) + (7 * 1000000 / aj4) + (19 * 1000000 / aj3) + (37 * 1000000 / aj2) + (61 * 1000000 / aj1)) - 2.7;
		#endif

		#ifdef USE_T10C
			/*x[0] = 30.162525; 
			x[1] = 0.10003946;
			x[2] = 22.81192;
			x[3] = 15.871827;
			x[4] = 0.10000233;
			x[5] = 0.5149511;
			x[6] = 7.505953;
			x[7] = 21.264076;
			x[8] = 21.383036;
			x[9] = 0.10000795;*/
			// double aux, u[8];

			// truss (x, sumConst, u, 0);
			
			// for (m = 0; m < 8; m++) {
			// 	//printf("\nu = %lf ", u[m]);
			// 	if(u[m] <= 0)
			// 		aux = -u[m];
			// 	else
			// 		aux = u[m];
			// 	g[m] = aux - 2;
			// }

			// for (m = 0; m < N_VAR; m++) {
			// 	//printf(" sumConst = %lf\n", sumConst[m]);
			// 	if(sumConst[m] <= 0)
			// 		aux = -sumConst[m];
			// 	else
			// 		aux = sumConst[m]; 
			// 	g[m + 8] = aux - 25000;
			// }

			// PEDRO
			// Passes the constraints values to 'g' array
			// (the first position from 'valuesArray' is the objective Function)
			for (m = 0; m < N_CON; m++){
				g[m] = valuesArray[i][m+1];
			}
		#endif

		#ifdef USE_T10D
			/*x[0] = 33.5; 
			x[1] = 1.62;
			x[2] = 22.9;
			x[3] = 14.2;
			x[4] = 1.62;
			x[5] = 1.62;
			x[6] = 7.97;
			x[7] = 22.9;
			x[8] = 22.0;
			x[9] = 1.62;*/

			// double aux, u[8];

			// double tab[] = {1.62, 1.80, 1.99, 2.13, 2.38, 2.62, 2.93, 3.13, 3.38, 3.47, 3.55, 3.63, 3.88, 4.22, 4.49, 4.59, 4.80, 4.97, 5.12, 5.74, 7.97, 11.50, 13.50, 14.20, 15.50, 16.90, 18.80, 19.90, 22.00, 26.50, 30.00, 33.50};
			// int pos;
			// for (j = 0; j < N_VAR; j++) {
			// 	pos = round(x[j]);
			// 	x[j] = tab[pos];
			// }

			// truss (x, sumConst, u, 0);

			// for (m = 0; m < 8; m++) {
			// 	//printf("\nu = %lf ", u[m]);
			// 	if(u[m] <= 0)
			// 		aux = -u[m];
			// 	else
			// 		aux = u[m];
			// 	g[m] = aux - 2;
			// }

			// for (m = 0; m < N_VAR; m++) {
			// 	//printf(" sumConst = %lf\n", sumConst[m]);
			// 	if(sumConst[m] <= 0)
			// 		aux = -sumConst[m];
			// 	else
			// 		aux = sumConst[m]; 
			// 	g[m + 8] = aux - 25000;
			// }

			// PEDRO
			// Passes the constraints values to 'g' array
			// (the first position from 'valuesArray' is the objective Function)
			for (m = 0; m < N_CON; m++){
				g[m] = valuesArray[i][m+1];
			}
		#endif

		#ifdef USE_T25C

			double aux, u[18];

			truss (x, sumConst, u, 0);
			
			for (m = 0; m < 18; m++) {
				//printf("\nu = %lf ", u[m]);
				if(u[m] <= 0)
					aux = -u[m];
				else
					aux = u[m];
				g[m] = aux - 0.35;
			}

			for (m = 0; m < N_VAR; m++) {
				//printf(" sumConst = %lf\n", sumConst[m]);
				if(sumConst[m] <= 0)
					aux = -sumConst[m];
				else
					aux = sumConst[m]; 
				g[m + 18] = aux - 40000;
			}

			// PEDRO
			// Passes the constraints values to 'g' array
			// (the first position from 'valuesArray' is the objective Function)
			// for (m = 0; m < N_CON; m++){
			for (m = 0; m < 29; m++){
				g[m] = valuesArray[i][m+1];
			}

		#endif

		#ifdef USE_T25D

			// double aux, u[18];

			// double tab[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6,0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.8, 3.0, 3.2, 3.4};
			// int pos;
			// for (j = 0; j < N_VAR; j++) {
			// 	pos = round(x[j]);
			// 	x[j] = tab[pos];
			// }

			// truss (x, sumConst, u, 0);
			
			// for (m = 0; m < 18; m++) {
			// 	//printf("\nu = %lf ", u[m]);
			// 	if(u[m] <= 0)
			// 		aux = -u[m];
			// 	else
			// 		aux = u[m];
			// 	g[m] = aux - 0.35;
			// }

			// for (m = 0; m < N_VAR; m++) {
			// 	//printf(" sumConst = %lf\n", sumConst[m]);
			// 	if(sumConst[m] <= 0)
			// 		aux = -sumConst[m];
			// 	else
			// 		aux = sumConst[m]; 
			// 	g[m + 18] = aux - 40000;
			// }

			// PEDRO
			// Passes the constraints values to 'g' array
			// (the first position from 'valuesArray' is the objective Function)
			for (m = 0; m < N_CON; m++){
				g[m] = valuesArray[i][m+1];
			}

		#endif

		#ifdef USE_T52C

			double aux, u[32];

			truss (x, sumConst, u, 0);
		
			for (m = 0; m < N_VAR; m++) {
				if(sumConst[m] <= 0)
					aux = -sumConst[m];
				else
					aux = sumConst[m]; 
				g[m] = aux - 180;
			}

		#endif

		#ifdef USE_T52D

			double aux, u[32];

			double tab[] = {11.102791, 14.103566, 19.604806, 25.006202, 30.707597, 39.109767, 44.211008, 56.313953, 60.214884, 76.619070, 78.519535, 99.424651, 100.024806, 113.062946, 126.631473, 145.736434, 156.338760, 162.040155, 180.044651, 199.049302, 213.052868, 238.059070, 262.064961, 263.065271, 288.071473, 293.072713, 309.076589, 313.077674, 338.083876, 347.086047, 355.088062, 362.975349, 384.095194, 387.095969, 388.096279, 418.103721, 422.104651, 444.925271, 459.113798, 480.119070, 497.123256, 512.126977, 574.142326, 722.179070, 797.197674, 853.211628, 930.232248, 1085.269147, 1150.285271, 1350.334884, 1390.344806, 1420.352248, 1550.384496, 1600.396899, 1690.419225, 1880.466357, 1990.493643, 2200.545736, 2290.568062, 2450.607752, 2650.657364, 2800.694574, 3000.744186, 3350.831008};
			int pos;
			for (j = 0; j < N_VAR; j++) {
				pos = round(x[j]);
				x[j] = tab[pos];
			}

			truss (x, sumConst, u, 0);
			
			for (m = 0; m < N_VAR; m++) {
				if(sumConst[m] <= 0)
					aux = -sumConst[m];
				else
					aux = sumConst[m]; 
				g[m] = aux - 180;
			}

		#endif
		
		#ifdef USE_T60C
			// int cont, cont2, casos = 3;
			// double aux, u[45], auxG[N_CON / 3];
			
			// for (cont = 0; cont < casos; cont++) {
			// 	truss (x, sumConst, u, cont);
					
			// 	for (m = 0; m < 45; m++) {
			// 		if(u[m] <= 0)
			// 			aux = -u[m];
			// 		else
			// 			aux = u[m];
			// 		if (m == 6)
			// 			auxG[0] = aux - 1.75;
			// 		else if (m == 7)
			// 			auxG[1] = aux - 1.75;
			// 		else if (m == 22)
			// 			auxG[2] = aux - 2.25;
			// 		else if (m == 23)
			// 			auxG[3] = aux - 2.25;
			// 		else if (m == 33)
			// 			auxG[4] = aux - 2.75;
			// 		else if (m == 34)
			// 			auxG[5] = aux - 2.75;

			// 	}

			// 	for (m = 0; m < N_VAR; m++) {
			// 		if(sumConst[m] <= 0)
			// 			aux = -sumConst[m];
			// 		else
			// 			aux = sumConst[m]; 
			// 		auxG[m + 6] = aux - 10000;
			// 	}

			// 	for (cont2 = 0; cont2 < (N_CON / 3); cont2++) {
			// 		if (cont == 0)
			// 			g[cont2] = auxG[cont2]; 
			// 		else if (cont == 1)
			// 			g[cont2 + 66] = auxG[cont2]; 	
			// 		else 
			// 			g[cont2 + 132] = auxG[cont2]; 	
			// 	}
			// }

			// PEDRO
			// Passes the constraints values to 'g' array
			// (the first position from 'valuesArray' is the objective Function)
			for (m = 0; m < N_CON; m++){
				g[m] = valuesArray[i][m+1];
			}

		#endif


		#ifdef USE_T72C
			// int cont, cont2, casos = 2;
			// double aux, u[48], auxG[N_CON / 2];
			
			// for (cont = 0; cont < casos; cont++) {
			// 	truss (x, sumConst, u, cont);
					
			// 	for (m = 0; m < 48; m++) {
			// 		if(u[m] <= 0)
			// 			aux = -u[m];
			// 		else
			// 			aux = u[m];
			// 		auxG[m] = aux - 0.25;
			// 	}

			// 	for (m = 0; m < N_VAR; m++) {
			// 		if(sumConst[m] <= 0)
			// 			aux = -sumConst[m];
			// 		else
			// 			aux = sumConst[m]; 
			// 		auxG[m + 48] = aux - 25000;
			// 	}
			// 	for (cont2 = 0; cont2 < (N_CON / 2); cont2++) {
			// 		if (cont == 0)
			// 			g[cont2] = auxG[cont2]; 
			// 		else
			// 			g[cont2 + 120] = auxG[cont2]; 	
			// 	}	
			// }

			// PEDRO
			// Passes the constraints values to 'g' array
			// (the first position from 'valuesArray' is the objective Function)
			for (m = 0; m < N_CON; m++){
				g[m] = valuesArray[i][m+1];
			}
		#endif

		#ifdef USE_T942C
			double aux, u[696];

			truss (x, sumConst, u, 0);
			
			for (m = 0; m < 696; m++) {
				//printf("\nu = %lf ", u[m]);
				if(u[m] <= 0)
					aux = -u[m];
				else
					aux = u[m];
				g[m] = aux - 15;
			}

			for (m = 0; m < N_VAR; m++) {
				//printf(" sumConst = %lf\n", sumConst[m]);
				if(sumConst[m] <= 0)
					aux = -sumConst[m];
				else
					aux = sumConst[m]; 
				g[m + 696] = aux - 25000;
			}
		#endif

		// for (int it = 0; it < N_CON; it++){
		// 	printf("g[%d]: %f", it, g[it]);
		// }
		// printf("dentro de constraint");
		// exit(1);

		sumConst[i] = 0.;

		for (int m = 0; m < N_CON; m++) {
			// printf("g[%d]: %f", m, g[m]);
			if (g[m] > 0) { 
				violation[m] += g[m];
				numViolation[m]++;
				swarm->particles[i].v[m] = g[m]; 
				violationAcum[i][m] += g[m];
				sumConst[i] += g[m];
			} else {
				swarm->particles[i].v[m] = 0.;
			} 	
		} 
	}
	// printf("ultima linha constraint");
	// exit(1);
}

/**
 * APM function - penalty method
 */
void APM (swarm_t * swarm, double functionAvg, double * k, double ** violationAcum, double valuesArray[MAX_PAR][N_CON+1]) {
	int i, j, numViolation[N_CON];
	double violation[N_CON], sumConst[MAX_PAR], functionAvgK, function, aux, quadrado = 0.;
	

	if (functionAvg > 0)
		functionAvgK = functionAvg;
	else
		functionAvgK = -functionAvg;

	for (i = 0; i < N_CON; i++) {
		violation[i] = 0.; 
		numViolation[i] = 0;
	}
	for (i = 0; i < MAX_PAR; i++) {
		sumConst[i] = 1;
		for(j = 0; j < N_CON; j++) {
			violationAcum[i][j] = 1;
		}
	}

	// printf("pre constraint");
	constraint (swarm, violationAcum, violation, sumConst, numViolation, valuesArray);
	// printf("pós constraint");

	for (i = 0; i < N_CON; i++) { 
		violation[i] = violation[i] / MAX_PAR; 
	}
	
	// Calculate k
	for (i = 0; i < N_CON; i++) {
		quadrado += pow(violation[i], 2);
	}

	for (i = 0; i < N_CON; i++) {
		if (quadrado != 0.)	{
			k[i] = (functionAvgK * violation[i]) / quadrado;
			//printf("\nk = %lf quadrado = %lf violation = %lf AvgK = %lf\n", k[i], quadrado, violation[i], functionAvgK);
		}
		else
			k[i] = 0.;
	}

	/*for (i = 0; i < N_CON; i++)
		printf("k = %.12lf ", k[i]);*/

	// particles
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].fitness; 

		aux = 0.;
		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].v[j] * k[j];
		}

		if(aux != 0.) {
			if(function < functionAvg) 
				function = functionAvg;	
			function += aux;
		}

		swarm->particles[i].fitnessAPM = function;
	}

	// pBest
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].pBest.fitness; 
		aux = 0.;
		 
		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].pBest.v[j] * k[j];
		}

		if(aux != 0.) {
			if(function < functionAvg) 
				function = functionAvg;	
			function += aux;
		}

		swarm->particles[i].pBest.fitnessAPM = function;
		  
	}

	// gBest
	function = swarm->gBest.fitness; 
	aux = 0.;
		 
	for (j = 0; j < N_CON; j++) {
		aux += swarm->gBest.v[j] * k[j];
	}

	if(aux != 0) {
		if(function < functionAvg) 
			function = functionAvg;	
		function += aux;
	}

	swarm->gBest.fitnessAPM = function;	 
}

/**
 * Sporadic APM function - penalty method
 */
void sporadic_APM (swarm_t * swarm, double functionAvg, double * k,  double ** violationAcum, int generation, double valuesArray[MAX_PAR][N_CON+1]) {
	int i, j, numViolation[N_CON], frequency = MAX_TIME / 10;
	double violation[N_CON], sumConst[MAX_PAR], functionAvgK, function, aux;
	
	if(functionAvg > 0) 
		functionAvgK = functionAvg;
	else 
		functionAvgK = -functionAvg;
	
	for (i = 0; i < N_CON; i++) {
		violation[i] = 0.;
		numViolation[i] = 0; 
	}

	constraint (swarm, violationAcum, violation, sumConst, numViolation, valuesArray);

	// Calculate k
	if(generation % frequency == 0) {

		for (i = 0; i < N_CON; i++) 
			violation[i] = violation[i] / MAX_PAR;

		double quadrado = 0.;	
		for (i = 0; i < N_CON; i++) {
			quadrado += pow(violation[i], 2);
		}
		for (i = 0; i < N_CON; i++) {
			if (quadrado != 0.)	
				k[i] = functionAvgK * violation[i] / quadrado;
			else
				k[i] = 0.;
		}
	}

	// particles
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].fitness; 
		aux = 0.;

		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].v[j] * k[j];
		}

		if(aux != 0) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].fitnessAPM = function;
	}
	// pBest
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].pBest.fitness; 
		aux = 0.;
		 
		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].pBest.v[j] * k[j];
		}

		if(aux != 0.) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].pBest.fitnessAPM = function;
		  
	}

	// gBest
	function = swarm->gBest.fitness; 
	aux = 0.;
		 
	for (j = 0; j < N_CON; j++) {
		aux += swarm->gBest.v[j] * k[j];
	}

	if(aux != 0) {
		if(function < functionAvg) 
			function = functionAvg;	
	}
	function += aux;
	swarm->gBest.fitnessAPM = function;

}

/**
 * Sporadic Acumulation APM function - penalty method
 */
void sporadic_acumulation_APM (swarm_t * swarm, double functionAvg, double * k,  double ** violationAcum, int generation, double valuesArray[MAX_PAR][N_CON+1]) {
	int i, j, numViolation[N_CON], frequency = MAX_TIME / 10;
	double violation[N_CON], sumConst[MAX_PAR], functionAvgK, function, aux;
	
	if(functionAvg > 0) 
		functionAvgK = functionAvg;
	else 
		functionAvgK = -functionAvg;

	for (i = 0; i < N_CON; i++) {
		numViolation[i] = 0;
	}

	constraint (swarm, violationAcum, violation, sumConst, numViolation, valuesArray);	

	// Calculate k
	if(generation % frequency == 0) {
		for (i = 0; i < N_CON; i++) {
			violation[i] = 0.; 
		}	
		for (i = 0; i < MAX_PAR; i++) {
			for (j = 0; j < N_CON; j++) {
				violation[j] += violationAcum[i][j];
			}
		} 
		for (i = 0; i < N_CON; i++) 
			violation[i] = violation[i] / MAX_PAR;

		double quadrado = 0.;	
		for (i = 0; i < N_CON; i++) {
			quadrado += pow(violation[i], 2);
		}
		for (i = 0; i < N_CON; i++) {
			if (quadrado != 0.)	
				k[i] = functionAvgK * violation[i] / quadrado;
			else
				k[i] = 0.;
		}
	}

	// particles
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].fitness; 
		aux = 0.;

		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].v[j] * k[j];
		}

		if(aux != 0) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].fitnessAPM = function;
	}
	// pBest
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].pBest.fitness; 
		aux = 0.;
		 
		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].pBest.v[j] * k[j];
		}

		if(aux != 0.) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].pBest.fitnessAPM = function;
		  
	}

	// gBest
	function = swarm->gBest.fitness; 
	aux = 0.;
		 
	for (j = 0; j < N_CON; j++) {
		aux += swarm->gBest.v[j] * k[j];
	}

	if(aux != 0) {
		if(function < functionAvg) 
			function = functionAvg;	
	}
	function += aux;
	swarm->gBest.fitnessAPM = function;

	if(generation % frequency == 0) {
		for (i = 0; i < MAX_PAR; i++)
			for (j = 0; j < N_CON; j++) 
				 violationAcum[i][j] = 0.; 
	}

}

/**
 * Monotonic APM function - penalty method
 */
void monotonic_APM (swarm_t * swarm, double functionAvg, double * k_aux,  double ** violationAcum, double valuesArray[MAX_PAR][N_CON+1]) {
	int i, j, numViolation[N_CON];
	double violation[N_CON], sumConst[MAX_PAR], k[N_CON], functionAvgK, function, aux;
	
	if(functionAvg > 0) 
		functionAvgK = functionAvg;
	else 
		functionAvgK = -functionAvg;

	for (i = 0; i < N_CON; i++) {
		violation[i] = 0.; 
		numViolation[i] = 0;
	}
	
	constraint (swarm, violationAcum, violation, sumConst, numViolation, valuesArray);	

	for (i = 0; i < N_CON; i++) 
			violation[i] = violation[i] / MAX_PAR;

	// Calculate k
	double quadrado = 0.;	
	for (i = 0; i < N_CON; i++) {
		quadrado += pow(violation[i], 2);
	}
	
	for (i = 0; i < N_CON; i++) {
			if (quadrado != 0.)	
				k[i] = functionAvgK * violation[i] / quadrado;
			else
				k[i] = 0.;
		
		if(k[i] > k_aux[i])
			k_aux[i] = k[i];
	}

	// particles
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].fitness; 
		aux = 0.;

		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].v[j] * k_aux[j];
		}

		if(aux != 0) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].fitnessAPM = function;
	}

	// pBest
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].pBest.fitness; 
		aux = 0.;
		 
		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].pBest.v[j] * k_aux[j];
		}

		if(aux != 0.) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].pBest.fitnessAPM = function;
		  
	}

	// gBest
	function = swarm->gBest.fitness; 
	aux = 0.;
		 
	for (j = 0; j < N_CON; j++) {
		aux += swarm->gBest.v[j] * k_aux[j];
	}

	if(aux != 0) {
		if(function < functionAvg) 
			function = functionAvg;	
	}
	function += aux;
	swarm->gBest.fitnessAPM = function;	
}

/**
 * Monotonic Sporadic APM function - penalty method
 */
void monotonic_sporadic_APM (swarm_t * swarm, double functionAvg, double * k_aux,  double ** violationAcum, int generation, double valuesArray[MAX_PAR][N_CON+1]) {
	int i, j, numViolation[N_CON], frequency = MAX_TIME / 10;
	double violation[N_CON], sumConst[MAX_PAR], k[N_CON], functionAvgK, function, aux;
	
	if(functionAvg > 0) 
		functionAvgK = functionAvg;
	else 
		functionAvgK = -functionAvg;

	for (i = 0; i < N_CON; i++) {
		violation[i] = 0.; 
		numViolation[i] = 0;
	}
	
	constraint (swarm, violationAcum, violation, sumConst, numViolation, valuesArray);	

	for (i = 0; i < N_CON; i++) 
			violation[i] = violation[i] / MAX_PAR;

	// Calculate k
	double quadrado = 0.;	
	for (i = 0; i < N_CON; i++) {
		quadrado += pow(violation[i], 2);
	}

	if(generation % frequency == 0) {
	
		for (i = 0; i < N_CON; i++) {
			if (quadrado != 0.)	
					k[i] = functionAvgK * violation[i] / quadrado;
				else
					k[i] = 0.;
		
			if(k[i] > k_aux[i])
				k_aux[i] = k[i];
		}
	}

	// particles
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].fitness; 
		aux = 0.;

		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].v[j] * k_aux[j];
		}

		if(aux != 0) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].fitnessAPM = function;
	}

	// pBest
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].pBest.fitness; 
		aux = 0.;
		 
		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].pBest.v[j] * k_aux[j];
		}

		if(aux != 0.) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].pBest.fitnessAPM = function;
		  
	}

	// gBest
	function = swarm->gBest.fitness; 
	aux = 0.;
		 
	for (j = 0; j < N_CON; j++) {
		aux += swarm->gBest.v[j] * k_aux[j];
	}

	if(aux != 0) {
		if(function < functionAvg) 
			function = functionAvg;	
	}
	function += aux;
	swarm->gBest.fitnessAPM = function;	
}

/**
 * Damping APM function - penalty method
 */
void damping_APM (swarm_t * swarm, double functionAvg, double * k_aux, double ** violationAcum, double valuesArray[MAX_PAR][N_CON+1]) {
	int i, j, numViolation[N_CON], frequency = MAX_TIME / 10;
	double violation[N_CON], sumConst[MAX_PAR], k[N_CON], functionAvgK, function, aux, theta = 0.5;

	if(functionAvg > 0) 
		functionAvgK = functionAvg;
	else 
		functionAvgK = -functionAvg;

	for (i = 0; i < N_CON; i++) {
		violation[i] = 0.; 
		numViolation[i] = 0;
	}
	
	constraint (swarm, violationAcum, violation, sumConst, numViolation, valuesArray);	

	for (i = 0; i < N_CON; i++) 
			violation[i] = violation[i] / MAX_PAR;

	// Calculate k
	double quadrado = 0.;	
	for (i = 0; i < N_CON; i++) {
		quadrado += pow(violation[i], 2);
	}

	for (i = 0; i < N_CON; i++) {
		if (quadrado != 0.)	
			k[i] = functionAvgK * violation[i] / quadrado;
		else
			k[i] = 0.;

		k_aux[i] = theta * k[i] + (1 - theta) * k_aux[i];

	}
	
	// particles
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].fitness; 
		aux = 0.;

		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].v[j] * k_aux[j];
		}

		if(aux != 0) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].fitnessAPM = function;
	}

	// pBest
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].pBest.fitness; 
		aux = 0.;
		 
		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].pBest.v[j] * k_aux[j];
		}

		if(aux != 0.) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].pBest.fitnessAPM = function;
		  
	}

	// gBest
	function = swarm->gBest.fitness; 
	aux = 0.;
		 
	for (j = 0; j < N_CON; j++) {
		aux += swarm->gBest.v[j] * k_aux[j];
	}

	if(aux != 0) {
		if(function < functionAvg) 
			function = functionAvg;	
	}
	function += aux;
	swarm->gBest.fitnessAPM = function;	
}

/**
 * worst APM function - penalty method
 */
void worst_APM (swarm_t * swarm, double functionAvg, double * k, double ** violationAcum, double valuesArray[MAX_PAR][N_CON+1]) {
	int i, j, numViolation[N_CON];
	double violation[N_CON], sumConst[MAX_PAR], functionAvgK, function, aux, worst, s;

	if(functionAvg > 0) 
		functionAvgK = functionAvg;
	else 
		functionAvgK = -functionAvg;

	for (i = 0; i < N_CON; i++) {
		violation[i] = 0.; 
		numViolation[i] = 0;
	}
	
	constraint (swarm, violationAcum, violation, sumConst, numViolation, valuesArray);	

	worst = functionAvg;
	for (i = 0; i < MAX_PAR; i++) {
		s = 0.;
		for (j = 0; j < N_CON; j++)
			s+= swarm->particles[i].v[j]; 
		if((s == 0.) && (worst < swarm->particles[i].fitnessAPM)) 
			worst = swarm->particles[i].fitnessAPM;
	}

	for (i = 0; i < N_CON; i++) 
		violation[i] = violation[i] / MAX_PAR;

	// Calculate k
	double quadrado = 0.;	
	for (i = 0; i < N_CON; i++) {
		quadrado += pow(violation[i], 2);
	}

	for (i = 0; i < N_CON; i++) {
		if (quadrado != 0.)	
			k[i] = functionAvgK * violation[i] / quadrado;
		else
			k[i] = 0.;
	}

	// particles
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].fitness; 
		aux = 0.;

		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].v[j] * k[j];
		}

		if(aux != 0) {
			if(function < worst) 
				function = worst;	
		}
		function += aux;
		swarm->particles[i].fitnessAPM = function;
	}
	// pBest
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].pBest.fitness; 
		aux = 0.;
		 
		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].pBest.v[j] * k[j];
		}

		if(aux != 0.) {
			if(function < worst) 
				function = worst;	
		}
		function += aux;
		swarm->particles[i].pBest.fitnessAPM = function;
		  
	}

	// gBest
	function = swarm->gBest.fitness; 
	aux = 0.;
		 
	for (j = 0; j < N_CON; j++) {
		aux += swarm->gBest.v[j] * k[j];
	}

	if(aux != 0) {
		if(function < worst) 
			function = worst;	
	}
	function += aux;
	swarm->gBest.fitnessAPM = function;
}

/**
 * worst 2 APM function - penalty method
 */
void worst_2_APM (swarm_t * swarm, double functionAvg, double * k, double ** violationAcum, double valuesArray[MAX_PAR][N_CON+1]) {
	int i, j, numViolation[N_CON];
	double violation[N_CON], sumConst[MAX_PAR], functionAvgK, function, aux, worst, s;

	if(functionAvg > 0) 
		functionAvgK = functionAvg;
	else 
		functionAvgK = -functionAvg;

	for (i = 0; i < N_CON; i++) {
		violation[i] = 0.; 
		numViolation[i] = 0;
	}
	
	constraint (swarm, violationAcum, violation, sumConst, numViolation, valuesArray);	

	worst = functionAvg;
	for (i = 0; i < MAX_PAR; i++) {
		s = 0.;
		for (j = 0; j < N_CON; j++)
			s+= swarm->particles[i].v[j]; 
		if((s == 0.) && (worst < swarm->particles[i].fitnessAPM)) 
			worst = swarm->particles[i].fitnessAPM;
	}

	for (i = 0; i < N_CON; i++) 
		violation[i] = violation[i] / MAX_PAR;

	// Calculate k
	double quadrado = 0.;	
	for (i = 0; i < N_CON; i++) {
		quadrado += pow(violation[i], 2);
	}

	for (i = 0; i < N_CON; i++) {
		if (quadrado != 0.)	
			k[i] = worst * violation[i] / quadrado;
		else
			k[i] = 0.;
	}

	// particles
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].fitness; 
		aux = 0.;

		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].v[j] * k[j];
		}

		if(aux != 0) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].fitnessAPM = function;
	}
	// pBest
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].pBest.fitness; 
		aux = 0.;
		 
		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].pBest.v[j] * k[j];
		}

		if(aux != 0.) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].pBest.fitnessAPM = function;
		  
	}

	// gBest
	function = swarm->gBest.fitness; 
	aux = 0.;
		 
	for (j = 0; j < N_CON; j++) {
		aux += swarm->gBest.v[j] * k[j];
	}

	if(aux != 0) {
		if(function < functionAvg) 
			function = functionAvg;	
	}
	function += aux;
	swarm->gBest.fitnessAPM = function;
}

/**
 * worst 3 APM function - penalty method
 */
void worst_3_APM (swarm_t * swarm, double functionAvg, double * k, double ** violationAcum, double valuesArray[MAX_PAR][N_CON+1]) {
	int i, j, numViolation[N_CON];
	double violation[N_CON], sumConst[MAX_PAR], functionAvgK, function, aux, worst, s;

	if(functionAvg > 0) 
		functionAvgK = functionAvg;
	else 
		functionAvgK = -functionAvg;

	for (i = 0; i < N_CON; i++) {
		violation[i] = 0.; 
		numViolation[i] = 0;
	}
	
	constraint (swarm, violationAcum, violation, sumConst, numViolation, valuesArray);	

	worst = functionAvg;
	for (i = 0; i < MAX_PAR; i++) {
		s = 0.;
		for (j = 0; j < N_CON; j++)
			s+= swarm->particles[i].v[j]; 
		if((s == 0.) && (worst < swarm->particles[i].fitnessAPM)) 
			worst = swarm->particles[i].fitnessAPM;
	}

	for (i = 0; i < N_CON; i++) 
		violation[i] = violation[i] / MAX_PAR;

	// Calculate k
	double quadrado = 0.;	
	for (i = 0; i < N_CON; i++) {
		quadrado += pow(violation[i], 2);
	}

	for (i = 0; i < N_CON; i++) {
		if (quadrado != 0.)	
			k[i] = worst * violation[i] / quadrado;
		else
			k[i] = 0.;
	}

	// particles
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].fitness; 
		aux = 0.;

		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].v[j] * k[j];
		}

		if(aux != 0) {
			if(function < worst) 
				function = worst;	
		}
		function += aux;
		swarm->particles[i].fitnessAPM = function;
	}
	// pBest
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].pBest.fitness; 
		aux = 0.;
		 
		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].pBest.v[j] * k[j];
		}

		if(aux != 0.) {
			if(function < worst) 
				function = worst;	
		}
		function += aux;
		swarm->particles[i].pBest.fitnessAPM = function;
		  
	}

	// gBest
	function = swarm->gBest.fitness; 
	aux = 0.;
		 
	for (j = 0; j < N_CON; j++) {
		aux += swarm->gBest.v[j] * k[j];
	}

	if(aux != 0) {
		if(function < worst) 
			function = worst;	
	}
	function += aux;
	swarm->gBest.fitnessAPM = function;
}

/**
 * med APM function - penalty method
 */
void med_APM (swarm_t * swarm, double functionAvg, double * k, double ** violationAcum, double valuesArray[MAX_PAR][N_CON+1]) {
	int i, j, numViolation[N_CON];
	double violation[N_CON], sumConst[MAX_PAR], functionAvgK, function, aux, worst, s;

	if(functionAvg > 0) 
		functionAvgK = functionAvg;
	else 
		functionAvgK = -functionAvg;

	for (i = 0; i < N_CON; i++) {
		violation[i] = 0.; 
		numViolation[i] = 0;
	}
	
	constraint (swarm, violationAcum, violation, sumConst, numViolation, valuesArray);	

	for (i = 0; i < N_CON; i++) {
		if (!numViolation[i])
			numViolation[i]++;	
		violation[i] = violation[i] / numViolation[i];
	}

	// Calculate k
	double quadrado = 0.;	
	for (i = 0; i < N_CON; i++) {
		quadrado += pow(violation[i], 2);
	}

	for (i = 0; i < N_CON; i++) {
		if (quadrado != 0.)	
			k[i] = functionAvgK * violation[i] / quadrado;
		else
			k[i] = 0.;
	}

	// particles
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].fitness; 
		aux = 0.;

		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].v[j] * k[j];
		}

		if(aux != 0) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].fitnessAPM = function;
	}
	// pBest
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].pBest.fitness; 
		aux = 0.;
		 
		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].pBest.v[j] * k[j];
		}

		if(aux != 0.) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].pBest.fitnessAPM = function;
		  
	}

	// gBest
	function = swarm->gBest.fitness; 
	aux = 0.;
		 
	for (j = 0; j < N_CON; j++) {
		aux += swarm->gBest.v[j] * k[j];
	}

	if(aux != 0) {
		if(function < functionAvg) 
			function = functionAvg;	
	}
	function += aux;
	swarm->gBest.fitnessAPM = function;
}


/**
 * med 2 APM function - penalty method
 */
void med_2_APM (swarm_t * swarm, double functionAvg, double * k, double ** violationAcum, double valuesArray[MAX_PAR][N_CON+1]) {
	int i, j, numViolation[N_CON], ind;
	double violation[N_CON], sumConst[MAX_PAR], functionAvgK, function, aux, worst, s, sumViolation;

	if(functionAvg > 0) 
		functionAvgK = functionAvg;
	else 
		functionAvgK = -functionAvg;

	for (i = 0; i < N_CON; i++) {
		violation[i] = 0.; 
		numViolation[i] = 0;
	}
	
	constraint (swarm, violationAcum, violation, sumConst, numViolation, valuesArray);	

	ind = 0;
	for (i = 0; i < MAX_PAR; i++) {
		sumViolation = 0.;
		for (j = 0; j < N_CON; j++) {
			sumViolation += swarm->particles[i].v[j];
		}
		if (sumViolation) 
			ind++;
	}	

	if (!ind)
		ind++;
	functionAvg  = functionAvg * MAX_PAR / ind;
	functionAvgK = functionAvgK * MAX_PAR / ind;

	for (i = 0; i < N_CON; i++) {
		if (!numViolation[i])
			numViolation[i]++;	
		violation[i] = violation[i] / numViolation[i];
	}

	// Calculate k
	double quadrado = 0.;	
	for (i = 0; i < N_CON; i++) {
		quadrado += pow(violation[i], 2);
	}

	for (i = 0; i < N_CON; i++) {
		if (quadrado != 0.)	
			k[i] = functionAvgK * violation[i] / quadrado;
		else
			k[i] = 0.;
	}

	// particles
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].fitness; 
		aux = 0.;

		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].v[j] * k[j];
		}

		if(aux != 0) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].fitnessAPM = function;
	}
	// pBest
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].pBest.fitness; 
		aux = 0.;
		 
		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].pBest.v[j] * k[j];
		}

		if(aux != 0.) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].pBest.fitnessAPM = function;
		  
	}

	// gBest
	function = swarm->gBest.fitness; 
	aux = 0.;
		 
	for (j = 0; j < N_CON; j++) {
		aux += swarm->gBest.v[j] * k[j];
	}

	if(aux != 0) {
		if(function < functionAvg) 
			function = functionAvg;	
	}
	function += aux;
	swarm->gBest.fitnessAPM = function;
}

/**
 * med 3 APM function - penalty method
 */
void med_3_APM (swarm_t * swarm, double functionAvg, double * k, double ** violationAcum, double valuesArray[MAX_PAR][N_CON+1]) {
	int i, j, numViolation[N_CON], ind;
	double violation[N_CON], sumConst[MAX_PAR], functionAvgK, function, aux, worst, s, sumViolation;

	if(functionAvg > 0) 
		functionAvgK = functionAvg;
	else 
		functionAvgK = -functionAvg;

	for (i = 0; i < N_CON; i++) {
		violation[i] = 0.; 
		numViolation[i] = 0;
	}
	
	constraint (swarm, violationAcum, violation, sumConst, numViolation, valuesArray);	

	ind = 0;
	for (i = 0; i < MAX_PAR; i++) {
		sumViolation = 0.;
		for (j = 0; j < N_CON; j++) {
			sumViolation += swarm->particles[i].v[j];
		}
		if (sumViolation) 
			ind++;
	}	

	if (!ind)
		ind++;
	functionAvgK = functionAvgK * MAX_PAR / ind;

	for (i = 0; i < N_CON; i++) {
		if (!numViolation[i])
			numViolation[i]++;	
		violation[i] = violation[i] / numViolation[i];
	}

	// Calculate k
	double quadrado = 0.;	
	for (i = 0; i < N_CON; i++) {
		quadrado += pow(violation[i], 2);
	}

	for (i = 0; i < N_CON; i++) {
		if (quadrado != 0.)	
			k[i] = functionAvgK * violation[i] / quadrado;
		else
			k[i] = 0.;
	}

	// particles
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].fitness; 
		aux = 0.;

		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].v[j] * k[j];
		}

		if(aux != 0) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].fitnessAPM = function;
	}
	// pBest
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].pBest.fitness; 
		aux = 0.;
		 
		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].pBest.v[j] * k[j];
		}

		if(aux != 0.) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].pBest.fitnessAPM = function;
		  
	}

	// gBest
	function = swarm->gBest.fitness; 
	aux = 0.;
		 
	for (j = 0; j < N_CON; j++) {
		aux += swarm->gBest.v[j] * k[j];
	}

	if(aux != 0) {
		if(function < functionAvg) 
			function = functionAvg;	
	}
	function += aux;
	swarm->gBest.fitnessAPM = function;
}

/**
 * med 4 APM function - penalty method
 */
void med_4_APM (swarm_t * swarm, double functionAvg, double * k, double ** violationAcum, double valuesArray[MAX_PAR][N_CON+1]) {
	int i, j, numViolation[N_CON], ind;
	double violation[N_CON], sumConst[MAX_PAR], functionAvgK, function, aux, worst, s, sumViolation;

	if(functionAvg > 0) 
		functionAvgK = functionAvg;
	else 
		functionAvgK = -functionAvg;

	for (i = 0; i < N_CON; i++) {
		violation[i] = 0.; 
		numViolation[i] = 0;
	}
	
	constraint (swarm, violationAcum, violation, sumConst, numViolation, valuesArray);	

	ind = 0;
	for (i = 0; i < MAX_PAR; i++) {
		sumViolation = 0.;
		for (j = 0; j < N_CON; j++) {
			sumViolation += swarm->particles[i].v[j];
		}
		if (sumViolation) 
			ind++;
	}	

	if (!ind)
		ind++;
	functionAvg  = functionAvg * MAX_PAR / ind;
	
	for (i = 0; i < N_CON; i++) {
		if (!numViolation[i])
			numViolation[i]++;	
		violation[i] = violation[i] / numViolation[i];
	}

	// Calculate k
	double quadrado = 0.;	
	for (i = 0; i < N_CON; i++) {
		quadrado += pow(violation[i], 2);
	}

	for (i = 0; i < N_CON; i++) {
		if (quadrado != 0.)	
			k[i] = functionAvgK * violation[i] / quadrado;
		else
			k[i] = 0.;
	}

	// particles
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].fitness; 
		aux = 0.;

		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].v[j] * k[j];
		}

		if(aux != 0) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].fitnessAPM = function;
	}
	// pBest
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].pBest.fitness; 
		aux = 0.;
		 
		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].pBest.v[j] * k[j];
		}

		if(aux != 0.) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].pBest.fitnessAPM = function;
		  
	}

	// gBest
	function = swarm->gBest.fitness; 
	aux = 0.;
		 
	for (j = 0; j < N_CON; j++) {
		aux += swarm->gBest.v[j] * k[j];
	}

	if(aux != 0) {
		if(function < functionAvg) 
			function = functionAvg;	
	}
	function += aux;
	swarm->gBest.fitnessAPM = function;
}

/**
 * med 5 APM function - penalty method
 */
void med_5_APM (swarm_t * swarm, double functionAvg, double * k, double ** violationAcum, double valuesArray[MAX_PAR][N_CON+1]) {
	int i, j, numViolation[N_CON], ind;
	double violation[N_CON], sumConst[MAX_PAR], functionAvgK, function, aux, worst, s, sumViolation;

	if(functionAvg > 0) 
		functionAvgK = functionAvg;
	else 
		functionAvgK = -functionAvg;

	for (i = 0; i < N_CON; i++) {
		violation[i] = 0.; 
		numViolation[i] = 0;
	}
	
	constraint (swarm, violationAcum, violation, sumConst, numViolation, valuesArray);	

	ind = 0;
	for (i = 0; i < MAX_PAR; i++) {
		sumViolation = 0.;
		for (j = 0; j < N_CON; j++) {
			sumViolation += swarm->particles[i].v[j];
		}
		if (sumViolation) 
			ind++;
	}	

	if (!ind)
		ind++;
	functionAvg  = functionAvg * MAX_PAR / ind;
	
	for (i = 0; i < N_CON; i++)
		violation[i] = violation[i] / MAX_PAR;

	// Calculate k
	double quadrado = 0.;	
	for (i = 0; i < N_CON; i++) {
		quadrado += pow(violation[i], 2);
	}

	for (i = 0; i < N_CON; i++) {
		if (quadrado != 0.)	
			k[i] = functionAvgK * violation[i] / quadrado;
		else
			k[i] = 0.;
	}

	// particles
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].fitness; 
		aux = 0.;

		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].v[j] * k[j];
		}

		if(aux != 0) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].fitnessAPM = function;
	}
	// pBest
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].pBest.fitness; 
		aux = 0.;
		 
		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].pBest.v[j] * k[j];
		}

		if(aux != 0.) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].pBest.fitnessAPM = function;
		  
	}

	// gBest
	function = swarm->gBest.fitness; 
	aux = 0.;
		 
	for (j = 0; j < N_CON; j++) {
		aux += swarm->gBest.v[j] * k[j];
	}

	if(aux != 0) {
		if(function < functionAvg) 
			function = functionAvg;	
	}
	function += aux;
	swarm->gBest.fitnessAPM = function;
}

/**
 * med 6 APM function - penalty method
 */
void med_6_APM (swarm_t * swarm, double functionAvg, double * k, double ** violationAcum, double valuesArray[MAX_PAR][N_CON+1]) {
	int i, j, numViolation[N_CON], ind;
	double violation[N_CON], sumConst[MAX_PAR], functionAvgK, function, aux, worst, s, sumViolation;

	if(functionAvg > 0) 
		functionAvgK = functionAvg;
	else 
		functionAvgK = -functionAvg;

	for (i = 0; i < N_CON; i++) {
		violation[i] = 0.; 
		numViolation[i] = 0;
	}
	
	constraint (swarm, violationAcum, violation, sumConst, numViolation, valuesArray);	

	ind = 0;
	for (i = 0; i < MAX_PAR; i++) {
		sumViolation = 0.;
		for (j = 0; j < N_CON; j++) {
			sumViolation += swarm->particles[i].v[j];
		}
		if (sumViolation) 
			ind++;
	}	

	if (!ind)
		ind++;
	functionAvgK = functionAvgK * MAX_PAR / ind;	
	
	for (i = 0; i < N_CON; i++)
		violation[i] = violation[i] / MAX_PAR;

	// Calculate k
	double quadrado = 0.;	
	for (i = 0; i < N_CON; i++) {
		quadrado += pow(violation[i], 2);
	}

	for (i = 0; i < N_CON; i++) {
		if (quadrado != 0.)	
			k[i] = functionAvgK * violation[i] / quadrado;
		else
			k[i] = 0.;
	}

	// particles
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].fitness; 
		aux = 0.;

		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].v[j] * k[j];
		}

		if(aux != 0) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].fitnessAPM = function;
	}
	// pBest
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].pBest.fitness; 
		aux = 0.;
		 
		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].pBest.v[j] * k[j];
		}

		if(aux != 0.) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].pBest.fitnessAPM = function;
		  
	}

	// gBest
	function = swarm->gBest.fitness; 
	aux = 0.;
		 
	for (j = 0; j < N_CON; j++) {
		aux += swarm->gBest.v[j] * k[j];
	}

	if(aux != 0) {
		if(function < functionAvg) 
			function = functionAvg;	
	}
	function += aux;
	swarm->gBest.fitnessAPM = function;
}


/**
 * med 7 APM function - penalty method
 */
void med_7_APM (swarm_t * swarm, double functionAvg, double * k, double ** violationAcum, double valuesArray[MAX_PAR][N_CON+1]) {
	int i, j, numViolation[N_CON], ind;
	double violation[N_CON], sumConst[MAX_PAR], functionAvgK, function, aux, worst, s, sumViolation;

	if(functionAvg > 0) 
		functionAvgK = functionAvg;
	else 
		functionAvgK = -functionAvg;

	for (i = 0; i < N_CON; i++) {
		violation[i] = 0.; 
		numViolation[i] = 0;
	}
	
	constraint (swarm, violationAcum, violation, sumConst, numViolation, valuesArray);	

	ind = 0;
	for (i = 0; i < MAX_PAR; i++) {
		sumViolation = 0.;
		for (j = 0; j < N_CON; j++) {
			sumViolation += swarm->particles[i].v[j];
		}
		if (sumViolation) 
			ind++;
	}	

	if (!ind)
		ind++;
	functionAvg  = functionAvg * MAX_PAR / ind;
	functionAvgK = functionAvgK * MAX_PAR / ind;	
	
	for (i = 0; i < N_CON; i++)
		violation[i] = violation[i] / MAX_PAR;

	// Calculate k
	double quadrado = 0.;	
	for (i = 0; i < N_CON; i++) {
		quadrado += pow(violation[i], 2);
	}

	for (i = 0; i < N_CON; i++) {
		if (quadrado != 0.)	
			k[i] = functionAvgK * violation[i] / quadrado;
		else
			k[i] = 0.;
	}

	// particles
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].fitness; 
		aux = 0.;

		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].v[j] * k[j];
		}

		if(aux != 0) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].fitnessAPM = function;
	}
	// pBest
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].pBest.fitness; 
		aux = 0.;
		 
		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].pBest.v[j] * k[j];
		}

		if(aux != 0.) {
			if(function < functionAvg) 
				function = functionAvg;	
		}
		function += aux;
		swarm->particles[i].pBest.fitnessAPM = function;
		  
	}

	// gBest
	function = swarm->gBest.fitness; 
	aux = 0.;
		 
	for (j = 0; j < N_CON; j++) {
		aux += swarm->gBest.v[j] * k[j];
	}

	if(aux != 0) {
		if(function < functionAvg) 
			function = functionAvg;	
	}
	function += aux;
	swarm->gBest.fitnessAPM = function;
}

/**
 * med worst APM function - penalty method
 */
void med_worst_APM (swarm_t * swarm, double functionAvg, double * k, double ** violationAcum, double valuesArray[MAX_PAR][N_CON+1]) {
	int i, j, numViolation[N_CON], ind;
	double violation[N_CON], sumConst[MAX_PAR], functionAvgK, function, aux, worst, s, sumViolation;

	if(functionAvg > 0) 
		functionAvgK = functionAvg;
	else 
		functionAvgK = -functionAvg;

	for (i = 0; i < N_CON; i++) {
		violation[i] = 0.; 
		numViolation[i] = 0;
	}
	
	constraint (swarm, violationAcum, violation, sumConst, numViolation, valuesArray);	

	worst = functionAvg;
	for (i = 0; i < MAX_PAR; i++) {
		s = 0.;
		for (j = 0; j < N_CON; j++)
			s+= swarm->particles[i].v[j]; 
		if((s == 0.) && (worst < swarm->particles[i].fitnessAPM)) 
			worst = swarm->particles[i].fitnessAPM;
	}
		
	for (i = 0; i < N_CON; i++) {
		if (!numViolation[i])
			numViolation[i]++;	
		violation[i] = violation[i] / numViolation[i];
	}

	// Calculate k
	double quadrado = 0.;	
	for (i = 0; i < N_CON; i++) {
		quadrado += pow(violation[i], 2);
	}

	for (i = 0; i < N_CON; i++) {
		if (quadrado != 0.)	
			k[i] = functionAvgK * violation[i] / quadrado;
		else
			k[i] = 0.;
	}

	// particles
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].fitness; 
		aux = 0.;

		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].v[j] * k[j];
		}

		if(aux != 0) {
			if(function < worst) 
				function = worst;	
		}
		function += aux;
		swarm->particles[i].fitnessAPM = function;
	}
	// pBest
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].pBest.fitness; 
		aux = 0.;
		 
		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].pBest.v[j] * k[j];
		}

		if(aux != 0.) {
			if(function < worst) 
				function = worst;	
		}
		function += aux;
		swarm->particles[i].pBest.fitnessAPM = function;
		  
	}

	// gBest
	function = swarm->gBest.fitness; 
	aux = 0.;
		 
	for (j = 0; j < N_CON; j++) {
		aux += swarm->gBest.v[j] * k[j];
	}

	if(aux != 0) {
		if(function < worst) 
			function = worst;	
	}
	function += aux;
	swarm->gBest.fitnessAPM = function;
}

/**
 * med worst 2 APM function - penalty method
 */
void med_worst_2_APM (swarm_t * swarm, double functionAvg, double * k, double ** violationAcum, double valuesArray[MAX_PAR][N_CON+1]) {
	int i, j, numViolation[N_CON], ind;
	double violation[N_CON], sumConst[MAX_PAR], functionAvgK, function, aux, worst, s, sumViolation;

	if(functionAvg > 0) 
		functionAvgK = functionAvg;
	else 
		functionAvgK = -functionAvg;

	for (i = 0; i < N_CON; i++) {
		violation[i] = 0.; 
		numViolation[i] = 0;
	}
	
	constraint (swarm, violationAcum, violation, sumConst, numViolation, valuesArray);	

	worst = functionAvg;
	for (i = 0; i < MAX_PAR; i++) {
		s = 0.;
		for (j = 0; j < N_CON; j++)
			s+= swarm->particles[i].v[j]; 
		if((s == 0.) && (worst < swarm->particles[i].fitnessAPM)) 
			worst = swarm->particles[i].fitnessAPM;
	}
		
	for (i = 0; i < N_CON; i++) {
		if (!numViolation[i])
			numViolation[i]++;	
		violation[i] = violation[i] / numViolation[i];
	}

	// Calculate k
	double quadrado = 0.;	
	for (i = 0; i < N_CON; i++) {
		quadrado += pow(violation[i], 2);
	}

	for (i = 0; i < N_CON; i++) {
		if (quadrado != 0.)	
			k[i] = worst * violation[i] / quadrado;
		else
			k[i] = 0.;
	}

	// particles
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].fitness; 
		aux = 0.;

		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].v[j] * k[j];
		}

		if(aux != 0) {
			if(function < worst) 
				function = worst;	
		}
		function += aux;
		swarm->particles[i].fitnessAPM = function;
	}
	// pBest
	for (i = 0; i < MAX_PAR; i++) {
		function = swarm->particles[i].pBest.fitness; 
		aux = 0.;
		 
		for (j = 0; j < N_CON; j++) {
			aux += swarm->particles[i].pBest.v[j] * k[j];
		}

		if(aux != 0.) {
			if(function < worst) 
				function = worst;	
		}
		function += aux;
		swarm->particles[i].pBest.fitnessAPM = function;
		  
	}

	// gBest
	function = swarm->gBest.fitness; 
	aux = 0.;
		 
	for (j = 0; j < N_CON; j++) {
		aux += swarm->gBest.v[j] * k[j];
	}

	if(aux != 0) {
		if(function < worst) 
			function = worst;	
	}	
	function += aux;
	swarm->gBest.fitnessAPM = function;
}

/**
 * Fitness function - executes the objective function.
 */
void fitness (swarm_t * swarm, double * k, double ** violationAcum, int generation, F101Truss10Bar * truss) {

	int i, j, m, cont1, cont2, cont3;
	double x[N_VAR], functionAvg = 0., function = 0., length[N_VAR];
	double xTEST[truss->getDimension()];

	// PEDRO
	// if (truss->getNumberConstraints() != N_CON){
	// 	printf("Error . Truss constraints differs from N_CON. N_CON: %d\t truss.getNumberConstraints(): %d\n", N_CON, truss->getNumberConstraints());
	// 	printf("N_VAR: %d\tTruss dimension: %d\n", N_VAR, truss->getDimension());
	// 	exit (1);
	// }
	// double valuesArray[MAX_PAR][truss->getNumberConstraints() + 1]; // objectiveFunction size + number of constraints
	double valuesArray[MAX_PAR][N_CON + 1]; // objectiveFunction size + number of constraints
	// printf("declarou, number constraints: %d \n", truss->getNumberConstraints());
	for (i = 0; i < MAX_PAR; i++) {
		for (j = 0; j < N_VAR; j++) {
			x[j] = swarm->particles[i].position[j];
			//printf("x(fit) = %lf  ", x[j]);
		}
		//printf("\n");
		
		// printf("Array de x em fitness\n");
		// for ( j = 0; j< N_VAR; j++){
		// 	printf("%f ", x[j]);
		// }
		// printf("\n");

		// Functions
		#ifdef USE_STRING
			/*x[0] = 11.852177;
			x[1] = 0.34747463;
			x[2] = 0.051301897;*/
			function = (x[0] + 2) * x[1] * x[2] * x[2];
		#endif

		#ifdef USE_REDUCER
			/*x[0] = 3.5;
			x[1] = 0.7;
			x[2] = 17;
			x[3] = 7.3000035;
			x[4] = 7.7115326;
			x[5] = 3.350216;
			x[6] = 5.286655;*/
			function = 0.7854 * x[0] * pow(x[1],2) * (3.3333 * pow(x[2],2) + 14.9334 * x[2] - 43.0934) - 1.508 * x[0] * (pow(x[5], 2) + pow(x[6], 2)) + 7.4777 * (pow(x[5], 3) + pow(x[6], 3)) + 0.7854 * (x[3] * pow(x[5], 2) + x[4] * pow(x[6], 2));
		#endif

		#ifdef USE_WELDED
			/*x[0] = 0.2442949; 
			x[1] = 6.2116738;
			x[2] = 8.3015486;
			x[3] = 0.2443003;*/
			function = (1.10471 * x[0] * x[0] * x[1]) + 0.04811 * x[2] * x[3] * (14.0 + x[1]);
		#endif

		#ifdef USE_WELDED_MOD1_7
			/*x[0] = 0.2442949; 
			x[1] = 6.2116738;
			x[2] = 8.3015486;
			x[3] = 0.2443003;*/
			function = (1.10471 * x[0] * x[0] * x[1]) + 0.04811 * x[2] * x[3] * (14.0 + x[1]);
		#endif

		#ifdef USE_PRESSURE
			double aux;
			aux = round(x[0]);
			x[0] = aux * 0.0625;
			aux = round(x[1]);
			x[1] = aux * 0.0625;
			/*x[0] = 0.8125;
			x[1] = 0.4375;
			x[2] = 42.086994;
			x[3] = 176.779128;*/
			function = (0.6224 * x[0] * x[2] * x[3]) + (1.7781 * x[1] * x[2] * x[2]) + (3.1661 * x[0] * x[0] * x[3]) + (19.84 * x[0] * x[0] * x[2]); 
		#endif
		
		#ifdef USE_CANTILEVER
			double B[] = {2.4, 2.6, 2.8, 3.1}; 
			double H[] = {45.0, 50.0, 55.0, 60.0};
			int aux, cont;
			x[0] = (int) x[0];
			x[1] = (int) x[1];
			aux = round(x[2]);
			x[2] = B[aux];
			aux = round(x[4]);
			x[4] = B[aux];
			aux = round(x[3]);
			x[3] = H[aux];
			aux = round(x[5]);
			x[5] = H[aux];

			/*x[0] = 3; 
			x[1] = 60;
			x[2] = 3.1;
			x[3] = 60;
			x[4] = 2.6;
			x[5] = 50;
			x[6] = 2.2094;
			x[7] = 44.0428;
			x[8] = 2.0944;
			x[9] = 31.9867;*/

			function = 0.;
			for (cont = 0; cont < 5; cont++ ){
				function += x[cont * 2] * x[cont * 2 + 1];
			}
			function *= 100;
		#endif

		#ifdef USE_T10C
			double length[N_VAR];
			length[0] = length[1] = length[2] = length[3] = length[4] = length[5] = 360.0 / 10;
			length[6] = length[7] = length[8] = length[9] = 509.116882454 / 10;  

			/*x[0] = 30.162525; 
			x[1] = 0.10003946;
			x[2] = 22.81192;
			x[3] = 15.871827;
			x[4] = 0.10000233;
			x[5] = 0.5149511;
			x[6] = 7.505953;
			x[7] = 21.264076;
			x[8] = 21.383036;
			x[9] = 0.10000795;*/

			
			function = 0.;
			for (j = 0; j < N_VAR; j++) {
				function += x[j] * length[j];
			}
			//function /= 10;

			// MODIFICADO POR PEDRO
			// Evaluate using Eureka
			truss->evaluation(x, valuesArray[i]);
			// printf("function (objFunc) Erica: %f\t function (objFunc) Eureka: %f\n", function, valuesArray[i][0]);			// function = valuesArray[0]; // Function gets the objective function value
			function = valuesArray[i][0]; // Avaliação do Eureka ao invés da Erica
			// printf("Constraints: \t");
			// for (j =0; j < truss->getNumberConstraints() + 2; ++j){
			// 	printf("%f ", valuesArray[i][j]);
			// 	// printf("%f ", valuesArrayTest[j]);
			// }
			// exit(4);

		#endif

		#ifdef USE_T10D
			length[0] = length[1] = length[2] = length[3] = length[4] = length[5] = 360.0 / 10;
			length[6] = length[7] = length[8] = length[9] = 509.116882454 / 10;  

			/*x[0] = 33.5; 
			x[1] = 1.62;
			x[2] = 22.9;
			x[3] = 14.2;
			x[4] = 1.62;
			x[5] = 1.62;
			x[6] = 7.97;
			x[7] = 22.9;
			x[8] = 22.0;
			x[9] = 1.62;*/

			x[0] = 30.999442;
			x[1] = 0.214494;
			x[2] = 28.288087;
			x[3] = 24.951834;
			x[4] = 0.001618;
			x[5] = 0.440247;
			x[6] = 20.079916;
			x[7] = 27.867945;
			x[8] = 28.277859;
			x[9] = 0.000000;

			double tab[] = {1.62, 1.80, 1.99, 2.13, 2.38, 2.62, 2.93, 3.13, 3.38, 3.47, 3.55, 3.63, 3.88, 4.22, 4.49, 4.59, 4.80, 4.97, 5.12, 5.74, 7.97, 11.50, 13.50, 14.20, 15.50, 16.90, 18.80, 19.90, 22.00, 26.50, 30.00, 33.50};
			int pos;
			for (j = 0; j < N_VAR; j++) {
				pos = round(x[j]);
				x[j] = tab[pos];
			}

			function = 0.;
			for (j = 0; j < N_VAR; j++) {
				function += x[j] * length[j];
			}
			//function /= 10;

			// MODIFICADO POR PEDRO
			// Evaluate using Eureka

			printf("generation: %d", generation);
			int itP;
			printf("Printing xtests (design variables)\n");
			for(itP=0 ;itP<10 ;itP++){
				printf("x[%d] = %f\n", itP, x[itP]);
			}
			printf("function (objFunc) Erica: %f\t function (objFunc) Eureka: %f\n", function, valuesArray[i][0]);			// function = valuesArray[0]; // Function gets the objective function value
			exit(5);
			truss->evaluation(x, valuesArray[i]);
			// printf("function (objFunc) Erica: %f\t function (objFunc) Eureka: %f\n", function, valuesArray[i][0]);			// function = valuesArray[0]; // Function gets the objective function value
			function = valuesArray[i][0]; // Avaliação do Eureka ao invés da Erica
		#endif

		#ifdef USE_T25C
			
			length[0] = 75.000000; 
			length[1] = 130.503831; 
			length[2] = 130.503831;
			length[3] = 130.503831; 
			length[4] = 130.503831; 
			length[5] = 106.800047; 
			length[6] = 106.800047; 
			length[7] = 106.800047; 
			length[8] = 106.800047; 
			length[9] = 75.000000; 
			length[10] = 75.000000; 
			length[11] = 75.000000; 
			length[12] = 75.000000; 
			length[13] = 181.142209; 
			length[14] = 181.142209; 
			length[15] = 181.142209; 
			length[16] = 181.142209; 
			length[17] = 181.142209; 
			length[18] = 181.142209; 
			length[19] = 181.142209; 
			length[20] = 181.142209; 
			length[21] = 133.463478; 
			length[22] = 133.463478; 
			length[23] = 133.463478; 
			length[24] = 133.463478;

			swarm->particles[i].position[2] = swarm->particles[i].position[3] = swarm->particles[i].position[4] = swarm->particles[i].position[1];
			swarm->particles[i].position[6] = swarm->particles[i].position[7] = swarm->particles[i].position[8] = swarm->particles[i].position[5];
			swarm->particles[i].position[10] = swarm->particles[i].position[9];
			swarm->particles[i].position[12] = swarm->particles[i].position[11];
			swarm->particles[i].position[14] = swarm->particles[i].position[15] = swarm->particles[i].position[16] = swarm->particles[i].position[13];
			swarm->particles[i].position[18] = swarm->particles[i].position[19] = swarm->particles[i].position[20] = swarm->particles[i].position[17];
			swarm->particles[i].position[22] = swarm->particles[i].position[23] = swarm->particles[i].position[24] = swarm->particles[i].position[21];

			x[2] = x[3] = x[4] = x[1];
			x[6] = x[7] = x[8] = x[5];
			x[10] = x[9];
			x[12] = x[11];
			x[14] = x[15] = x[16] = x[13];
			x[18] = x[19] = x[20] = x[17];
			x[22] = x[23] = x[24] = x[21];

			function = 0.;
			
			for (j = 0; j < N_VAR; j++) {
				function += x[j] * length[j];
			}
			function /= 10;

			// MODIFICADO POR PEDRO
			// Evaluate using Eureka
			truss->evaluation(x, valuesArray[i]);
			// printf("function (objFunc) Erica: %f\t function (objFunc) Eureka: %f\n", function, valuesArray[i][0]);			// function = valuesArray[0]; // Function gets the objective function value
			function = valuesArray[i][0]; // Avaliação do Eureka ao invés da Erica
		#endif

		#ifdef USE_T25D
			
			// // MODIFICADO POR PEDRO
			// // Evaluate using Eureka
			// printf("X array: \t");
			// for (j = 0; j < truss->getDimension(); ++j){
			// 	printf("%f ", x[j]);
			// }
			// truss->evaluation(x, valuesArray[i]);
			// // printf("function (objFunc) Erica: %f\t function (objFunc) Eureka: %f\n", function, valuesArray[i][0]);			// function = valuesArray[0]; // Function gets the objective function value
			// // function = valuesArray[i][0]; // Avaliação do Eureka ao invés da Erica
			// printf("\n");
			// // printf("function (objFunc): %f\n",  valuesArray[i][0]);			// function = valuesArray[0]; // Function gets the objective function value
			// printf("Constraints: \t");
			// for (j =0; j < truss->getNumberConstraints() + 2; ++j){
			// 	printf("%f ", valuesArray[i][j]);
			// }
			// printf("\n");
			// // exit(3);

			length[0] = 75.000000; 
			length[1] = 130.503831; 
			length[2] = 130.503831;
			length[3] = 130.503831; 
			length[4] = 130.503831; 
			length[5] = 106.800047; 
			length[6] = 106.800047; 
			length[7] = 106.800047; 
			length[8] = 106.800047; 
			length[9] = 75.000000; 
			length[10] = 75.000000; 
			length[11] = 75.000000; 
			length[12] = 75.000000; 
			length[13] = 181.142209; 
			length[14] = 181.142209; 
			length[15] = 181.142209; 
			length[16] = 181.142209; 
			length[17] = 181.142209; 
			length[18] = 181.142209; 
			length[19] = 181.142209; 
			length[20] = 181.142209; 
			length[21] = 133.463478; 
			length[22] = 133.463478; 
			length[23] = 133.463478; 
			length[24] = 133.463478;



			swarm->particles[i].position[2] = swarm->particles[i].position[3] = swarm->particles[i].position[4] = swarm->particles[i].position[1];
			swarm->particles[i].position[6] = swarm->particles[i].position[7] = swarm->particles[i].position[8] = swarm->particles[i].position[5];
			swarm->particles[i].position[10] = swarm->particles[i].position[9];
			swarm->particles[i].position[12] = swarm->particles[i].position[11];
			swarm->particles[i].position[14] = swarm->particles[i].position[15] = swarm->particles[i].position[16] = swarm->particles[i].position[13];
			swarm->particles[i].position[18] = swarm->particles[i].position[19] = swarm->particles[i].position[20] = swarm->particles[i].position[17];
			swarm->particles[i].position[22] = swarm->particles[i].position[23] = swarm->particles[i].position[24] = swarm->particles[i].position[21];

			// x[0] = 0.043417;
			// x[1] = 2.417307;
			// x[2] = 2.417307;
			// x[3] = 2.417307;
			// x[4] = 2.417307;
			// x[5] = 28.999588;
			// x[6] = 28.999588;
			// x[7] = 28.999588;
			// x[8] = 28.999588;
			// x[9] = 0.068685;
			// x[10] = 0.068685;
			// x[11] = 19.736034;
			// x[12] = 19.736034;
			// x[13] = 9.427932;
			// x[14] = 9.427932;
			// x[15] = 9.427932;
			// x[16] = 9.427932;
			// x[17] = 4.331147;
			// x[18] = 4.331147;
			// x[19] = 4.331147;
			// x[20] = 4.331147;
			// x[21] = 29.000000;
			// x[22] = 29.000000;
			// x[23] = 29.000000;
			// x[24] = 29.000000;


			x[2] = x[3] = x[4] = x[1];
			x[6] = x[7] = x[8] = x[5];
			x[10] = x[9];
			x[12] = x[11];
			x[14] = x[15] = x[16] = x[13];
			x[18] = x[19] = x[20] = x[17];
			x[22] = x[23] = x[24] = x[21];

			xTEST[0] = x[0];
			xTEST[1] = x[4];
			xTEST[2] = x[8];
			xTEST[3] = x[10];
			xTEST[4] = x[12];
			xTEST[5] = x[16];
			xTEST[6] = x[20];
			xTEST[7] = x[24];
			
			
		// x[0] = 0.100000
		// x[1] = 0.300000
		// x[2] = 3.400000
		// x[3] = 0.100000
		// x[4] = 2.100000
		// x[5] = 1.000000
		// x[6] = 0.500000
		// x[7] = 3.400000



			double tab[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.8, 3.0, 3.2, 3.4};
			int pos;
			for (j = 0; j < N_VAR; j++) {
				pos = round(x[j]);
				x[j] = tab[pos];
				if (j < 8){ // number of constraints
					pos = round(xTEST[j]);
					xTEST[j] = tab[pos];
				}
			}
			// printf("X APÓS DISCRETIZAR: \t");
			// for (j = 0; j < truss->getDimension(); ++j){
			// 	printf("%f ", xTEST[j]);
			// }
			// printf("\n");

			truss->evaluation(xTEST, valuesArray[i]);



			function = 0.;
			for (j = 0; j < N_VAR; j++) {
				function += x[j] * length[j];
			}
			function /= 10;

			// // MODIFICADO POR PEDRO
			// // Evaluate using Eureka
			// printf("X array: \t");
			// for (j = 0; j < truss->getDimension(); ++j){
			// 	printf("%f ", x[j]);
			// }
			// printf("\n");
			// // truss->evaluation(x, valuesArray[i]);
			// printf("function (objFunc) Erica: %f\t function (objFunc) Eureka: %f\n", function, valuesArray[i][0]);			// function = valuesArray[0]; // Function gets the objective function value
			// // function = valuesArray[i][0]; // Avaliação do Eureka ao invés da Erica
			// printf("\n");
			// // printf("function (objFunc): %f\n",  valuesArray[i][0]);			// function = valuesArray[0]; // Function gets the objective function value
			// printf("Constraints: \t");
			// for (j =0; j < truss->getNumberConstraints() + 2; ++j){
			// 	printf("%f ", valuesArray[i][j]);
			// }
			// printf("\n");
			// // exit(3);

			function = valuesArray[i][0]; // Avaliação do Eureka ao invés da Erica
			// printf("generation: %d", generation);
			// // int itP;
			// printf("Printing xtests (design variables)\n");
			// for(itP=0 ;itP<8 ;itP++){
			// 	printf("x[%d] = %f\n", itP, xTEST[itP]);
			// }
			// printf("function (objFunc) Erica: %f\t function (objFunc) Eureka: %f\n", function, valuesArray[i][0]);			// function = valuesArray[0]; // Function gets the objective function value
			// exit(5);
			// printf("function (objFunc) Erica: %f\t function (objFunc) Eureka: %f\n", function, valuesArray[i][0]);			// function = valuesArray[0]; // Function gets the objective function value
			// exit(3);

		#endif

		#ifdef USE_T52C

			length[0] = 3000.000000/1000;
			length[1] = 3000.000000/1000;
			length[2] = 3000.000000/1000;
			length[3] = 3000.000000/1000;
			length[4] = 3605.551275/1000;
			length[5] = 3605.551275/1000;
			length[6] = 3605.551275/1000;
			length[7] = 3605.551275/1000;
			length[8] = 3605.551275/1000;
			length[9] = 3605.551275/1000;
			length[10] = 2000.000000/1000;
			length[11] = 2000.000000/1000;
			length[12] = 2000.000000/1000;
			length[13] = 3000.000000/1000;
			length[14] = 3000.000000/1000;
			length[15] = 3000.000000/1000;
			length[16] = 3000.000000/1000;
			length[17] = 3605.551275/1000;
			length[18] = 3605.551275/1000;
			length[19] = 3605.551275/1000;
			length[20] = 3605.551275/1000;
			length[21] = 3605.551275/1000;
			length[22] = 3605.551275/1000;
			length[23] = 2000.000000/1000;
			length[24] = 2000.000000/1000;
			length[25] = 2000.000000/1000;
			length[26] = 3000.000000/1000;
			length[27] = 3000.000000/1000;
			length[28] = 3000.000000/1000;
			length[29] = 3000.000000/1000;
			length[30] = 3605.551275/1000;
			length[31] = 3605.551275/1000;
			length[32] = 3605.551275/1000;
			length[33] = 3605.551275/1000;
			length[34] = 3605.551275/1000;
			length[35] = 3605.551275/1000;
			length[36] = 2000.000000/1000;
			length[37] = 2000.000000/1000;
			length[38] = 2000.000000/1000;
			length[39] = 3000.000000/1000;
			length[40] = 3000.000000/1000;
			length[41] = 3000.000000/1000;
			length[42] = 3000.000000/1000;
			length[43] = 3605.551275/1000;
			length[44] = 3605.551275/1000;
			length[45] = 3605.551275/1000;
			length[46] = 3605.551275/1000;
			length[47] = 3605.551275/1000;
			length[48] = 3605.551275/1000;
			length[49] = 2000.000000/1000;
			length[50] = 2000.000000/1000;
			length[51] = 2000.000000/1000;

			swarm->particles[i].position[1] = swarm->particles[i].position[2] = swarm->particles[i].position[3] = swarm->particles[i].position[0];
			swarm->particles[i].position[5] = swarm->particles[i].position[6] = swarm->particles[i].position[7] = swarm->particles[i].position[8] = swarm->particles[i].position[9] = swarm->particles[i].position[4]; 
			swarm->particles[i].position[11] = swarm->particles[i].position[12] = swarm->particles[i].position[10];
			swarm->particles[i].position[14] = swarm->particles[i].position[15] = swarm->particles[i].position[16] = swarm->particles[i].position[13];
			swarm->particles[i].position[18] = swarm->particles[i].position[19] = swarm->particles[i].position[20] = swarm->particles[i].position[21] = swarm->particles[i].position[22] = swarm->particles[i].position[17];
			swarm->particles[i].position[24] = swarm->particles[i].position[25] = swarm->particles[i].position[23];
			swarm->particles[i].position[27] = swarm->particles[i].position[28] = swarm->particles[i].position[29] = swarm->particles[i].position[26];
			swarm->particles[i].position[31] = swarm->particles[i].position[32] = swarm->particles[i].position[33] = swarm->particles[i].position[34] = swarm->particles[i].position[35] = swarm->particles[i].position[30];
			swarm->particles[i].position[37] = swarm->particles[i].position[38] = swarm->particles[i].position[36];
			swarm->particles[i].position[40] = swarm->particles[i].position[41] = swarm->particles[i].position[42] = swarm->particles[i].position[39];
			swarm->particles[i].position[44] = swarm->particles[i].position[45] = swarm->particles[i].position[46] = swarm->particles[i].position[47] = swarm->particles[i].position[48] = swarm->particles[i].position[43];
			swarm->particles[i].position[50] = swarm->particles[i].position[51] = swarm->particles[i].position[49];

			x[1] = x[2] = x[3] = x[0];
			x[5] = x[6] = x[7] = x[8] = x[9] = x[4]; 
			x[11] = x[12] = x[10];
			x[14] = x[15] = x[16] = x[13];
			x[18] = x[19] = x[20] = x[21] = x[22] = x[17];
			x[24] = x[25] = x[23];
			x[27] = x[28] = x[29] = x[26];
			x[31] = x[32] = x[33] = x[34] = x[35] = x[30];
			x[37] = x[38] = x[36];
			x[40] = x[41] = x[42] = x[39];
			x[44] = x[45] = x[46] = x[47] = x[48] = x[43];
			x[50] = x[51] = x[49];

			function = 0.;
			
			for (j = 0; j < N_VAR; j++) {
				function += x[j] * length[j];
			}
			function /= 10;
		#endif

		#ifdef USE_T52D
			
			length[0] = 3000.000000/1000;
			length[1] = 3000.000000/1000;
			length[2] = 3000.000000/1000;
			length[3] = 3000.000000/1000;
			length[4] = 3605.551275/1000;
			length[5] = 3605.551275/1000;
			length[6] = 3605.551275/1000;
			length[7] = 3605.551275/1000;
			length[8] = 3605.551275/1000;
			length[9] = 3605.551275/1000;
			length[10] = 2000.000000/1000;
			length[11] = 2000.000000/1000;
			length[12] = 2000.000000/1000;
			length[13] = 3000.000000/1000;
			length[14] = 3000.000000/1000;
			length[15] = 3000.000000/1000;
			length[16] = 3000.000000/1000;
			length[17] = 3605.551275/1000;
			length[18] = 3605.551275/1000;
			length[19] = 3605.551275/1000;
			length[20] = 3605.551275/1000;
			length[21] = 3605.551275/1000;
			length[22] = 3605.551275/1000;
			length[23] = 2000.000000/1000;
			length[24] = 2000.000000/1000;
			length[25] = 2000.000000/1000;
			length[26] = 3000.000000/1000;
			length[27] = 3000.000000/1000;
			length[28] = 3000.000000/1000;
			length[29] = 3000.000000/1000;
			length[30] = 3605.551275/1000;
			length[31] = 3605.551275/1000;
			length[32] = 3605.551275/1000;
			length[33] = 3605.551275/1000;
			length[34] = 3605.551275/1000;
			length[35] = 3605.551275/1000;
			length[36] = 2000.000000/1000;
			length[37] = 2000.000000/1000;
			length[38] = 2000.000000/1000;
			length[39] = 3000.000000/1000;
			length[40] = 3000.000000/1000;
			length[41] = 3000.000000/1000;
			length[42] = 3000.000000/1000;
			length[43] = 3605.551275/1000;
			length[44] = 3605.551275/1000;
			length[45] = 3605.551275/1000;
			length[46] = 3605.551275/1000;
			length[47] = 3605.551275/1000;
			length[48] = 3605.551275/1000;
			length[49] = 2000.000000/1000;
			length[50] = 2000.000000/1000;
			length[51] = 2000.000000/1000;

			swarm->particles[i].position[1] = swarm->particles[i].position[2] = swarm->particles[i].position[3] = swarm->particles[i].position[0];
			swarm->particles[i].position[5] = swarm->particles[i].position[6] = swarm->particles[i].position[7] = swarm->particles[i].position[8] = swarm->particles[i].position[9] = swarm->particles[i].position[4]; 
			swarm->particles[i].position[11] = swarm->particles[i].position[12] = swarm->particles[i].position[10];
			swarm->particles[i].position[14] = swarm->particles[i].position[15] = swarm->particles[i].position[16] = swarm->particles[i].position[13];
			swarm->particles[i].position[18] = swarm->particles[i].position[19] = swarm->particles[i].position[20] = swarm->particles[i].position[21] = swarm->particles[i].position[22] = swarm->particles[i].position[17];
			swarm->particles[i].position[24] = swarm->particles[i].position[25] = swarm->particles[i].position[23];
			swarm->particles[i].position[27] = swarm->particles[i].position[28] = swarm->particles[i].position[29] = swarm->particles[i].position[26];
			swarm->particles[i].position[31] = swarm->particles[i].position[32] = swarm->particles[i].position[33] = swarm->particles[i].position[34] = swarm->particles[i].position[35] = swarm->particles[i].position[30];
			swarm->particles[i].position[37] = swarm->particles[i].position[38] = swarm->particles[i].position[36];
			swarm->particles[i].position[40] = swarm->particles[i].position[41] = swarm->particles[i].position[42] = swarm->particles[i].position[39];
			swarm->particles[i].position[44] = swarm->particles[i].position[45] = swarm->particles[i].position[46] = swarm->particles[i].position[47] = swarm->particles[i].position[48] = swarm->particles[i].position[43];
			swarm->particles[i].position[50] = swarm->particles[i].position[51] = swarm->particles[i].position[49];

			x[1] = x[2] = x[3] = x[0];
			x[5] = x[6] = x[7] = x[8] = x[9] = x[4]; 
			x[11] = x[12] = x[10];
			x[14] = x[15] = x[16] = x[13];
			x[18] = x[19] = x[20] = x[21] = x[22] = x[17];
			x[24] = x[25] = x[23];
			x[27] = x[28] = x[29] = x[26];
			x[31] = x[32] = x[33] = x[34] = x[35] = x[30];
			x[37] = x[38] = x[36];
			x[40] = x[41] = x[42] = x[39];
			x[44] = x[45] = x[46] = x[47] = x[48] = x[43];
			x[50] = x[51] = x[49];

			double tab[] = {11.102791, 14.103566, 19.604806, 25.006202, 30.707597, 39.109767, 44.211008, 56.313953, 60.214884, 76.619070, 78.519535, 99.424651, 100.024806, 113.062946, 126.631473, 145.736434, 156.338760, 162.040155, 180.044651, 199.049302, 213.052868, 238.059070, 262.064961, 263.065271, 288.071473, 293.072713, 309.076589, 313.077674, 338.083876, 347.086047, 355.088062, 362.975349, 384.095194, 387.095969, 388.096279, 418.103721, 422.104651, 444.925271, 459.113798, 480.119070, 497.123256, 512.126977, 574.142326, 722.179070, 797.197674, 853.211628, 930.232248, 1085.269147, 1150.285271, 1350.334884, 1390.344806, 1420.352248, 1550.384496, 1600.396899, 1690.419225, 1880.466357, 1990.493643, 2200.545736, 2290.568062, 2450.607752, 2650.657364, 2800.694574, 3000.744186, 3350.831008};

			int pos;
			for (j = 0; j < N_VAR; j++) {
				pos = round(x[j]);
				x[j] = tab[pos];
			}

			function = 0.;
			for (j = 0; j < N_VAR; j++) {
				function += x[j] * length[j];
			}
			function /= 10;
		#endif

		#ifdef USE_T60C

				// MODIFICADO POR PEDRO
			// Evaluate using Eureka
			// x[0] = 0.5;
			// x[1] = 1;
			// x[2] = 2.5;
			// x[3] = 3;
			// x[4] = 3.5;
			// x[5] = 4;
			// x[6] = 4.5;
			// x[7] = 5;

			// printf("X array: \t");
			// for (j = 0; j < truss->getDimension(); ++j){
			// 	printf("%f ", x[j]);
			// }
			// truss->evaluation(x, valuesArray[i]);
			// // printf("function (objFunc) Erica: %f\t function (objFunc) Eureka: %f\n", function, valuesArray[i][0]);			// function = valuesArray[0]; // Function gets the objective function value
			// // function = valuesArray[i][0]; // Avaliação do Eureka ao invés da Erica
			// printf("\n");
			// // printf("function (objFunc): %f\n",  valuesArray[i][0]);			// function = valuesArray[0]; // Function gets the objective function value
			// printf("Constraints: \t");
			// for (j =0; j < truss->getNumberConstraints() + 2; ++j){
			// 	printf("%f ", valuesArray[i][j]);
			// }
			// printf("\n");
			// exit(3);

			length[0] = 51.763809;
			length[1] = 51.763808;
			length[2] = 51.763809;
			length[3] = 51.763809;
			length[4] = 51.763808;
			length[5] = 51.763809;
			length[6] = 51.763809;
			length[7] = 51.763808;
			length[8] = 51.763809;
			length[9] = 51.763809;
			length[10] = 51.763808;
			length[11] = 51.763809;
			length[12] = 46.587427;
			length[13] = 46.587433;
			length[14] = 46.587427;
			length[15] = 46.587427;
			length[16] = 46.587433;
			length[17] = 46.587427;
			length[18] = 46.587427;
			length[19] = 46.587433;
			length[20] = 46.587427;
			length[21] = 46.587427;
			length[22] = 46.587433;
			length[23] = 46.587427;
			length[24] = 50.115293;
			length[25] = 50.115296;
			length[26] = 50.115294;
			length[27] = 50.115293;
			length[28] = 50.115296;
			length[29] = 50.115294;
			length[30] = 50.115293;
			length[31] = 50.115296;
			length[32] = 50.115294;
			length[33] = 50.115293;
			length[34] = 50.115296;
			length[35] = 50.115294;
			length[36] = 50.115294;
			length[37] = 50.115296;
			length[38] = 50.115293;
			length[39] = 50.115294;
			length[40] = 50.115296;
			length[41] = 50.115293;
			length[42] = 50.115294;
			length[43] = 50.115296;
			length[44] = 50.115293;
			length[45] = 50.115294;
			length[46] = 50.115296;
			length[47] = 50.115293;
			length[48] = 10.000000;
			length[49] = 9.999997;
			length[50] = 9.999997;
			length[51] = 10.000000;
			length[52] = 9.999997;
			length[53] = 9.999997;
			length[54] = 10.000000;
			length[55] = 9.999997;
			length[56] = 9.999997;
			length[57] = 10.000000;
			length[58] = 9.999997;
			length[59] = 9.999997;
	
			int cont, cont1, cont2;
			for (cont = 49; cont < 60; cont++) {
				swarm->particles[i].position[cont] = swarm->particles[i].position[48];
			}
			cont1 = 0;
			for (cont2 = 12; cont2 < 24; cont2++) {
				swarm->particles[i].position[cont2] = swarm->particles[i].position[cont1];
				cont1++;
			}
			cont1 = 24;
			for (cont2 = 36; cont2 < 48; cont2++) {
				swarm->particles[i].position[cont2] = swarm->particles[i].position[cont1];
				cont1++;
			}

			for (cont = 49; cont < 60; cont++) {
				x[cont] = x[48];
			}
			cont1 = 0;
			for (cont2 = 12; cont2 < 24; cont2++) {
				x[cont2] = x[cont1];
				cont1++;
			}
			cont1 = 24;
			for (cont2 = 36; cont2 < 48; cont2++) {
				x[cont2] = x[cont1];
				cont1++;
			}

			function = 0.;
			
			for (j = 0; j < N_VAR; j++) {
				function += x[j] * length[j];
			}
			function /= 10;

			// this->grouping[0] = 1;
			// this->grouping[1] = 2;
			// this->grouping[2] = 3;
			// this->grouping[3] = 4;
			// this->grouping[4] = 5;
			// this->grouping[5] = 6;
			// this->grouping[6] = 7;
			// this->grouping[7] = 8;
			// this->grouping[8] = 9;
			// this->grouping[9] = 10;
			// this->grouping[10] = 11;
			// this->grouping[11] = 12;
			// this->grouping[12] = 1;
			// this->grouping[13] = 2;
			// this->grouping[14] = 3;
			// this->grouping[15] = 4;
			// this->grouping[16] = 5;
			// this->grouping[17] = 6;
			// this->grouping[18] = 7;
			// this->grouping[19] = 8;
			// this->grouping[20] = 9;
			// this->grouping[21] = 10;
			// this->grouping[22] = 11;
			// this->grouping[23] = 12;
			// this->grouping[24] = 13;
			// this->grouping[25] = 14;
			// this->grouping[26] = 15;
			// this->grouping[27] = 16;
			// this->grouping[28] = 17;
			// this->grouping[29] = 18;
			// this->grouping[30] = 19;
			// this->grouping[31] = 20;
			// this->grouping[32] = 21;
			// this->grouping[33] = 22;
			// this->grouping[34] = 23;
			// this->grouping[35] = 24;
			// this->grouping[36] = 13;
			// this->grouping[37] = 14;
			// this->grouping[38] = 15;
			// this->grouping[39] = 16;
			// this->grouping[40] = 17;
			// this->grouping[41] = 18;
			// this->grouping[42] = 19;
			// this->grouping[43] = 20;
			// this->grouping[44] = 21;
			// this->grouping[45] = 22;
			// this->grouping[46] = 23;
			// this->grouping[47] = 24;
			// this->grouping[48] = 0;
			// this->grouping[49] = 0;
			// this->grouping[50] = 0;
			// this->grouping[51] = 0;
			// this->grouping[52] = 0;
			// this->grouping[53] = 0;
			// this->grouping[54] = 0;
			// this->grouping[55] = 0;
			// this->grouping[56] = 0;
			// this->grouping[57] = 0;
			// this->grouping[58] = 0;
			// this->grouping[59] = 0;
			// Best solution (below)
			// x[0] = 2.102395;
			// x[1] = 0.503932;
			// x[2] = 2.000652;
			// x[3] = 1.932965;
			// x[4] = 0.573366;
			// x[5] = 1.893816;
			// x[6] = 1.925602;
			// x[7] = 1.060730;
			// x[8] = 1.764067;
			// x[9] = 1.615935;
			// x[10] = 0.506942;
			// x[11] = 2.129686;
			// x[12] = 2.102395;
			// x[13] = 0.503932;
			// x[14] = 2.000652;
			// x[15] = 1.932965;
			// x[16] = 0.573366;
			// x[17] = 1.893816;
			// x[18] = 1.925602;
			// x[19] = 1.060730;
			// x[20] = 1.764067;
			// x[21] = 1.615935;
			// x[22] = 0.506942;
			// x[23] = 2.129686;
			// x[24] = 1.264224;
			// x[25] = 1.154981;
			// x[26] = 0.508007;
			// x[27] = 0.709542;
			// x[28] = 0.989056;
			// x[29] = 1.140337;
			// x[30] = 1.139156;
			// x[31] = 1.072108;
			// x[32] = 1.077415;
			// x[33] = 0.623711;
			// x[34] = 1.090640;
			// x[35] = 1.259983;
			// x[36] = 1.264224;
			// x[37] = 1.154981;
			// x[38] = 0.508007;
			// x[39] = 0.709542;
			// x[40] = 0.989056;
			// x[41] = 1.140337;
			// x[42] = 1.139156;
			// x[43] = 1.072108;
			// x[44] = 1.077415;
			// x[45] = 0.623711;
			// x[46] = 1.090640;
			// x[47] = 1.259983;
			// x[48] = 1.168084;
			// x[49] = 1.168084;
			// x[50] = 1.168084;
			// x[51] = 1.168084;
			// x[52] = 1.168084;
			// x[53] = 1.168084;
			// x[54] = 1.168084;
			// x[55] = 1.168084;
			// x[56] = 1.168084;
			// x[57] = 1.168084;
			// x[58] = 1.168084;
			// x[59] = 1.168084;

			xTEST[0] = x[48]; // 49-60
			xTEST[1] = x[0]; // 1 e 13
			xTEST[2] = x[1]; // 2 e 14
			xTEST[3] = x[2]; // 3 e 15
			xTEST[4] = x[3]; // 4 e 16
			xTEST[5] = x[4]; // 5 e 17
			xTEST[6] = x[5]; // 6 e 18
			xTEST[7] = x[6]; // 7 e 19
			xTEST[8] = x[7]; // 8 e 20
			xTEST[9] = x[8]; // 9 e 21
			xTEST[10] = x[9]; // 10 e 22
			xTEST[11] = x[10]; // 11 e 23
			xTEST[12] = x[11]; // 12 e 24
			xTEST[13] = x[24]; // 25 e 37
			xTEST[14] = x[25]; // 26 e 38
			xTEST[15] = x[26]; // 27 e 39
			xTEST[16] = x[27]; // 28 e 40
			xTEST[17] = x[28]; // 29 e 41
			xTEST[18] = x[29]; // 30 e 42
			xTEST[19] = x[30]; // 31 e 43
			xTEST[20] = x[31]; // 32 e 44
			xTEST[21] = x[32]; // 33 e 45
			xTEST[22] = x[33]; // 34 e 46
			xTEST[23] = x[34]; // 35 e 47
			xTEST[24] = x[35]; // 36 e 48

			// printf("Grouped project variables for 60bars");
			// int itP;
			// printf("Printing xtests (design variables)");
			// for(itP=0 ;itP<25 ;itP++){
			// 	printf("x[%d] = %f ", itP, xTEST[itP]);
			// }
			// printf("function (objFunc) Erica: %f\t function (objFunc) Eureka: %f\n", function, valuesArray[i][0]);			// function = valuesArray[0]; // Function gets the objective function value
			// exit(5);



			// MODIFICADO POR PEDRO
			// Evaluate using Eureka
			// truss->evaluation(x, valuesArray[i]);
			truss->evaluation(xTEST, valuesArray[i]);
			// printf("generation: %d", generation);
			// int itP;
			// printf("Printing xtests (design variables)");
			// for(itP=0 ;itP<25 ;itP++){
			// 	printf("x[%d] = %f\n", itP, xTEST[itP]);
			// }
			// printf("function (objFunc) Erica: %f\t function (objFunc) Eureka: %f\n", function, valuesArray[i][0]);			// function = valuesArray[0]; // Function gets the objective function value
			// exit(5);
			// if(generation > 230){
			
			// 	printf("generation: %d", generation);
			
			// 	int itP;
			// 	printf("Printing xtests (design variables)");
			// 	for(itP=0 ;itP<25 ;itP++){
			// 		printf("x[%d] = %f ", itP, xTEST[itP]);
			// 	}
			// 	printf("function (objFunc) Erica: %f\t function (objFunc) Eureka: %f\n", function, valuesArray[i][0]);			// function = valuesArray[0]; // Function gets the objective function value
			// 	exit(5);
			// }
			function = valuesArray[i][0]; // Avaliação do Eureka ao invés da Erica
		#endif


		#ifdef USE_T72C
			length[0] = 60.000000;
			length[1] = 60.000000;
			length[2] = 60.000000;
			length[3] = 60.000000;
			length[4] = 134.164079;
			length[5] = 134.164079;
			length[6] = 134.164079;
			length[7] = 134.164079;
			length[8] = 134.164079;
			length[9] = 134.164079;
			length[10] = 134.164079;
			length[11] = 134.164079;
			length[12] = 120.000000;
			length[13] = 120.000000;
			length[14] = 120.000000;
			length[15] = 120.000000;
			length[16] = 169.705627;
			length[17] = 169.705627;
			length[18] = 60.000000;
			length[19] = 60.000000;
			length[20] = 60.000000;
			length[21] = 60.000000;
			length[22] = 134.164079;
			length[23] = 134.164079;
			length[24] = 134.164079;
			length[25] = 134.164079;
			length[26] = 134.164079;
			length[27] = 134.164079;
			length[28] = 134.164079;
			length[29] = 134.164079;
			length[30] = 120.000000;
			length[31] = 120.000000;
			length[32] = 120.000000;
			length[33] = 120.000000;
			length[34] = 169.705627;
			length[35] = 169.705627;
			length[36] = 60.000000;
			length[37] = 60.000000;
			length[38] = 60.000000;
			length[39] = 60.000000;
			length[40] = 134.164079;
			length[41] = 134.164079;
			length[42] = 134.164079;
			length[43] = 134.164079;
			length[44] = 134.164079;
			length[45] = 134.164079;
			length[46] = 134.164079;
			length[47] = 134.164079;
			length[48] = 120.000000;
			length[49] = 120.000000;
			length[50] = 120.000000;
			length[51] = 120.000000;
			length[52] = 169.705627;
			length[53] = 169.705627;
			length[54] = 60.000000;
			length[55] = 60.000000;
			length[56] = 60.000000;
			length[57] = 60.000000;
			length[58] = 134.164079;
			length[59] = 134.164079;
			length[60] = 134.164079;
			length[61] = 134.164079;
			length[62] = 134.164079;
			length[63] = 134.164079;
			length[64] = 134.164079;
			length[65] = 134.164079;
			length[66] = 120.000000;
			length[67] = 120.000000;
			length[68] = 120.000000;
			length[69] = 120.000000;
			length[70] = 169.705627;
			length[71] = 169.705627;

			swarm->particles[i].position[1] = swarm->particles[i].position[2] = swarm->particles[i].position[3] = swarm->particles[i].position[0];
			swarm->particles[i].position[5] = swarm->particles[i].position[6] = swarm->particles[i].position[7] = swarm->particles[i].position[8] = swarm->particles[i].position[9] = swarm->particles[i].position[10] = swarm->particles[i].position[11] = swarm->particles[i].position[4]; 
			swarm->particles[i].position[13] = swarm->particles[i].position[14] = swarm->particles[i].position[15] = swarm->particles[i].position[12];
			swarm->particles[i].position[17] = swarm->particles[i].position[16];
			swarm->particles[i].position[19] = swarm->particles[i].position[20] = swarm->particles[i].position[21] = swarm->particles[i].position[18];
			swarm->particles[i].position[23] = swarm->particles[i].position[24] = swarm->particles[i].position[25] = swarm->particles[i].position[26] = swarm->particles[i].position[27] = swarm->particles[i].position[28] = swarm->particles[i].position[29] = swarm->particles[i].position[22];
			swarm->particles[i].position[31] = swarm->particles[i].position[32] = swarm->particles[i].position[33] = swarm->particles[i].position[30];
			swarm->particles[i].position[35] = swarm->particles[i].position[34];
			swarm->particles[i].position[37] = swarm->particles[i].position[38] = swarm->particles[i].position[39] = swarm->particles[i].position[36];
			swarm->particles[i].position[41] = swarm->particles[i].position[42] = swarm->particles[i].position[43] = swarm->particles[i].position[44] = swarm->particles[i].position[45] = swarm->particles[i].position[46] = swarm->particles[i].position[47] = swarm->particles[i].position[40];
			swarm->particles[i].position[49] = swarm->particles[i].position[50] = swarm->particles[i].position[51] = swarm->particles[i].position[48];
			swarm->particles[i].position[53] = swarm->particles[i].position[52]; 
			swarm->particles[i].position[55] = swarm->particles[i].position[56] = swarm->particles[i].position[57] = swarm->particles[i].position[54];
			swarm->particles[i].position[59] = swarm->particles[i].position[60] = swarm->particles[i].position[61] = swarm->particles[i].position[62] = swarm->particles[i].position[63] = swarm->particles[i].position[64] = swarm->particles[i].position[65] = swarm->particles[i].position[58];
			swarm->particles[i].position[67] = swarm->particles[i].position[68] = swarm->particles[i].position[69] = swarm->particles[i].position[66];
			swarm->particles[i].position[71] = swarm->particles[i].position[70];
			
			x[1] = x[2] = x[3] = x[0];
			x[5] = x[6] = x[7] = x[8] = x[9] = x[10] = x[11] = x[4]; 
			x[13] = x[14] = x[15] = x[12];
			x[17] = x[16];
			x[19] = x[20] = x[21] = x[18];
			x[23] = x[24] = x[25] = x[26] = x[27] = x[28] = x[29] = x[22];
			x[31] = x[32] = x[33] = x[30];
			x[35] = x[34];
			x[37] = x[38] = x[39] = x[36];
			x[41] = x[42] = x[43] = x[44] = x[45] = x[46] = x[47] = x[40];
			x[49] = x[50] = x[51] = x[48];
			x[53] = x[52]; 
			x[55] = x[56] = x[57] = x[54];
			x[59] = x[60] = x[61] = x[62] = x[63] = x[64] = x[65] = x[58];
			x[67] = x[68] = x[69] = x[66];
			x[71] = x[70];

			function = 0.;
			for (j = 0; j < N_VAR; j++) {
				function += x[j] * length[j];
			}
			function /= 10;

			xTEST[0] = x[3];
			xTEST[1] = x[11];
			xTEST[2] = x[15];
			xTEST[3] = x[17];
			xTEST[4] = x[21];
			xTEST[5] = x[29];
			xTEST[6] = x[33];
			xTEST[7] = x[35];
			xTEST[8] = x[39];
			xTEST[9] = x[47];
			xTEST[10] = x[51];
			xTEST[11] = x[53];
			xTEST[12] = x[57];
			xTEST[13] = x[65];
			xTEST[14] = x[69];
			xTEST[15] = x[71];

			// MODIFICADO POR PEDRO
			// Evaluate using Eureka
			truss->evaluation(xTEST, valuesArray[i]);
			// printf("function (objFunc) Erica: %f\t function (objFunc) Eureka: %f\n", function, valuesArray[i][0]);			// function = valuesArray[0]; // Function gets the objective function value
			function = valuesArray[i][0]; // Avaliação do Eureka ao invés da Erica
		#endif


		#ifdef USE_T942C

			length[0] = 237.587878;
			length[1] = 237.587878;
			length[2] = 168.000000;
			length[3] = 168.000000;
			length[4] = 168.000000;
			length[5] = 168.000000;
			length[6] = 168.000000;
			length[7] = 168.000000;
			length[8] = 168.000000;
			length[9] = 168.000000;
			length[10] = 144.000000;
			length[11] = 144.000000;
			length[12] = 144.000000;
			length[13] = 144.000000;
			length[14] = 144.000000;
			length[15] = 144.000000;
			length[16] = 144.000000;
			length[17] = 144.000000;
			length[18] = 221.269067;
			length[19] = 221.269067;
			length[20] = 221.269067;
			length[21] = 221.269067;
			length[22] = 221.269067;
			length[23] = 221.269067;
			length[24] = 221.269067;
			length[25] = 221.269067;
			length[26] = 221.269067;
			length[27] = 221.269067;
			length[28] = 221.269067;
			length[29] = 221.269067;
			length[30] = 221.269067;
			length[31] = 221.269067;
			length[32] = 221.269067;
			length[33] = 221.269067;
			length[34] = 168.000000;
			length[35] = 168.000000;
			length[36] = 168.000000;
			length[37] = 168.000000;
			length[38] = 168.000000;
			length[39] = 168.000000;
			length[40] = 168.000000;
			length[41] = 168.000000;
			length[42] = 168.000000;
			length[43] = 168.000000;
			length[44] = 168.000000;
			length[45] = 168.000000;
			length[46] = 144.000000;
			length[47] = 144.000000;
			length[48] = 144.000000;
			length[49] = 144.000000;
			length[50] = 144.000000;
			length[51] = 144.000000;
			length[52] = 144.000000;
			length[53] = 144.000000;
			length[54] = 144.000000;
			length[55] = 144.000000;
			length[56] = 144.000000;
			length[57] = 144.000000;
			length[58] = 221.269067;
			length[59] = 221.269067;
			length[60] = 221.269067;
			length[61] = 221.269067;
			length[62] = 221.269067;
			length[63] = 221.269067;
			length[64] = 221.269067;
			length[65] = 221.269067;
			length[66] = 221.269067;
			length[67] = 221.269067;
			length[68] = 221.269067;
			length[69] = 221.269067;
			length[70] = 221.269067;
			length[71] = 221.269067;
			length[72] = 221.269067;
			length[73] = 221.269067;
			length[74] = 221.269067;
			length[75] = 221.269067;
			length[76] = 221.269067;
			length[77] = 221.269067;
			length[78] = 221.269067;
			length[79] = 221.269067;
			length[80] = 221.269067;
			length[81] = 221.269067;
			length[82] = 168.000000;
			length[83] = 168.000000;
			length[84] = 168.000000;
			length[85] = 168.000000;
			length[86] = 186.676190;
			length[87] = 186.676190;
			length[88] = 186.676190;
			length[89] = 186.676190;
			length[90] = 302.152279;
			length[91] = 302.152279;
			length[92] = 302.152279;
			length[93] = 302.152279;
			length[94] = 302.152279;
			length[95] = 302.152279;
			length[96] = 302.152279;
			length[97] = 302.152279;
			length[98] = 208.968897;
			length[99] = 208.968897;
			length[100] = 208.968897;
			length[101] = 208.968897;
			length[102] = 208.968897;
			length[103] = 208.968897;
			length[104] = 208.968897;
			length[105] = 208.968897;
			length[106] = 173.170436;
			length[107] = 173.170436;
			length[108] = 173.170436;
			length[109] = 173.170436;
			length[110] = 173.170436;
			length[111] = 173.170436;
			length[112] = 173.170436;
			length[113] = 173.170436;
			length[114] = 173.170436;
			length[115] = 173.170436;
			length[116] = 173.170436;
			length[117] = 173.170436;
			length[118] = 173.170436;
			length[119] = 173.170436;
			length[120] = 173.170436;
			length[121] = 173.170436;
			length[122] = 144.000000;
			length[123] = 144.000000;
			length[124] = 144.000000;
			length[125] = 144.000000;
			length[126] = 144.000000;
			length[127] = 144.000000;
			length[128] = 144.000000;
			length[129] = 144.000000;
			length[130] = 225.219893;
			length[131] = 225.219893;
			length[132] = 225.219893;
			length[133] = 225.219893;
			length[134] = 225.219893;
			length[135] = 225.219893;
			length[136] = 225.219893;
			length[137] = 225.219893;
			length[138] = 225.219893;
			length[139] = 225.219893;
			length[140] = 225.219893;
			length[141] = 225.219893;
			length[142] = 225.219893;
			length[143] = 225.219893;
			length[144] = 225.219893;
			length[145] = 225.219893;
			length[146] = 225.219893;
			length[147] = 225.219893;
			length[148] = 225.219893;
			length[149] = 225.219893;
			length[150] = 225.219893;
			length[151] = 225.219893;
			length[152] = 225.219893;
			length[153] = 225.219893;
			length[154] = 225.219893;
			length[155] = 225.219893;
			length[156] = 225.219893;
			length[157] = 225.219893;
			length[158] = 225.219893;
			length[159] = 225.219893;
			length[160] = 225.219893;
			length[161] = 225.219893;
			length[162] = 144.000000;
			length[163] = 144.000000;
			length[164] = 144.000000;
			length[165] = 144.000000;
			length[166] = 144.000000;
			length[167] = 144.000000;
			length[168] = 144.000000;
			length[169] = 144.000000;
			length[170] = 173.170436;
			length[171] = 173.170436;
			length[172] = 173.170436;
			length[173] = 173.170436;
			length[174] = 173.170436;
			length[175] = 173.170436;
			length[176] = 173.170436;
			length[177] = 173.170436;
			length[178] = 173.170436;
			length[179] = 173.170436;
			length[180] = 173.170436;
			length[181] = 173.170436;
			length[182] = 173.170436;
			length[183] = 173.170436;
			length[184] = 173.170436;
			length[185] = 173.170436;
			length[186] = 144.000000;
			length[187] = 144.000000;
			length[188] = 144.000000;
			length[189] = 144.000000;
			length[190] = 144.000000;
			length[191] = 144.000000;
			length[192] = 144.000000;
			length[193] = 144.000000;
			length[194] = 225.219893;
			length[195] = 225.219893;
			length[196] = 225.219893;
			length[197] = 225.219893;
			length[198] = 225.219893;
			length[199] = 225.219893;
			length[200] = 225.219893;
			length[201] = 225.219893;
			length[202] = 225.219893;
			length[203] = 225.219893;
			length[204] = 225.219893;
			length[205] = 225.219893;
			length[206] = 225.219893;
			length[207] = 225.219893;
			length[208] = 225.219893;
			length[209] = 225.219893;
			length[210] = 225.219893;
			length[211] = 225.219893;
			length[212] = 225.219893;
			length[213] = 225.219893;
			length[214] = 225.219893;
			length[215] = 225.219893;
			length[216] = 225.219893;
			length[217] = 225.219893;
			length[218] = 225.219893;
			length[219] = 225.219893;
			length[220] = 225.219893;
			length[221] = 225.219893;
			length[222] = 225.219893;
			length[223] = 225.219893;
			length[224] = 225.219893;
			length[225] = 225.219893;
			length[226] = 144.000000;
			length[227] = 144.000000;
			length[228] = 144.000000;
			length[229] = 144.000000;
			length[230] = 144.000000;
			length[231] = 144.000000;
			length[232] = 144.000000;
			length[233] = 144.000000;
			length[234] = 173.170436;
			length[235] = 173.170436;
			length[236] = 173.170436;
			length[237] = 173.170436;
			length[238] = 173.170436;
			length[239] = 173.170436;
			length[240] = 173.170436;
			length[241] = 173.170436;
			length[242] = 173.170436;
			length[243] = 173.170436;
			length[244] = 173.170436;
			length[245] = 173.170436;
			length[246] = 173.170436;
			length[247] = 173.170436;
			length[248] = 173.170436;
			length[249] = 173.170436;
			length[250] = 173.170436;
			length[251] = 173.170436;
			length[252] = 173.170436;
			length[253] = 173.170436;
			length[254] = 173.170436;
			length[255] = 173.170436;
			length[256] = 173.170436;
			length[257] = 173.170436;
			length[258] = 144.000000;
			length[259] = 144.000000;
			length[260] = 144.000000;
			length[261] = 144.000000;
			length[262] = 144.000000;
			length[263] = 144.000000;
			length[264] = 144.000000;
			length[265] = 144.000000;
			length[266] = 144.000000;
			length[267] = 144.000000;
			length[268] = 144.000000;
			length[269] = 144.000000;
			length[270] = 225.219893;
			length[271] = 225.219893;
			length[272] = 225.219893;
			length[273] = 225.219893;
			length[274] = 225.219893;
			length[275] = 225.219893;
			length[276] = 225.219893;
			length[277] = 225.219893;
			length[278] = 225.219893;
			length[279] = 225.219893;
			length[280] = 225.219893;
			length[281] = 225.219893;
			length[282] = 225.219893;
			length[283] = 225.219893;
			length[284] = 225.219893;
			length[285] = 225.219893;
			length[286] = 225.219893;
			length[287] = 225.219893;
			length[288] = 225.219893;
			length[289] = 225.219893;
			length[290] = 225.219893;
			length[291] = 225.219893;
			length[292] = 225.219893;
			length[293] = 225.219893;
			length[294] = 225.219893;
			length[295] = 225.219893;
			length[296] = 225.219893;
			length[297] = 225.219893;
			length[298] = 225.219893;
			length[299] = 225.219893;
			length[300] = 225.219893;
			length[301] = 225.219893;
			length[302] = 225.219893;
			length[303] = 225.219893;
			length[304] = 225.219893;
			length[305] = 225.219893;
			length[306] = 225.219893;
			length[307] = 225.219893;
			length[308] = 225.219893;
			length[309] = 225.219893;
			length[310] = 225.219893;
			length[311] = 225.219893;
			length[312] = 225.219893;
			length[313] = 225.219893;
			length[314] = 225.219893;
			length[315] = 225.219893;
			length[316] = 225.219893;
			length[317] = 225.219893;
			length[318] = 144.000000;
			length[319] = 144.000000;
			length[320] = 144.000000;
			length[321] = 144.000000;
			length[322] = 144.000000;
			length[323] = 144.000000;
			length[324] = 144.000000;
			length[325] = 144.000000;
			length[326] = 144.000000;
			length[327] = 144.000000;
			length[328] = 144.000000;
			length[329] = 144.000000;
			length[330] = 173.170436;
			length[331] = 173.170436;
			length[332] = 173.170436;
			length[333] = 173.170436;
			length[334] = 173.170436;
			length[335] = 173.170436;
			length[336] = 173.170436;
			length[337] = 173.170436;
			length[338] = 186.676190;
			length[339] = 186.676190;
			length[340] = 186.676190;
			length[341] = 186.676190;
			length[342] = 293.264386;
			length[343] = 293.264386;
			length[344] = 293.264386;
			length[345] = 293.264386;
			length[346] = 293.264386;
			length[347] = 293.264386;
			length[348] = 293.264386;
			length[349] = 293.264386;
			length[350] = 208.968897;
			length[351] = 208.968897;
			length[352] = 208.968897;
			length[353] = 208.968897;
			length[354] = 208.968897;
			length[355] = 208.968897;
			length[356] = 208.968897;
			length[357] = 208.968897;
			length[358] = 186.676190;
			length[359] = 186.676190;
			length[360] = 186.676190;
			length[361] = 186.676190;
			length[362] = 186.676190;
			length[363] = 186.676190;
			length[364] = 186.676190;
			length[365] = 186.676190;
			length[366] = 173.170436;
			length[367] = 173.170436;
			length[368] = 173.170436;
			length[369] = 173.170436;
			length[370] = 173.170436;
			length[371] = 173.170436;
			length[372] = 173.170436;
			length[373] = 173.170436;
			length[374] = 173.170436;
			length[375] = 173.170436;
			length[376] = 173.170436;
			length[377] = 173.170436;
			length[378] = 173.170436;
			length[379] = 173.170436;
			length[380] = 173.170436;
			length[381] = 173.170436;
			length[382] = 168.000000;
			length[383] = 168.000000;
			length[384] = 168.000000;
			length[385] = 168.000000;
			length[386] = 168.000000;
			length[387] = 168.000000;
			length[388] = 168.000000;
			length[389] = 168.000000;
			length[390] = 144.000000;
			length[391] = 144.000000;
			length[392] = 144.000000;
			length[393] = 144.000000;
			length[394] = 144.000000;
			length[395] = 144.000000;
			length[396] = 144.000000;
			length[397] = 144.000000;
			length[398] = 225.219893;
			length[399] = 225.219893;
			length[400] = 225.219893;
			length[401] = 225.219893;
			length[402] = 225.219893;
			length[403] = 225.219893;
			length[404] = 225.219893;
			length[405] = 225.219893;
			length[406] = 225.219893;
			length[407] = 225.219893;
			length[408] = 225.219893;
			length[409] = 225.219893;
			length[410] = 225.219893;
			length[411] = 225.219893;
			length[412] = 225.219893;
			length[413] = 225.219893;
			length[414] = 225.219893;
			length[415] = 225.219893;
			length[416] = 225.219893;
			length[417] = 225.219893;
			length[418] = 225.219893;
			length[419] = 225.219893;
			length[420] = 144.000000;
			length[421] = 144.000000;
			length[422] = 225.219893;
			length[423] = 225.219893;
			length[424] = 225.219893;
			length[425] = 225.219893;
			length[426] = 225.219893;
			length[427] = 225.219893;
			length[428] = 225.219893;
			length[429] = 225.219893;
			length[430] = 144.000000;
			length[431] = 144.000000;
			length[432] = 144.000000;
			length[433] = 144.000000;
			length[434] = 144.000000;
			length[435] = 144.000000;
			length[436] = 144.000000;
			length[437] = 144.000000;
			length[438] = 144.000000;
			length[439] = 144.000000;
			length[440] = 144.000000;
			length[441] = 144.000000;
			length[442] = 144.000000;
			length[443] = 144.000000;
			length[444] = 144.000000;
			length[445] = 144.000000;
			length[446] = 221.269067;
			length[447] = 221.269067;
			length[448] = 221.269067;
			length[449] = 221.269067;
			length[450] = 221.269067;
			length[451] = 221.269067;
			length[452] = 221.269067;
			length[453] = 221.269067;
			length[454] = 221.269067;
			length[455] = 221.269067;
			length[456] = 221.269067;
			length[457] = 221.269067;
			length[458] = 221.269067;
			length[459] = 221.269067;
			length[460] = 221.269067;
			length[461] = 221.269067;
			length[462] = 173.170436;
			length[463] = 173.170436;
			length[464] = 173.170436;
			length[465] = 173.170436;
			length[466] = 173.170436;
			length[467] = 173.170436;
			length[468] = 173.170436;
			length[469] = 173.170436;
			length[470] = 173.170436;
			length[471] = 173.170436;
			length[472] = 173.170436;
			length[473] = 173.170436;
			length[474] = 173.170436;
			length[475] = 173.170436;
			length[476] = 173.170436;
			length[477] = 173.170436;
			length[478] = 173.170436;
			length[479] = 173.170436;
			length[480] = 173.170436;
			length[481] = 173.170436;
			length[482] = 173.170436;
			length[483] = 173.170436;
			length[484] = 173.170436;
			length[485] = 173.170436;
			length[486] = 168.000000;
			length[487] = 168.000000;
			length[488] = 168.000000;
			length[489] = 168.000000;
			length[490] = 168.000000;
			length[491] = 168.000000;
			length[492] = 168.000000;
			length[493] = 168.000000;
			length[494] = 168.000000;
			length[495] = 168.000000;
			length[496] = 168.000000;
			length[497] = 168.000000;
			length[498] = 144.000000;
			length[499] = 144.000000;
			length[500] = 144.000000;
			length[501] = 144.000000;
			length[502] = 144.000000;
			length[503] = 144.000000;
			length[504] = 144.000000;
			length[505] = 144.000000;
			length[506] = 144.000000;
			length[507] = 144.000000;
			length[508] = 144.000000;
			length[509] = 144.000000;
			length[510] = 225.219893;
			length[511] = 225.219893;
			length[512] = 225.219893;
			length[513] = 225.219893;
			length[514] = 225.219893;
			length[515] = 225.219893;
			length[516] = 225.219893;
			length[517] = 225.219893;
			length[518] = 225.219893;
			length[519] = 225.219893;
			length[520] = 225.219893;
			length[521] = 225.219893;
			length[522] = 225.219893;
			length[523] = 225.219893;
			length[524] = 225.219893;
			length[525] = 225.219893;
			length[526] = 225.219893;
			length[527] = 225.219893;
			length[528] = 225.219893;
			length[529] = 225.219893;
			length[530] = 225.219893;
			length[531] = 225.219893;
			length[532] = 144.000000;
			length[533] = 144.000000;
			length[534] = 225.219893;
			length[535] = 225.219893;
			length[536] = 225.219893;
			length[537] = 225.219893;
			length[538] = 225.219893;
			length[539] = 225.219893;
			length[540] = 225.219893;
			length[541] = 225.219893;
			length[542] = 225.219893;
			length[543] = 225.219893;
			length[544] = 225.219893;
			length[545] = 225.219893;
			length[546] = 225.219893;
			length[547] = 225.219893;
			length[548] = 225.219893;
			length[549] = 225.219893;
			length[550] = 225.219893;
			length[551] = 225.219893;
			length[552] = 225.219893;
			length[553] = 225.219893;
			length[554] = 225.219893;
			length[555] = 225.219893;
			length[556] = 225.219893;
			length[557] = 225.219893;
			length[558] = 144.000000;
			length[559] = 144.000000;
			length[560] = 144.000000;
			length[561] = 144.000000;
			length[562] = 144.000000;
			length[563] = 144.000000;
			length[564] = 144.000000;
			length[565] = 144.000000;
			length[566] = 144.000000;
			length[567] = 144.000000;
			length[568] = 144.000000;
			length[569] = 144.000000;
			length[570] = 144.000000;
			length[571] = 144.000000;
			length[572] = 144.000000;
			length[573] = 144.000000;
			length[574] = 144.000000;
			length[575] = 144.000000;
			length[576] = 144.000000;
			length[577] = 144.000000;
			length[578] = 144.000000;
			length[579] = 144.000000;
			length[580] = 144.000000;
			length[581] = 144.000000;
			length[582] = 221.269067;
			length[583] = 221.269067;
			length[584] = 221.269067;
			length[585] = 221.269067;
			length[586] = 225.219893;
			length[587] = 225.219893;
			length[588] = 225.219893;
			length[589] = 221.269067;
			length[590] = 221.269067;
			length[591] = 221.269067;
			length[592] = 221.269067;
			length[593] = 221.269067;
			length[594] = 221.269067;
			length[595] = 221.269067;
			length[596] = 221.269067;
			length[597] = 221.269067;
			length[598] = 221.269067;
			length[599] = 221.269067;
			length[600] = 221.269067;
			length[601] = 221.269067;
			length[602] = 221.269067;
			length[603] = 221.269067;
			length[604] = 221.269067;
			length[605] = 221.269067;
			length[606] = 173.170436;
			length[607] = 173.170436;
			length[608] = 173.170436;
			length[609] = 173.170436;
			length[610] = 173.170436;
			length[611] = 173.170436;
			length[612] = 173.170436;
			length[613] = 173.170436;
			length[614] = 173.170436;
			length[615] = 173.170436;
			length[616] = 173.170436;
			length[617] = 173.170436;
			length[618] = 173.170436;
			length[619] = 173.170436;
			length[620] = 173.170436;
			length[621] = 173.170436;
			length[622] = 173.170436;
			length[623] = 173.170436;
			length[624] = 173.170436;
			length[625] = 173.170436;
			length[626] = 173.170436;
			length[627] = 173.170436;
			length[628] = 173.170436;
			length[629] = 173.170436;
			length[630] = 168.000000;
			length[631] = 168.000000;
			length[632] = 168.000000;
			length[633] = 168.000000;
			length[634] = 168.000000;
			length[635] = 168.000000;
			length[636] = 168.000000;
			length[637] = 168.000000;
			length[638] = 168.000000;
			length[639] = 168.000000;
			length[640] = 168.000000;
			length[641] = 168.000000;
			length[642] = 144.000000;
			length[643] = 144.000000;
			length[644] = 144.000000;
			length[645] = 144.000000;
			length[646] = 144.000000;
			length[647] = 144.000000;
			length[648] = 144.000000;
			length[649] = 144.000000;
			length[650] = 144.000000;
			length[651] = 144.000000;
			length[652] = 144.000000;
			length[653] = 144.000000;
			length[654] = 225.219893;
			length[655] = 225.219893;
			length[656] = 225.219893;
			length[657] = 225.219893;
			length[658] = 225.219893;
			length[659] = 225.219893;
			length[660] = 225.219893;
			length[661] = 225.219893;
			length[662] = 225.219893;
			length[663] = 225.219893;
			length[664] = 225.219893;
			length[665] = 225.219893;
			length[666] = 225.219893;
			length[667] = 225.219893;
			length[668] = 225.219893;
			length[669] = 225.219893;
			length[670] = 225.219893;
			length[671] = 225.219893;
			length[672] = 225.219893;
			length[673] = 225.219893;
			length[674] = 225.219893;
			length[675] = 225.219893;
			length[676] = 144.000000;
			length[677] = 144.000000;
			length[678] = 225.219893;
			length[679] = 225.219893;
			length[680] = 225.219893;
			length[681] = 225.219893;
			length[682] = 225.219893;
			length[683] = 225.219893;
			length[684] = 225.219893;
			length[685] = 225.219893;
			length[686] = 225.219893;
			length[687] = 225.219893;
			length[688] = 225.219893;
			length[689] = 225.219893;
			length[690] = 225.219893;
			length[691] = 225.219893;
			length[692] = 225.219893;
			length[693] = 225.219893;
			length[694] = 225.219893;
			length[695] = 225.219893;
			length[696] = 225.219893;
			length[697] = 225.219893;
			length[698] = 225.219893;
			length[699] = 225.219893;
			length[700] = 225.219893;
			length[701] = 225.219893;
			length[702] = 144.000000;
			length[703] = 144.000000;
			length[704] = 144.000000;
			length[705] = 144.000000;
			length[706] = 144.000000;
			length[707] = 144.000000;
			length[708] = 144.000000;
			length[709] = 144.000000;
			length[710] = 144.000000;
			length[711] = 144.000000;
			length[712] = 144.000000;
			length[713] = 144.000000;
			length[714] = 144.000000;
			length[715] = 144.000000;
			length[716] = 144.000000;
			length[717] = 144.000000;
			length[718] = 144.000000;
			length[719] = 144.000000;
			length[720] = 144.000000;
			length[721] = 144.000000;
			length[722] = 144.000000;
			length[723] = 144.000000;
			length[724] = 144.000000;
			length[725] = 144.000000;
			length[726] = 221.269067;
			length[727] = 221.269067;
			length[728] = 221.269067;
			length[729] = 221.269067;
			length[730] = 225.219893;
			length[731] = 225.219893;
			length[732] = 225.219893;
			length[733] = 221.269067;
			length[734] = 221.269067;
			length[735] = 221.269067;
			length[736] = 221.269067;
			length[737] = 221.269067;
			length[738] = 221.269067;
			length[739] = 221.269067;
			length[740] = 221.269067;
			length[741] = 221.269067;
			length[742] = 221.269067;
			length[743] = 221.269067;
			length[744] = 221.269067;
			length[745] = 221.269067;
			length[746] = 221.269067;
			length[747] = 221.269067;
			length[748] = 221.269067;
			length[749] = 221.269067;
			length[750] = 173.170436;
			length[751] = 173.170436;
			length[752] = 173.170436;
			length[753] = 173.170436;
			length[754] = 173.170436;
			length[755] = 173.170436;
			length[756] = 173.170436;
			length[757] = 173.170436;
			length[758] = 173.170436;
			length[759] = 173.170436;
			length[760] = 173.170436;
			length[761] = 173.170436;
			length[762] = 173.170436;
			length[763] = 173.170436;
			length[764] = 173.170436;
			length[765] = 173.170436;
			length[766] = 173.170436;
			length[767] = 173.170436;
			length[768] = 173.170436;
			length[769] = 173.170436;
			length[770] = 173.170436;
			length[771] = 173.170436;
			length[772] = 173.170436;
			length[773] = 173.170436;
			length[774] = 168.000000;
			length[775] = 168.000000;
			length[776] = 168.000000;
			length[777] = 168.000000;
			length[778] = 168.000000;
			length[779] = 168.000000;
			length[780] = 168.000000;
			length[781] = 168.000000;
			length[782] = 168.000000;
			length[783] = 168.000000;
			length[784] = 168.000000;
			length[785] = 168.000000;
			length[786] = 144.000000;
			length[787] = 144.000000;
			length[788] = 144.000000;
			length[789] = 144.000000;
			length[790] = 144.000000;
			length[791] = 144.000000;
			length[792] = 144.000000;
			length[793] = 144.000000;
			length[794] = 144.000000;
			length[795] = 144.000000;
			length[796] = 144.000000;
			length[797] = 144.000000;
			length[798] = 225.219893;
			length[799] = 225.219893;
			length[800] = 225.219893;
			length[801] = 225.219893;
			length[802] = 225.219893;
			length[803] = 225.219893;
			length[804] = 225.219893;
			length[805] = 225.219893;
			length[806] = 225.219893;
			length[807] = 225.219893;
			length[808] = 225.219893;
			length[809] = 225.219893;
			length[810] = 225.219893;
			length[811] = 225.219893;
			length[812] = 225.219893;
			length[813] = 225.219893;
			length[814] = 225.219893;
			length[815] = 225.219893;
			length[816] = 225.219893;
			length[817] = 225.219893;
			length[818] = 225.219893;
			length[819] = 225.219893;
			length[820] = 144.000000;
			length[821] = 144.000000;
			length[822] = 225.219893;
			length[823] = 225.219893;
			length[824] = 225.219893;
			length[825] = 225.219893;
			length[826] = 225.219893;
			length[827] = 225.219893;
			length[828] = 225.219893;
			length[829] = 225.219893;
			length[830] = 225.219893;
			length[831] = 225.219893;
			length[832] = 225.219893;
			length[833] = 225.219893;
			length[834] = 225.219893;
			length[835] = 225.219893;
			length[836] = 225.219893;
			length[837] = 225.219893;
			length[838] = 225.219893;
			length[839] = 225.219893;
			length[840] = 225.219893;
			length[841] = 225.219893;
			length[842] = 225.219893;
			length[843] = 225.219893;
			length[844] = 225.219893;
			length[845] = 225.219893;
			length[846] = 144.000000;
			length[847] = 144.000000;
			length[848] = 144.000000;
			length[849] = 144.000000;
			length[850] = 144.000000;
			length[851] = 144.000000;
			length[852] = 144.000000;
			length[853] = 144.000000;
			length[854] = 144.000000;
			length[855] = 144.000000;
			length[856] = 144.000000;
			length[857] = 144.000000;
			length[858] = 144.000000;
			length[859] = 144.000000;
			length[860] = 144.000000;
			length[861] = 144.000000;
			length[862] = 144.000000;
			length[863] = 144.000000;
			length[864] = 144.000000;
			length[865] = 144.000000;
			length[866] = 144.000000;
			length[867] = 144.000000;
			length[868] = 144.000000;
			length[869] = 144.000000;
			length[870] = 221.269067;
			length[871] = 221.269067;
			length[872] = 221.269067;
			length[873] = 221.269067;
			length[874] = 225.219893;
			length[875] = 225.219893;
			length[876] = 225.219893;
			length[877] = 221.269067;
			length[878] = 221.269067;
			length[879] = 221.269067;
			length[880] = 221.269067;
			length[881] = 221.269067;
			length[882] = 221.269067;
			length[883] = 221.269067;
			length[884] = 221.269067;
			length[885] = 221.269067;
			length[886] = 221.269067;
			length[887] = 221.269067;
			length[888] = 221.269067;
			length[889] = 221.269067;
			length[890] = 221.269067;
			length[891] = 221.269067;
			length[892] = 221.269067;
			length[893] = 221.269067;
			length[894] = 173.170436;
			length[895] = 173.170436;
			length[896] = 173.170436;
			length[897] = 173.170436;
			length[898] = 173.170436;
			length[899] = 173.170436;
			length[900] = 173.170436;
			length[901] = 173.170436;
			length[902] = 168.000000;
			length[903] = 168.000000;
			length[904] = 168.000000;
			length[905] = 168.000000;
			length[906] = 186.676190;
			length[907] = 186.676190;
			length[908] = 186.676190;
			length[909] = 186.676190;
			length[910] = 293.264386;
			length[911] = 293.264386;
			length[912] = 293.264386;
			length[913] = 293.264386;
			length[914] = 293.264386;
			length[915] = 293.264386;
			length[916] = 293.264386;
			length[917] = 293.264386;
			length[918] = 254.629142;
			length[919] = 254.629142;
			length[920] = 254.629142;
			length[921] = 254.629142;
			length[922] = 254.629142;
			length[923] = 254.629142;
			length[924] = 254.629142;
			length[925] = 254.629142;
			length[926] = 166.709328;
			length[927] = 166.709328;
			length[928] = 166.709328;
			length[929] = 166.709328;
			length[930] = 166.709328;
			length[931] = 166.709328;
			length[932] = 166.709328;
			length[933] = 166.709328;
			length[934] = 236.676995;
			length[935] = 236.676995;
			length[936] = 236.676995;
			length[937] = 236.676995;
			length[938] = 236.676995;
			length[939] = 236.676995;
			length[940] = 236.676995;
			length[941] = 236.676995;

			/*swarm->particles[i].position[1] = swarm->particles[i].position[2] = swarm->particles[i].position[3] = swarm->particles[i].position[0];
			swarm->particles[i].position[5] = swarm->particles[i].position[6] = swarm->particles[i].position[7] = swarm->particles[i].position[8] = swarm->particles[i].position[9] = swarm->particles[i].position[10] = swarm->particles[i].position[11] = swarm->particles[i].position[4]; 
			swarm->particles[i].position[13] = swarm->particles[i].position[14] = swarm->particles[i].position[15] = swarm->particles[i].position[12];
			swarm->particles[i].position[17] = swarm->particles[i].position[16];
			swarm->particles[i].position[19] = swarm->particles[i].position[20] = swarm->particles[i].position[21] = swarm->particles[i].position[18];
			swarm->particles[i].position[23] = swarm->particles[i].position[24] = swarm->particles[i].position[25] = swarm->particles[i].position[26] = swarm->particles[i].position[27] = swarm->particles[i].position[28] = swarm->particles[i].position[29] = swarm->particles[i].position[22];
			swarm->particles[i].position[31] = swarm->particles[i].position[32] = swarm->particles[i].position[33] = swarm->particles[i].position[30];
			swarm->particles[i].position[35] = swarm->particles[i].position[34];
			swarm->particles[i].position[37] = swarm->particles[i].position[38] = swarm->particles[i].position[39] = swarm->particles[i].position[36];
			swarm->particles[i].position[41] = swarm->particles[i].position[42] = swarm->particles[i].position[43] = swarm->particles[i].position[44] = swarm->particles[i].position[45] = swarm->particles[i].position[46] = swarm->particles[i].position[47] = swarm->particles[i].position[40];
			swarm->particles[i].position[49] = swarm->particles[i].position[50] = swarm->particles[i].position[51] = swarm->particles[i].position[48];
			swarm->particles[i].position[53] = swarm->particles[i].position[52]; 
			swarm->particles[i].position[55] = swarm->particles[i].position[56] = swarm->particles[i].position[57] = swarm->particles[i].position[54];
			swarm->particles[i].position[59] = swarm->particles[i].position[60] = swarm->particles[i].position[61] = swarm->particles[i].position[62] = swarm->particles[i].position[63] = swarm->particles[i].position[64] = swarm->particles[i].position[65] = swarm->particles[i].position[58];
			swarm->particles[i].position[67] = swarm->particles[i].position[68] = swarm->particles[i].position[69] = swarm->particles[i].position[66];
			swarm->particles[i].position[71] = swarm->particles[i].position[70];
			
			x[1] = x[2] = x[3] = x[0];
			x[5] = x[6] = x[7] = x[8] = x[9] = x[10] = x[11] = x[4]; 
			x[13] = x[14] = x[15] = x[12];
			x[17] = x[16];
			x[19] = x[20] = x[21] = x[18];
			x[23] = x[24] = x[25] = x[26] = x[27] = x[28] = x[29] = x[22];
			x[31] = x[32] = x[33] = x[30];
			x[35] = x[34];
			x[37] = x[38] = x[39] = x[36];
			x[41] = x[42] = x[43] = x[44] = x[45] = x[46] = x[47] = x[40];
			x[49] = x[50] = x[51] = x[48];
			x[53] = x[52]; 
			x[55] = x[56] = x[57] = x[54];
			x[59] = x[60] = x[61] = x[62] = x[63] = x[64] = x[65] = x[58];
			x[67] = x[68] = x[69] = x[66];
			x[71] = x[70];*/

			function = 0.;
			for (j = 0; j < N_VAR; j++) {
				function += x[j] * length[j];
			}
			function /= 10;
		#endif


		/*for (j = 0; j < N_VAR; j++)
			printf("x[%d] = %lf\t", j, x[j]);
		printf("function = %lf\n", function);*/

		swarm->particles[i].fitness = function;
		functionAvg += function;
	}
	functionAvg /= MAX_PAR;


	// CHANGE VARIANTS
	// printf("pré apm (dentro de fitness)\n");
	// APM (swarm, functionAvg, k, violationAcum, valuesArray);
	med_3_APM (swarm, functionAvg, k, violationAcum, valuesArray);
	// printf("pos APM!");
	// exit(1);
	//sporadic_APM (swarm, functionAvg, k, violationAcum, generation);
	//sporadic_acumulation_APM (swarm, functionAvg, k, violationAcum, generation);
	//monotonic_APM (swarm, functionAvg, k_aux, violationAcum);
	//monotonic_sporadic_APM (swarm, functionAvg, k, violationAcum, generation);
	//damping_APM (swarm, functionAvg, k, violationAcum);
	//worst_APM (swarm, functionAvg, k, violationAcum);
	//worst_2_APM (swarm, functionAvg, k, violationAcum); // Verificar essa função, erro no calculo do 'k' utilizando 'worst'. Com 'functionAvgK' funciona.
	//worst_3_APM (swarm, functionAvg, k, violationAcum); // Mesmo problema do 'worst_2_APM'
	//med_APM (swarm, functionAvg, k, violationAcum);
	//med_2_APM (swarm, functionAvg, k, violationAcum);
	//med_4_APM (swarm, functionAvg, k, violationAcum);
	//med_5_APM (swarm, functionAvg, k, violationAcum);swarm->particles[i].position
	//med_6_APM (swarm, functionAvg, k, violationAcum);
	//med_7_APM (swarm, functionAvg, k, violationAcum);
	//med_worst_APM (swarm, functionAvg, k, violationAcum);
	//med_worst_2_APM (swarm, functionAvg, k, violationAcum); //Mesmo problema do 'worst_2_APM'
}

/**
 * Calc_best function - calculate pBest and gBest.
 */
void calc_best (swarm_t * swarm) {
	int i, j, m;	
	double sumViolation1, sumViolation2;

	for (j = 0; j < MAX_PAR; j++) {
		sumViolation1 = sumViolation2 = 0.;
		for (i = 0; i < N_CON; i++) {
			sumViolation1 += swarm->particles[j].v[i];
			sumViolation2 += swarm->particles[j].pBest.v[i];
		}
		//printf("sumpBest - %lf %lf\n", sumViolation1, sumViolation2);
		//printf("\n%lf\n%lf\n%lf\n", swarm->particles[i].fitnessAPM, swarm->particles[i].pBest.fitnessAPM, swarm->gBest.fitnessAPM);
		if (((sumViolation1 == 0.) && (sumViolation2 == 0.)) || ((sumViolation1 != 0.) && (sumViolation2 != 0.))) {
			if (swarm->particles[j].fitnessAPM < swarm->particles[j].pBest.fitnessAPM) {
				//printf("entroupBest %d\n", j);
				swarm->particles[j].pBest.fitness = swarm->particles[j].fitness;
				swarm->particles[j].pBest.fitnessAPM = swarm->particles[j].fitnessAPM;
				for (m = 0; m < N_VAR; m++) {
					swarm->particles[j].pBest.position[m] = swarm->particles[j].position[m];
				}
				for (m = 0; m < N_CON; m++) {
					swarm->particles[j].pBest.v[m] = swarm->particles[j].v[m];
				}
			}
		}
		else if (sumViolation1 == 0.) {
			swarm->particles[j].pBest.fitness = swarm->particles[j].fitness;
			swarm->particles[j].pBest.fitnessAPM = swarm->particles[j].fitnessAPM;
			for (m = 0; m < N_VAR; m++) {
				swarm->particles[j].pBest.position[m] = swarm->particles[j].position[m];
			}
			for (m = 0; m < N_CON; m++) {
				swarm->particles[j].pBest.v[m] = swarm->particles[j].v[m];
			}
		}
		
			
		//printf("pBest: %.12lf\n", swarm->particles[j].pBest.fitness);
		sumViolation1 = sumViolation2 = 0.;
		for (i = 0; i < N_CON; i++) {
			sumViolation1 += swarm->particles[j].pBest.v[i];
			sumViolation2 += swarm->gBest.v[i];
		}
		if (((sumViolation1 == 0.) && (sumViolation2 == 0.)) || ((sumViolation1 != 0.) && (sumViolation2 != 0.))) {
			if (swarm->particles[j].pBest.fitnessAPM < swarm->gBest.fitnessAPM) {
				swarm->gBest.fitness = swarm->particles[j].pBest.fitness;
				swarm->gBest.fitnessAPM = swarm->particles[j].pBest.fitnessAPM;
				for (m = 0; m < N_VAR; m++) {
					swarm->gBest.position[m] = swarm->particles[j].pBest.position[m];
				}
				for (m = 0; m < N_CON; m++) {
					swarm->gBest.v[m] = swarm->particles[j].pBest.v[m];
				}
			}
		}
		else if (sumViolation1 == 0.) {
			swarm->gBest.fitness = swarm->particles[j].pBest.fitness;
			swarm->gBest.fitnessAPM = swarm->particles[j].pBest.fitnessAPM;
			for (m = 0; m < N_VAR; m++) {
				swarm->gBest.position[m] = swarm->particles[j].pBest.position[m];
			}
			for (m = 0; m < N_CON; m++) {
				swarm->gBest.v[m] = swarm->particles[j].pBest.v[m];
			}
		}
	}
	//printf(" \ngBest: %.12lf\n\n", swarm->gBest.fitness);
}

/**
 * Update function - update velocity and position.
 */
void update (swarm_t * swarm, double lowerBound[], double upperBound[], double alfa) {
	int j, m, c1, c2;
	double rand1, rand2, rand3, rand4, rand5, vCrazy, p, newVelocity, newPosition;
	c1 = c2 = 2.05;
	/*rand5 = rand () % 10000;
	rand5 /= 10000;
	vCrazy = rand5;*/
	vCrazy = 0.0001;
	
	for (j = 0; j < MAX_PAR; j++) { 	
		
		for (m = 0; m < N_VAR; m++) {
			// velocity
			// rands
			rand1 = rand () % 10000;
			rand1 /= 10000;
			rand2 = rand () % 10000;
			rand2 /= 10000;
			rand3 = rand () % 10000;
			rand3 /= 10000;
			if (rand3 <= 0.05)
				rand3 = -1;
			else
				rand3 = 1;
			rand4 = rand () % 10000;
			rand4 /= 10000;
			if (rand4 <= 0.5)
				p = 1;
			else
				p = 0;
			if (rand4 >= 0.5)
				rand4 = -1;
			else 
				rand4 = 1;	
		
			//printf("\nc1 - %d c2 - %d rand1 - %lf rand2 - %lf rand3 - %lf rand4 - %lf alfa - %lf\n", c1, c2, rand1, rand2, rand3, rand4, alfa);
			
			newVelocity = rand2 * rand3 * swarm->particles[j].velocity[m] + (1 - rand2) * c1 * rand1 * (swarm->particles[j].pBest.position[m] - swarm->particles[j].position[m]) + (1 - rand2) * c2 * (1 - rand1) * (swarm->gBest.position[m] - swarm->particles[j].position[m]);
			newVelocity += p * rand4 * vCrazy;
			//printf("newVelocity = %lf\n", newVelocity);
			swarm->particles[j].velocity[m] = newVelocity;

			// position
			newPosition = swarm->particles[j].position[m] + swarm->particles[j].velocity[m];
			if (newPosition > upperBound[m])
				newPosition = upperBound[m];
			else if (newPosition < lowerBound[m])
				newPosition = lowerBound[m];
			//printf("newPosition = %.12lf\n", newPosition);
			swarm->particles[j].position[m] = newPosition;
		}
	}
}

double valueToReach(){
  double fOptima = NULL;

	#ifdef USE_STRING
		fOptima = 0.012665 + 1e-8;
	#endif

	#ifdef USE_REDUCER
		fOptima = 2994.471066 + 1e-7;
	#endif

	#ifdef USE_WELDED_MOD1_7
		fOptima = 1.7248523110932348 + 1e-8;
	#endif

	#ifdef USE_PRESSURE
		fOptima = 6059.701660 + 1e-6;
	#endif

	return fOptima;
}

/**
 * Deallocate function - deallocates the swarm used.
 */
void deallocate (swarm_t * swarm) {
	int i;

	for (i = 0; i < MAX_PAR; i++) {
		free(swarm->particles[i].velocity);
		free(swarm->particles[i].position);
		free(swarm->particles[i].v);
		free(swarm->particles[i].pBest.position);
		free(swarm->particles[i].pBest.v);
	}
	free (swarm->gBest.position);
	free (swarm->gBest.v);
	free (swarm);
}

void deallocateBest (best_particle_t * best) {
	int i;

	free(best->v);
	free(best->position);
	free(best);
	// free(best->v);

	
}

int main () {
	// srand(time(NULL));
	srand(12);	
	// printf("truss.getDimensio(): %d\n", trussTestObj.getDimension());
	// int arr [5] ={1,2,3,4,5};
	
	// printf("antes de declarar truss");
	F101Truss10Bar *truss= new F101Truss10Bar();
	// F103Truss25Bar *truss= new F103Truss25Bar();
	// F105Truss60Bar *truss= new F105Truss60Bar();
	// F105Truss60Bar *truss= new F107Truss72Bar();
	// printf("truss.getDimensio(): %d\n", trussTestObj1->getDimension());
	
	// printf("dimension: %d\n", truss->getDimension());
	// printf("constraints: %d\n", truss->getNumberConstraints());
	// exit(3);

	


	// // Bloco abaixo funciona
	// testNamespace::testClass testClassObj;
	// printf("bar: %d", testClassObj.bar(10));

	


	// // Bloco abaixo funciona
	// printf("foo: %d", foo(10));
	// return 0;

	int i, j, m, run, index, total = 0;
	double lowerBound[N_VAR], upperBound[N_VAR], vectorAux[MAX_RUN], alfa, ** violationAcum, * k, feasibleAverage , feasibleBest, feasibleWorst, feasibleMedian, feasiblePosition[N_VAR], bestFeasiblePosition[N_VAR], dp;
	FILE * output = fopen ("outs/t10c_apm.txt", "a+");
	if (output == NULL) {
		printf("\nError..\n");
		return 0;
	}

	FILE * outputPP = fopen ("outs/t10c_apm_pp.txt", "a+");
	if (outputPP == NULL) {
		printf("\nError..\n");
		return 0;
	}

	FILE * outputPP2 = fopen ("outs/t10c_apm_pp2.txt", "a+");
	if (outputPP2 == NULL) {
		printf("\nError..\n");
		return 0;
	}

	FILE * outputAll = fopen ("outs/t10c_apm_all.txt", "a+");
	if (outputAll == NULL) {
		printf("\nError..\n");
		return 0;
	}

	FILE * outputIndependentExecs = fopen ("outs/independentExecs.txt", "a+");
	if (outputIndependentExecs == NULL) {
		printf("\nError..\n");
		return 0;
	}


	feasibleAverage = 0.;
	feasibleMedian = 0.;
	feasibleBest = INIT;
	feasibleWorst = -INIT;


	double sumFevalVtr = 0;
	int contVtr = 0;

	for (run = 0; run < MAX_RUN; run++) {	
		index = -1;
		int vtr = 0;

		swarm_t * swarm = allocate(); // Allocates the swarm.	

		// best
		best_particle_t * best;

		// printf("\n onde voce para\n");
		best = (best_particle_t * ) malloc(sizeof(best_particle_t));

		best->v = (float *) malloc(N_CON * sizeof(float));
		best->position = (double *) malloc(N_VAR * sizeof(double));
		best->fitness = INIT;
		best->fitnessAPM = INIT;

		for (m = 0; m < N_CON; m++)
			best->v[m] = INIT;
		for (m = 0; m < N_VAR; m++)
			best->position[m] = INIT;

		// Allocate
		violationAcum = (double **) malloc (MAX_PAR * sizeof(double * ));
		for (i = 0; i < MAX_PAR; i++)
			violationAcum[i] = (double *) malloc (N_CON * sizeof(double));
		k = (double * ) malloc (N_CON * sizeof (double));
	
		// Zero
		for(j = 0; j < MAX_PAR; j++)
			for(i = 0; i < N_CON; i++)
				violationAcum[j][i] = 0.;
		for(i = 0; i < N_CON; i++)
			k[i] = 0.;

		
		boundary(lowerBound, upperBound); // Boundary condition
		// printf("\n voltou aqui (pos boundary)\n");
		initialize(swarm, lowerBound, upperBound); // Initialize the particles
		// printf("\n vvoltou aqui (pos initialize)\n");

		fitness(swarm, k, violationAcum, MAX_TIME / 10, truss); // Fitness function	
		// printf("\n vvoltou aqui, pos fitness\n");

		// printf("saiu");
		// return 0;
		
		for (i = 0; i < MAX_TIME; i++) {
			alfa = (0.9 - 0.4) * ((float) (MAX_TIME - i) / MAX_TIME) + 0.4;

			calc_best(swarm); // Calculate pBest and gBest
			update(swarm, lowerBound, upperBound, alfa); // Update velocity and position
			fitness(swarm, k, violationAcum, i, truss); // Fitness function

			// best	
			for (m = 0; m < N_CON; m++) {
				if (swarm->gBest.v[m] > 0.0000001)
					break;
			}
			if (m == N_CON && swarm->gBest.fitnessAPM <= best->fitnessAPM) {
				index = i;
				for (m = 0; m < N_CON; m++)
					best->v[m] = swarm->gBest.v[m];
				for (m = 0; m < N_VAR; m++)
					best->position[m] = swarm->gBest.position[m];
				best->fitness = swarm->gBest.fitness;
				best->fitnessAPM = swarm->gBest.fitnessAPM;
			}

			if(best->fitnessAPM - valueToReach() <= 0){
				// printf("vtr obtained..\n");
				vtr = 1;
				// // printf("\n\n");
				// // printf("BEST %d - run #%d\n", index, run);
				// for (j = 0; j < N_VAR; j++){
				// 	printf("x[%d] = %.12lf\t\t", j, best->position[j]);
				// 	// bestFeasiblePosition[j] = best->position[j];
				// }	
				// printf("\nBest = %.12lf ", best->fitness);
				// printf("BestAPM = %.12lf ", best->fitnessAPM);
				// printf("\nsaiu if vtr");
				break;
			}


			// printf("continuou executando.\n");
			// Print gBest
			/*printf("\n\nIteração %d - ", i);
			for (j = 0; j < N_VAR; j++)	
				printf("x[%d] = %.12lf\t", j, swarm->gBest.position[j]);
			printf("\ngBest = %.12lf ", swarm->gBest.fitness);
			printf("gBestAPM = %.12lf ", swarm->gBest.fitnessAPM);
			for (j = 0; j < N_CON; j++)
				printf("v = %.12lf ", swarm->gBest.v[j]);*/
		}
			sumFevalVtr += index*MAX_PAR;
			contVtr+= 1;
			printf("\n\n");
			if(vtr){
				printf("vtr obtained..\n");
			}
			printf("BEST %d - feval %d - run #%d\n", index, index*MAX_PAR, run);
			for (j = 0; j < N_VAR; j++){
				printf("x[%d] = %.12lf\t\t", j, best->position[j]);
				// bestFeasiblePosition[j] = best->position[j];
			}	
			printf("\nBest = %.12lf ", best->fitness);
			printf("BestAPM = %.12lf ", best->fitnessAPM);
		
		// Imprime o melhor indíviduo de cada execução independente
		if (run == 0){ // Primeira iteração
			#ifdef USE_STRING
				fprintf (outputIndependentExecs, "STRING\n%lf ", best->fitness);
			#endif

			#ifdef USE_REDUCER
				fprintf (outputIndependentExecs, "REDUCER\n%lf ", best->fitness);
			#endif

			#ifdef USE_WELDED
				fprintf (outputIndependentExecs, "WELDED\n%lf ", best->fitness);
			#endif

			#ifdef USE_WELDED_MOD1_7
				fprintf (outputIndependentExecs, "WELDED_MOD1_7\n%lf ", best->fitness);
			#endif

			#ifdef USE_PRESSURE
				fprintf (outputIndependentExecs, "PRESSURE\n%lf ", best->fitness);
			#endif

			#ifdef USE_CANTILEVER
				fprintf (outputIndependentExecs, "CANTILEVER\n%lf ", best->fitness);
			#endif

			#ifdef USE_T10C
				fprintf (outputIndependentExecs, "T10C\n%lf ", best->fitness);
			#endif

			#ifdef USE_T10D
				fprintf (outputIndependentExecs, "T10D\n%lf ", best->fitness);
			#endif

			#ifdef USE_T25C
				fprintf (outputIndependentExecs, "T25C\n%lf ", best->fitness);
			#endif

			#ifdef USE_T25D
				fprintf (outputIndependentExecs, "T25D\n%lf ", best->fitness);
			#endif

			#ifdef USE_T60C
				fprintf (outputIndependentExecs, "T60C\n%lf ", best->fitness);
			#endif

			#ifdef USE_T72C
				fprintf (outputIndependentExecs, "T72C\n%lf ", best->fitness);
			#endif
		} else if (run == MAX_RUN -1){ // Última iteração
			fprintf (outputIndependentExecs, "%lf\tnRun: #%d\n", best->fitness, MAX_TIME * MAX_PAR);
		}
		else { // Demais iterações
			fprintf (outputIndependentExecs, "%lf ", best->fitness);
		}


		vectorAux[run] = best->fitnessAPM;
		if (best->fitnessAPM != INIT) {
			total++;
			feasibleAverage += best->fitnessAPM;
			if (best->fitnessAPM < feasibleBest){
				feasibleBest = best->fitnessAPM;
				for (j = 0; j < N_VAR; j++)
					bestFeasiblePosition[j] = best->position[j];
			}
			if (best->fitnessAPM > feasibleWorst)
				feasibleWorst = best->fitnessAPM;
			for (j = 0; j < N_VAR; j++)
				feasiblePosition[j] = best->position[j];
		}
		printf("\n");
		for (j = 0; j < N_CON; j++)
			printf("v = %e ", best->v[j]);

		printf("\n");
	

		//Deallocate
		for (i = 0; i < MAX_PAR; i++)	
			free (violationAcum[i]);
		
		free (violationAcum);
		free (k);
		// free (best);
		deallocate(swarm);
		// deallocateBest(best);
	}

	// Deallocate truss
	delete truss;

	// Calculate of the dp (desvio padrão)
	int cont, cont2 = 0;
	double auxDP = 0., vector[total];
	for (cont = 0; cont < MAX_RUN; cont++) {
		if (vectorAux[cont] != INIT) {
			auxDP += pow((vectorAux[cont] - feasibleAverage / total), 2);
			vector[cont2] = vectorAux[cont];
			cont2++;
			//printf("\nauxDP = %.14lf", auxDP);
		}
	}
	dp /= (total - 1);
	dp = sqrt(auxDP);

	// Printing (retirar)
	/*printf("\n\n");
	for(cont = 0; cont < MAX_RUN; cont++)
		printf("vectorAux[%d]: %lf\t", cont, vectorAux[cont] );
	printf("\n");
	for(cont2 = 0; cont2 < total; cont2++)
		printf("vector[%d]: %lf\t", cont2, vector[cont2]);*/

	// Calculate of the median 
	int number, term1, term2;
	double aux;
	printf("\n");
	for(cont = 1; cont < total; cont++) {
		aux = vector[cont];
		cont2 = cont - 1;
		while ((cont2 >= 0) && (aux < vector[cont2])) {
			vector[cont2 + 1] = vector[cont2];
			cont2--;
		}
		vector[cont2 + 1] = aux;
	}

	if ((total % 2) == 0) {
		term1 = total / 2;
		term2 = total / 2 + 1;
		feasibleMedian = (vector[term1 - 1] + vector[term2 - 1]) / 2;
		
	} else {
		number = (total + 1) / 2;
		feasibleMedian = vector[number - 1];
	}
	
	// Printing
	if (total != 0)	{
		//fprintf (output, "problema \t melhor \t mediana \t media \t pior \t desvio_padrao \t factibilidade\t nexecucao \n");
		//fprintf (outputAll, "problema \t melhor \t mediana \t media \t pior \t desvio_padrao \t variaveis \t factibilidade \t nexecucao \n");
		#ifdef USE_STRING
			fprintf (output, "STRING\t&\t %lf \t&\t %lf \t&\t %lf \t&\t %lf \t&\t %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			fprintf (outputPP, "STRING\t%lf \n", feasibleAverage / total);
			fprintf (outputPP2, "STRING\t%lf \n", feasibleBest);
			fprintf (outputAll, "STRING\t %lf \t %lf \t %lf \t %lf \t %e \t projectVariables(best individual): ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			for (j = 0; j < N_VAR; j++)
				fprintf (outputAll, "%lf \t", bestFeasiblePosition[j]);
			fprintf(outputAll, " endProjectVariables(best individual)\t");
			printf("\nSTRING\tBest: %lf Median: %lf Average: %lf Worst: %lf Desvio Padrão: %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
		#endif

		#ifdef USE_REDUCER
			fprintf (output, "REDUCER\t&\t %lf \t&\t %lf \t&\t %lf \t&\t %lf \t&\t %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			fprintf (outputPP, "REDUCER\t%lf \n", feasibleAverage / total);
			fprintf (outputPP2, "REDUCER\t%lf \n", feasibleBest);
			fprintf (outputAll, "REDUCER\t %lf \t %lf \t %lf \t %lf \t %e \t projectVariables(best individual): ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			for (j = 0; j < N_VAR; j++)
				fprintf (outputAll, "%lf \t", bestFeasiblePosition[j]);
			fprintf(outputAll, " endProjectVariables(best individual)\t");
			printf("\nREDUCER\tBest: %lf Median: %lf Average: %lf Worst: %lf Desvio Padrão: %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
		#endif

		#ifdef USE_WELDED
			fprintf (output, "WELDED\t&\t %lf \t&\t %lf \t&\t %lf \t&\t %lf \t&\t %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			fprintf (outputPP, "WELDED\t%lf \n", feasibleAverage / total);
			fprintf (outputPP2, "WELDED\t%lf \n", feasibleBest);
			fprintf (outputAll, "WELDED\t %lf \t %lf \t %lf \t %lf \t %e \t projectVariables(best individual): ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			for (j = 0; j < N_VAR; j++)
				fprintf (outputAll, "%lf \t", bestFeasiblePosition[j]);
			fprintf(outputAll, " endProjectVariables(best individual)\t");
			printf("\nWELDED\tBest: %lf Median: %lf Average: %lf Worst: %lf Desvio Padrão: %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
		#endif

		#ifdef USE_WELDED_MOD1_7
			fprintf (output, "WELDED_MOD1_7\t&\t %lf \t&\t %lf \t&\t %lf \t&\t %lf \t&\t %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			fprintf (outputPP, "WELDED_MOD1_7\t%lf \n", feasibleAverage / total);
			fprintf (outputPP2, "WELDED_MOD1_7\t%lf \n", feasibleBest);
			fprintf (outputAll, "WELDED_MOD1_7\t %lf \t %lf \t %lf \t %lf \t %e \t projectVariables(best individual): ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			for (j = 0; j < N_VAR; j++)
				fprintf (outputAll, "%lf \t", bestFeasiblePosition[j]);
			fprintf(outputAll, " endProjectVariables(best individual)\t");
			printf("\nWELDED_MOD1_7\tBest: %lf Median: %lf Average: %lf Worst: %lf Desvio Padrão: %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
		#endif

		#ifdef USE_PRESSURE
			fprintf (output, "PRESSURE\t&\t %lf \t&\t %lf \t&\t %lf \t&\t %lf \t&\t %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			fprintf (outputPP, "PRESSURE\t%lf \n", feasibleAverage / total);
			fprintf (outputPP2, "PRESSURE\t%lf \n", feasibleBest);
			fprintf (outputAll, "PRESSURE\t %lf \t %lf \t %lf \t %lf \t %e \t projectVariables(best individual): ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			for (j = 0; j < N_VAR; j++)
				fprintf (outputAll, "%lf \t", bestFeasiblePosition[j]);
			fprintf(outputAll, " endProjectVariables(best individual)\t");
			printf("\nPRESSURE\tBest: %lf Median: %lf Average: %lf Worst: %lf Desvio Padrão: %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
		#endif

		#ifdef USE_CANTILEVER
			fprintf (output, "CANTILEVER\t&\t %lf \t&\t %lf \t&\t %lf \t&\t %lf \t&\t %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			fprintf (outputPP, "CANTILEVER\t%lf \n", feasibleAverage / total);
			fprintf (outputPP2, "CANTILEVER\t%lf \n", feasibleBest);
			fprintf (outputAll, "CANTILEVER\t %lf \t %lf \t %lf \t %lf \t %e \t projectVariables(best individual): ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			for (j = 0; j < N_VAR; j++)
				fprintf (outputAll, "%lf \t", bestFeasiblePosition[j]);
			fprintf(outputAll, " endProjectVariables(best individual)\t");
			printf("\nCANTILEVER\tBest: %lf Median: %lf Average: %lf Worst: %lf Desvio Padrão: %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
		#endif

		#ifdef USE_T10C
			fprintf (output, "T10C\t&\t %lf \t&\t %lf \t&\t %lf \t&\t %lf \t&\t %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			fprintf (outputPP, "T10C\t%lf \n", feasibleAverage / total);
			fprintf (outputPP2, "T10C\t%lf \n", feasibleBest);
			fprintf (outputAll, "T10C\t %lf \t %lf \t %lf \t %lf \t %e \t projectVariables(best individual): ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			for (j = 0; j < N_VAR; j++)
				fprintf (outputAll, "%lf \t", bestFeasiblePosition[j]);
			fprintf(outputAll, " endProjectVariables(best individual)\t");
			printf("\nT10C\tBest: %lf Median: %lf Average: %lf Worst: %lf Desvio Padrão: %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
		#endif

		#ifdef USE_T10D
			fprintf (output, "T10D\t&\t %lf \t&\t %lf \t&\t %lf \t&\t %lf \t&\t %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			fprintf (outputPP, "T10D\t%lf \n", feasibleAverage / total);
			fprintf (outputPP2, "T10D\t%lf \n", feasibleBest);
			fprintf (outputAll, "T10D\t %lf \t %lf \t %lf \t %lf \t %e \t projectVariables(best individual): ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			for (j = 0; j < N_VAR; j++)
				fprintf (outputAll, "%lf \t", bestFeasiblePosition[j]);
			fprintf(outputAll, " endProjectVariables(best individual)\t");
			printf("\nT10D\tBest: %lf Median: %lf Average: %lf Worst: %lf Desvio Padrão: %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
		#endif

		#ifdef USE_T25C
			fprintf (output, "T25C\t&\t %lf \t&\t %lf \t&\t %lf \t&\t %lf \t&\t %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			fprintf (outputPP, "T25C\t%lf \n", feasibleAverage / total);
			fprintf (outputPP2, "T25C\t%lf \n", feasibleBest);
			fprintf (outputAll, "T25C\t %lf \t %lf \t %lf \t %lf \t %e \t projectVariables(best individual): ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			for (j = 0; j < N_VAR; j++)
				fprintf (outputAll, "%lf \t", bestFeasiblePosition[j]);
			fprintf(outputAll, " endProjectVariables(best individual)\t");
			printf("\nT25C\tBest: %lf Median: %lf Average: %lf Worst: %lf Desvio Padrão: %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
		#endif

		#ifdef USE_T25D
			fprintf (output, "T25D\t&\t %lf \t&\t %lf \t&\t %lf \t&\t %lf \t&\t %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			fprintf (outputPP, "T25D\t%lf \n", feasibleAverage / total);
			fprintf (outputPP2, "T25D\t%lf \n", feasibleBest);
			fprintf (outputAll, "T25D\t %lf \t %lf \t %lf \t %lf \t %e \t projectVariables(best individual): ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			for (j = 0; j < N_VAR; j++)
				fprintf (outputAll, "%lf \t", bestFeasiblePosition[j]);
			fprintf(outputAll, " endProjectVariables(best individual)\t");
			printf("\nT25D\tBest: %lf Median: %lf Average: %lf Worst: %lf Desvio Padrão: %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
		#endif

		#ifdef USE_T52C
			fprintf (output, "T52C\t&\t %lf \t&\t %lf \t&\t %lf \t&\t %lf \t&\t %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			fprintf (outputPP, "T52C\t%lf \n", feasibleAverage / total);
			fprintf (outputPP2, "T52C\t%lf \n", feasibleBest);
			fprintf (outputAll, "T52C\t %lf \t %lf \t %lf \t %lf \t %e \t projectVariables(best individual): ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			for (j = 0; j < N_VAR; j++)
				fprintf (outputAll, "%lf \t", bestFeasiblePosition[j]);
			fprintf(outputAll, " endProjectVariables(best individual)\t");
			printf("\nT52C\tBest: %lf Median: %lf Average: %lf Worst: %lf Desvio Padrão: %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
		#endif

		#ifdef USE_T52D
			fprintf (output, "T52D\t&\t %lf \t&\t %lf \t&\t %lf \t&\t %lf \t&\t %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			fprintf (outputPP, "T52D\t%lf \n", feasibleAverage / total);
			fprintf (outputPP2, "T52D\t%lf \n", feasibleBest);
			fprintf (outputAll, "T52D\t %lf \t %lf \t %lf \t %lf \t %e \t projectVariables(best individual): ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			for (j = 0; j < N_VAR; j++)
				fprintf (outputAll, "%lf \t", bestFeasiblePosition[j]);
			fprintf(outputAll, " endProjectVariables(best individual)\t");
			printf("\nT52D\tBest: %lf Median: %lf Average: %lf Worst: %lf Desvio Padrão: %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
		#endif

		#ifdef USE_T60C
			fprintf (output, "T60C\t&\t %lf \t&\t %lf \t&\t %lf \t&\t %lf \t&\t %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			fprintf (outputPP, "T60C\t%lf \n", feasibleAverage / total);
			fprintf (outputPP2, "T60C\t%lf \n", feasibleBest);
			fprintf (outputAll, "T60C\t %lf \t %lf \t %lf \t %lf \t %e \t projectVariables(best individual): ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			for (j = 0; j < N_VAR; j++)
				fprintf (outputAll, "%lf \t", bestFeasiblePosition[j]);
			fprintf(outputAll, " endProjectVariables(best individual)\t");
			printf("\nT60C\tBest: %lf Median: %lf Average: %lf Worst: %lf Desvio Padrão: %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
		#endif

		#ifdef USE_T72C
			fprintf (output, "T72C\t&\t %lf \t&\t %lf \t&\t %lf \t&\t %lf \t&\t %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			fprintf (outputPP, "T72C\t%lf \n", feasibleAverage / total);
			fprintf (outputPP2, "T72C\t%lf \n", feasibleBest);
			fprintf (outputAll, "T72C\t %lf \t %lf \t %lf \t %lf \t %e \t projectVariables(best individual): ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			for (j = 0; j < N_VAR; j++)
				fprintf (outputAll, "%lf \t", bestFeasiblePosition[j]);
			fprintf(outputAll, " endProjectVariables(best individual)\t");
			printf("\nT72C\tBest: %lf Median: %lf Average: %lf Worst: %lf Desvio Padrão: %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
		#endif

		#ifdef USE_T942C
			fprintf (output, "T942C\t&\t %lf \t&\t %lf \t&\t %lf \t&\t %lf \t&\t %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			fprintf (outputPP, "T942C\t%lf \n", feasibleAverage / total);
			fprintf (outputPP2, "T942C\t%lf \n", feasibleBest);
			fprintf (outputAll, "T942C\t %lf \t %lf \t %lf \t %lf \t %e \t projectVariables(best individual): ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
			for (j = 0; j < N_VAR; j++)
				fprintf (outputAll, "%lf \t", bestFeasiblePosition[j]);
			fprintf(outputAll, " endProjectVariables(best individual)\t");
			printf("\nT942C\tBest: %lf Median: %lf Average: %lf Worst: %lf Desvio Padrão: %e ", feasibleBest, feasibleMedian, feasibleAverage / total, feasibleWorst, dp);
		#endif
	}

	printf("Feasibility: %d/%d  nRun: #%d\n\n", total, MAX_RUN, MAX_TIME * MAX_PAR);
	printf("Mean FevalVTR: %.4f\n", sumFevalVtr/30);
	fprintf(output, "\t&\t %d/%d \t %d\n\n", total, MAX_RUN, MAX_TIME * MAX_PAR);
	fprintf(outputAll, "%d/%d \t%d\n\n", total, MAX_RUN, MAX_TIME * MAX_PAR);
	fclose(output);
	fclose(outputPP);
	fclose(outputPP2);
	fclose(outputAll);

	return 0;
}
