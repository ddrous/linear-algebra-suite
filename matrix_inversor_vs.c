#include <stdio.h>
#include <stdlib.h>

//#include "stdafx.h"
#include <math.h>
void Afficher(float **X, int n)
{
	int l, m, p;
	for (l = 0; l<n; l++)
	{
		for (m = 0; m<n; m++)
			printf(" %f\t", X[l][m]);
		printf("|");
		for (p = n; p<2 * n; p++)
			printf("\t%f", X[l][p]);
		printf("\n\n");
	}
}
void Afficher2(float **X, int n)
{
	int l, m;
	for (l = 0; l<n; l++)
	{
		for (m = 0; m<n; m++)
			printf(" %f\t", X[l][m]);
		printf("\n\n");
	}
}
void produit_mat(float **A, float **B, int n)
{
	float **C;
	C = (float**)malloc(n * sizeof(float));
	for (int i = 0; i<n; i++)
	{
		C[i] = (float*)malloc((n) * sizeof(float));
	}
	for (int i = 0; i<n; i++)
	{
		for (int j = 0; j<n; j++)
		{
			C[i][j] = 0;
			for (int k = 0; k<n; k++)
			{
				C[i][j] = C[i][j] + A[i][k] * B[k][j];
			}
		}
	}
	Afficher2(C, n);
	for (int i = 0; i<n; i++)
	{
		free(C[i]);
	}
	free(C);
}
int main()
{
	int n, i, j, k, q, m, s, u, t = 1, r = 0, flag = 1, flag2 = 1, flag3 = 1, flag4 = 1, flag5 = 1;
	float pivot, aux, aux2, det = 1, **A, **A_prime, **B;
	printf("\n BIENVENUE DANS LE PROGRAMME DE CALCUL DE L'INVERSE D'UNE MATRICE CARREE PAR LA METHODE DE GAUSS-JORDAN\n");
	printf("\n Entrez la taille n de la matrice A sachant que n doit etre superieur ou egal a 1 :\t");
	scanf("%d", &n);
	printf("\n");
	A = (float**)malloc(n * sizeof(float));
	A_prime = (float**)malloc(n * sizeof(float));
	B = (float**)malloc(n * sizeof(float));
	for (i = 0; i<n; i++)
	{
		A[i] = (float*)malloc((2 * n) * sizeof(float));
		A_prime[i] = (float*)malloc(n * sizeof(float));
		B[i] = (float*)malloc((n) * sizeof(float));
	}
	for (i = 0; i<n; i++)
	{
		for (j = 0; j<n; j++)
		{
			printf("\n Entrez le coefficient A[%d][%d]\t\t", i + 1, j + 1);
			scanf("%f", &A[i][j]);
			A_prime[i][j] = A[i][j];
		}
		for (q = n; q<2 * n; q++)
			if (q == n + i)
				A[i][q] = 1;
			else
				A[i][q] = 0;
		printf("\n");
	}
	printf("\n");
	printf("\n ====> Matrice A <====\n");
	printf("\n");
	Afficher2(A, n);
	printf("\n");
	printf("\n Voici la matrice augmentee de A:\n");
	printf("\n");
	Afficher(A, n);
	//if((A[0][0]==1) || (A[0][0]!= 0))
	//  flag2=1;
	for (k = 0; k<n - 1; k++)
	{
		i = k;
		pivot = 0;
		while ((pivot == 0) && (i<n))
		{
			pivot = A[i][k];
			i = i + 1;
		}
		if ((pivot == 0) || (pivot == -0))
			det = 0;
		else
		{
			if (pivot != A[k][k])
			{
				r = r + 1;
				flag2 = 0;
			}
			for (j = k; j<2 * n; j++)
			{
				aux = A[k][j]; A[k][j] = A[i - 1][j]; A[i - 1][j] = aux; //Permutation pour trouver le nouveau pivot
			}
			if (flag2 == 0)
			{
				printf("\n\n Apres l'etape %d, la matrice augmentee de A devient:\n\n", t);
				Afficher(A, n);
				t = t + 1;
			}
			det = det*A[k][k];
			m = 2 * n + k - 1;
			if (A[k][k] == 1)
				flag3 = 0;
			for (j = k; j<2 * n; j++)
				A[k][m - j] = A[k][m - j] / A[k][k];
			if (flag3 == 1)
			{
				printf("\n\n Apres l'etape %d, la matrice augmentee de A devient:\n\n", t);
				Afficher(A, n);
				t = t + 1;
			}
			for (i = k + 1; i<n; i++)
			{
				flag4 = 1;
				if (A[i][k] == 0)
					flag4 = 0;
				for (j = k; j<2 * n; j++)//On commence par l'elements de la case A[i][2*n-1] pour ne pas changer la valeu de A[k][k]
					A[i][m - j] = A[i][m - j] - A[i][k] * A[k][m - j];//Elimination de Gauss proprement dite
			}
			if (flag4 == 1)
			{
				printf("\n\n Apres l'etape %d, la matrice augmentee de A devient:\n\n", t);
				Afficher(A, n);
				t = t + 1;
			}
		}
		flag2 = 1;
		flag3 = 1;
	}
	det = det*A[n - 1][n - 1];
	if (det == 0)
	{
		printf("\n La matice A n'est pas inversible et son determinant est nul\n");
		flag = 0;
	}
	else
	{
		if (A[n - 1][n - 1] != 1)
		{
			s = 3 * n - 2;
			for (j = n - 1; j<2 * n; j++)
				A[n - 1][s - j] = A[n - 1][s - j] / A[n - 1][n - 1];
			printf("\n\n Apres l'etape %d, la matrice augmentee de A devient:\n\n", t);
			Afficher(A, n);
			t = t + 1;
		}
	}
	if (flag == 1)
	{
		for (k = 0; k<n - 1; k++)//Debut de la remontee
		{
			for (i = k; i<n - 1; i++)
			{
				flag5 = 1;
				aux2 = A[n - 2 - i][n - 1 - k];
				if (aux2 == 0)
					flag5 = 0;
				for (j = 0; j<2 * n; j++)
				{
					A[n - 2 - i][j] = A[n - 2 - i][j] - A[n - 1 - k][j] * aux2;
				}
			}
			if (flag5 == 1)
			{
				printf("\n\n Apres l'etape %d, la matrice augmentee de A devient:\n\n", t);
				Afficher(A, n);
				t = t + 1;
			}
		}
		for (i = 0; i<n; i++)
		{
			for (j = 0; j<n; j++)
				B[i][j] = A[i][n + j];
		}
		u = pow(-1, r);
		det = u*det;
		printf("\n\n La matrice inverse de A est donnee par Inv(A) comme suit:\n\n");
		printf("\n ====> Matrice Inv(A) <====\n");
		Afficher2(B, n);
		printf("\n\n Et son determinant est donne par Det(A) = %f", det);
		printf("\n\n Juste pour verification, effectuons le produit de A par Inv(A) \n");
		printf("\n ====> Matrice A <====\n");
		Afficher2(A_prime, n);
		printf("\n ====> Produit Matriciel A x Inv(A) <====\n");
		produit_mat(A_prime, B, n);
	}
	else
		printf("\n Veillez entrez une matrice inversible. A vous de determiner les criteres d'inversibilites");
	for (i = 0; i<n; i++)
	{
		free(A[i]);
		free(B[i]);
	}
	free(A);
	free(B);
	printf("\n\n ====>Programmeur 'GenFreak' <====\n\t");
	
	getchar();
	getchar();
	return 0;
}
