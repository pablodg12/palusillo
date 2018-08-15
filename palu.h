#define palu_h
#include <stdio.h>

/*Palu, esta función recibe por inputos un arreglo unidimensional de doubles*/

/*Esta función es necesaria para generar la matriz identidad inicial*/
double* generate_identity(int size);

/*Función para generar la matriz de 0 inicial */
double* generate_zeros(int size);

/*Función para encontrar el argmax, luego de una permutación ignora la fila permutada anterior */
int argmax(int size, double* A,int j);

/*Función que realiza las permutaciones*/
double* row_perm(int size, double* A, int i, int j);

/*Función General de la descomposición Palu */
double* palu_decomp(double* A, double* b,int size);

/*Función para realizar el producto punto */
double* dot(int size, double* P, double* b);

/*Función que resuelve las matrices triangulares superiores e inferiores */
double* solve_triangular(int size, int upper, double* A, double* b);