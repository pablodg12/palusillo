#define palu_h
#include <stdio.h>

/*Palu, esta función recibe por inputos un arreglo unidimensional de doubles*/
/*Este codigo esta basado en la implementación del profesor Claudio Torres (tclaudioe) */

/*Función para encontrar el argmax, luego de una permutación ignora la fila permutada anterior */
int argmax(int size, double* A,int j);

/*Función que realiza las permutaciones*/
double* row_perm(int size, double* A, int i, int j);

/*Función General de la descomposición Palu */
double* solve_palu(double* A, double* b,int size);

/*Función que resuelve las matrices triangulares superiores e inferiores */
double* solve_ALU(int size, double* A, double* b);