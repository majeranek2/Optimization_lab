//Ten plik nie powinien by� edytowany

#pragma once

#include"solution.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN);

solution expansion(matrix(*ff)(matrix, matrix, matrix),  double x0,  double d,  double alpha, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
solution fib(matrix(*ff)(matrix, matrix, matrix), double A,  double B, vector< double> fi,  double epsilon, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
solution lag(matrix(*ff)(matrix, matrix, matrix),  double a,  double b,  double epsilon,  double gamma, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN);
solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
matrix HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, matrix alpha, matrix beta, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
matrix add_scalar_to_matrix_row(matrix& m, matrix scalar, int row);
matrix multiply_matrix_row_by_scalar(matrix& m, matrix scalar, int row);
matrix getBase(int n);
matrix Rosen_getQ(int n, matrix lambda, matrix d);
matrix Rosen_getVj(int n, int j, matrix Q, matrix d);

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double alpha, double beta, double gamma, double delta, double s, double c, double dc, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
matrix get_simplex(int n, double s = 1.0);
int get_min(matrix p, int n);
int get_max(matrix p, int n, int pmin);

//lab4
solution SD_1(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
solution SD_2(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN);
solution CG_1(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
solution CG_2(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN);
solution Newton_1(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
solution Newton_2(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
matrix df_dx1(matrix x);
matrix df_dx2(matrix x);
matrix f1(matrix x, matrix ud1, matrix ud2);
solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
