#pragma once

#include"ode_solver.h"

// LAB 0
matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);

// LAB 1
matrix ff1T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix ff1R(matrix, matrix ud1 = NAN, matrix ud2 = NAN);
matrix df1(double t, matrix, matrix ud1 = NAN, matrix ud2 = NAN);

// LAB 2
matrix ff2T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix ff2T_2(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix df2(double t, matrix Y, matrix ud1 = NAN, matrix ud2 = NAN);
matrix ff2R(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

// LAB 3
matrix ff3T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix ff3T_2(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix g1(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix g2(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix g3(matrix x, matrix a, matrix ud1 = NAN, matrix ud2 = NAN);
matrix S(matrix x, matrix a, matrix ud1 = NAN, matrix ud2 = NAN);
matrix F_zewn(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

//LAB 4
matrix ff4T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix ff4R(matrix x, matrix ud1, matrix ud2);
matrix gf4(matrix x, matrix ud1, matrix ud2);