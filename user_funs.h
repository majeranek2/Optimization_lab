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

