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


