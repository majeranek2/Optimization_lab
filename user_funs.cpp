﻿#include"user_funs.h"
#include <cmath>
#include <corecrt_math_defines.h>

matrix ff0T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
	int n = get_len(Y[0]);
	double teta_max = Y[1](0, 0);
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);
	double m = 1, l = 0.5, b = 0.5, g = 9.81;
	double I = m * pow(l, 2);
	dY(0) = Y(1);
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;
	return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = -cos(0.1 * m2d(x)) * exp(-(0.1 * m2d(x) - 2 * M_PI) * (0.1 * m2d(x) - 2 * M_PI)) + 0.002 * (0.1 * x) * (0.1 * x);
	return y;
}

matrix ff1R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(3, new double[3]{ 5, 1, 10 });
	matrix* Y = solve_ode(df1, 0, 1, 1000, Y0, ud1, x);
	int n = get_len(Y[0]);
	double max = Y[1](0, 2);
	for (int i = 1; i < n; ++i)
		if (max < Y[1](i, 2))
			max = Y[1](i, 2);
	y = abs(max - 50);
	return y;
}


matrix df1(double t, matrix Y, matrix ud1, matrix ud2) {
	matrix dY(3, 1);
	// Stałe i właściwości
	double g = 9.81;
	double b = 0.63;
	double a = 0.98;
	double P_A = 0.7;  // Pole przekroju zbiornika A w mm²
	double P_B = 1;  // Pole przekroju zbiornika B w mm²
	double D_B = 0.000365665;  // Pole przekroju wylotu zbiornika B w mm²

	// Obliczenie strumieni wypływu
	double F_out_A = Y(0, 0) > 0 ? a * b * m2d(ud2) * sqrt(2.0 * g * Y(0, 0) / P_A) : 0.0;
	double F_out_B = Y(1, 0) > 0 ? a * b * D_B * sqrt(2.0 * g * Y(1, 0) / P_B) : 0.0;

	// Dodatkowe wartości (do zastąpienia rzeczywistymi danymi)
	double Fin = 0.01;  // Przepływ do zbiornika B w mm³/s
	double Ta = 90;   // Temperatura napływającej wody w stopniach Celsiusza
	double tin = 10;  // Wartość czasowa w sekundach

	// Równania różniczkowe
	dY(0, 0) = F_out_A;
	dY(1, 0) = -F_out_A + F_out_B + Fin;
	dY(2, 0) = (Fin / Y(1, 0)) * (tin - Y(2, 0)) - (F_out_A / Y(1, 0)) * (Ta - Y(2, 0));

	return dY;
}

//matrix df1(double t, matrix Y, matrix ud1, matrix ud2) {
//	matrix dY(3, 1);
//	// Stałe i właściwości
//	double g = 9.81;
//	double b = 0.63;
//	double a = 0.98;
//	double P_A = 700000;  // Pole przekroju zbiornika A w mm²
//	double P_B = 1000000;  // Pole przekroju zbiornika B w mm²
//	double D_B = 36566.5;  // Pole przekroju wylotu zbiornika B w mm²
//
//	// Obliczenie strumieni wypływu
//	double F_out_A = Y(0, 0) > 0 ? -a * b * m2d(ud2) * sqrt(2.0 * g * Y(0, 0) / P_A) : 0.0;
//	double F_out_B = Y(1, 0) > 0 ? -a * b * D_B * sqrt(2.0 * g * Y(1, 0) / P_B) : 0.0;
//
//	// Dodatkowe wartości (do zastąpienia rzeczywistymi danymi)
//	double Fin = 10;  // Przepływ do zbiornika B w mm³/s
//	double Ta = 90;   // Temperatura napływającej wody w stopniach Celsiusza
//	double tin = 10;  // Wartość czasowa w sekundach
//
//	// Równania różniczkowe
//	dY(0, 0) =  F_out_A;
//	dY(1, 0) = -F_out_A + F_out_B + Fin;
//	dY(2, 0) = (Fin / Y(1, 0)) * (tin - Y(2, 0)) + (F_out_A / Y(1, 0)) * (Ta - Y(2, 0));
//
//	return dY;
//}
