#include"user_funs.h"
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
	double Va = 5;
	double Vb = 1;
	double Tb = 10;
	matrix Y0 = matrix(3, new double[3]{ Va, Vb, Tb });
	double tend = 1000;
	double t0 = 0;
	double dt = 1;
	matrix* Y = solve_ode(df1, t0, dt, tend, Y0, ud1, x);
	int n = get_len(Y[0]);
	double max = Y[1](0, 2);
	for (int i = 1; i < n; ++i)
		if (max < Y[1](i, 2))
			max = Y[1](i, 2);
	y = abs(max - 50);
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}


matrix df1(double t, matrix Y, matrix ud1, matrix ud2) {
	matrix dY(3, 1);
	double g = 9.81;
	double b = 0.63;
	double a = 0.98;
	double P_A = 0.7; 
	double P_B = 1; 
	double D_B = 0.000365665;  
	double Fin = 0.01;
	double Ta = 90;
	double tin = 10;

	double F_out_A = Y(0, 0) > 0 ? -a * b * m2d(ud2) * sqrt(2.0 * g * Y(0, 0) / P_A) : 0.0;
	double F_out_B = Y(1, 0) > 0 ? - a * b * D_B * sqrt(2.0 * g * Y(1, 0) / P_B) : 0.0;

	dY(0, 0) = F_out_A;
	dY(1, 0) = -F_out_A + F_out_B + Fin;
	dY(2, 0) = (Fin / Y(1, 0)) * (tin - Y(2, 0)) - (F_out_A / Y(1, 0)) * (Ta - Y(2, 0));

	return dY;
}

matrix ff2T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = x(0) * x(0) + x(1) * x(1) - cos(2.5 * M_PI * x(0)) - cos(2.5 * M_PI * x(1)) + 2;
	return y;
}

matrix ff2T_2(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	y = 2.5 * pow((x(0) * x(0) - x(1)), 2)+pow((1-x(0)),2);
	return y;
}
matrix ff2R(matrix x, matrix ud1, matrix ud2) {
	matrix Y0 = matrix(2, 1);
	matrix Y_ref(2, new double[2] {3.14, 0});
	matrix y;
	double tend = 100;
	double t0 = 0;
	double dt = 0.1;
	matrix* Y = solve_ode(df2, t0, dt, tend, Y0, Y_ref, x);
	int n = get_len(Y[0]);
	for (int i = 0; i < n; i++) {
		y = y + 10 * pow(Y_ref(0) - Y[1](i, 0), 2) + pow(Y_ref(1) - Y[1](i, 1), 2) + pow(x(0) * (Y_ref(0) - Y[1](i, 0)) + x(1) * (Y_ref(1) - Y[1](i, 1)), 2);
	}
	y = y * dt;
	return y;
}

matrix df2(double t, matrix Y, matrix ud1, matrix ud2) {
	matrix dY(2, 1);
	double l = 0.6;
	double mr = 1;
	double mc = 9.5;
	double b = 0.5;
	double I;
	I = 1 / 3 * mr * l * l + mc * l * l;
	//dY[0,0] = Y(1);
	//dY[1,0] = ((ud2(0) * (ud1(0) - Y(0)) + ud2(1) * (ud1(1) - Y(1)) - b * Y(1))) / I;
	dY(0, 0) = Y(1, 0);
	dY(1, 0) = ((ud2(0) * (ud1(0) - Y(0, 0)) + ud2(1) * (ud1(1) - Y(1, 0)) - b * Y(1, 0))) / I;

	return dY;
}