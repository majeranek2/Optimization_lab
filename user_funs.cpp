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

matrix ff3T_2(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	y = pow(x(0,0), 2) + pow(x(1,0), 2);
	// z ksiazki 3.29
	//y = pow(x(0, 0), 2) + pow(x(1, 0), 2) - ( (400.0 / ((100.0 * pow(x(0, 0), 2)) + pow(x(1, 0), 2) + 1.0)));

	// z ksiazki 3.30
	//y = 2.5 * pow(pow(x(0,0),2) - x(1,0), 2) + pow(1 - x(0, 0), 2);
	return y;
}

//matrix ff3T(matrix x, matrix ud1, matrix ud2) {
//	matrix y;
//	y = sin(M_PI * sqrt(pow(x(0, 0) / M_PI, 2) + pow(x(1, 0) / M_PI, 2)))
//		/ M_PI * sqrt(pow(x(0, 0) / M_PI, 2) + pow(x(1, 0) / M_PI, 2));
//	return y;
//}

//double licz(double v0x, double omg) 
//{
//	double x = 0;
//	double vx = v0x;
//	double vy = 0;
//	double x_prev = 0;
//	double m = 0.6;	//kg = 600 g
//	double r = 0.12;	//m=12 cm
//	double y0 = 100;	//m
//	double S = M_PI * r * r;
//	double ro = 1.2;	//kg/m3
//	double C = 0.47;
//	double t0 = 0;
//	double dt = 0.01;
//	double tk = 7; 
//	double g = 9.89;
//	double y, ax, ay;
//
//	while (abs(x - x_prev) >= 1e-6)
//	{
//		x += vx * dt;
//		y = y0 - vy * dt;
//		ax = (-0.5 * C * ro * S * vx * vx + ro * vy + omg* M_PI * r ^3) / m;
//		ay = (-m * g - 0.5 * C * ro * S * vy * vy - ro * vx - omg* M_PI * r^3) / m;
//		vx += ax * dt;
//		vy += ay * dt;
//		if (abs(x - x_prev) < 1e-6) {
//			break;
//		}
//		x_prev = x;
//	}
//	return x;
//
//}
//double znajdz()
//{
//	double max_xend = 0;
//	double best_v0x = 0;
//	double best_omega = 0;
//	double xend;
//
//	for (int v0x = -10; v0x < 11; v0x++) {
//		for (int omega = -23; omega < 24; omega++) {
//			xend = licz(v0x, omega);
//			if (4 <= xend <= 6 && abs(xend - 5) < abs(max_xend - 5)) {
//				max_xend = xend;
//				best_v0x = v0x;
//				best_omega = omega;
//				return max_xend, best_v0x, best_omega;
//
//				
//			}
//		}
//	}
//	return max_xend, best_v0x, best_omega;
//}


//matrix ff3R(matrix x, matrix ud1, matrix ud2)
//{
//	double max_xend = 0;
//	double best_v0x = 0;
//	double best_omega = 0;
//	max_xend, best_v0x, best_omega = znajdz();
//
//}
//matrix g1(matrix x, matrix ud1, matrix ud2) {
//	matrix y;
//	y = (x(0, 0)) + 1;
//	return y;
//}
//
//matrix g2(matrix x, matrix ud1, matrix ud2) {
//	matrix y;
//	y = (x(1, 0)) + 1;
//	return y;
//}
//
//matrix g3(matrix x, matrix a, matrix ud1, matrix ud2) {
//	matrix y;
//	y = sqrt(pow(x(0, 0), 2) + pow(x(1, 0), 2)) - a;
//	return y;
//}
//
//matrix S(matrix x, matrix a, matrix ud1, matrix ud2) {
//	matrix y;
//	matrix ZERO(0.0);
//	y = pow(max(ZERO, g1(x)), 2) +
//		pow(max(ZERO, g2(x)), 2) +
//		pow(max(ZERO, g3(x, a)), 2);
//	return y;
//}

matrix ff3T(matrix x, matrix ud1, matrix ud2) { // ud1: a, ud2: c [int 2, 2]
	double Y = sin(M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2))) 
		/ (M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2)));
	matrix g1(-x(0) + 1.0);
	matrix g2(-x(1) + 1.0);
	matrix g3(norm(x) - ud1(0));
	matrix ZERO(0.0);
	matrix y(Y);
	if (-x(0) + 1 > 0) {
		y = y + ud2(0) * pow(max(ZERO, g1), 2);
	}
	if (-x(1) + 1 > 0) {
		y = y + ud2(0) * pow(max(ZERO, g2), 2);
	}
	if (norm(x) - ud1(0) > 0) {
		y = y + ud2(0) * pow(max(ZERO, g3), 2);
	}
	return y;
}

//lab 4
matrix ff4T(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	y = pow((x(0) + 2 * x(1) - 7), 2) + pow((2 * x(0) + x(1) - 5), 2);
	return y;
}


matrix ff4R(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	int m = 100;	//liczba danych
	int n = get_len(x);		//liczba współrzędnych wektora gradientu
	static matrix X(n, m), Y(1, m);

	//wczytanie danych z pliku
	static bool pomm = true;
	if (pomm) {
		std::ifstream plik("XData.txt");
		if (!plik.is_open()) {
			std::cerr << "Error XData.txt" << std::endl;
			return 0;
		}
		plik >> X;
		plik.close();

		plik.open("YData.txt");
		if (!plik.is_open()) {
			std::cerr << "Error YData.txt" << std::endl;
			return 0;
		}
		plik >> Y;
		plik.close();

		pomm = false;
	}

	double h, pom;
	for (int i = 0; i < m; i++) {
		pom = (trans(x) * X[i])(); //wktor ztransponowany z paramtrami klasyfiktora * wktor ocen 
		h = 1.0 / (1.0 + exp(-pom));
		y = y - Y(0, i) * log(h) - (1 - Y(0, i)) * log(1 - h);
	}

	y = y / m;
	return y;
}





matrix gf4(matrix x, matrix ud1, matrix ud2) {
	int m = 100;	//liczba danych	
	int n = get_len(x);//liczba współrzędnych wektora gradientu
	matrix g(n, 1);
	static matrix X(n, m), Y(1, m);

	//wczytanie danych z pliku
	static bool pomm = true;
	if (pomm) {
		std::ifstream plik("XData.txt");
		if (!plik.is_open()) {
			std::cerr << "Error XData.txt" << std::endl;
			return 0;
		}
		plik >> X;
		plik.close();

		plik.open("YData.txt");
		if (!plik.is_open()) {
			std::cerr << "Error YData.txt" << std::endl;
			return 0;
		}
		plik >> Y;
		plik.close();

		pomm = false;
	}
	double h, pom;
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < m; i++) {
			pom = (trans(x) * X[i])();
			h = 1.0 / (1.0 + exp(-pom));
			g(j) = g(j) + X(j, i) * (h - Y(0, i));
		}
		g(j) = g(j) / m;
	}
	return g;
}