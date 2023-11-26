#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

solution expansion(matrix(*ff)(matrix, matrix, matrix),  double x0,  double d,  double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		ofstream Sout("wyniki_exp.csv");
		//double* p = new  double[2]{ 0, 0 };
		int i = 0;
		vector<double> x;
		 x.push_back(x0);			// x0
		 x.push_back(x[0] + d);	// x1
		if (ff1T(x[1]) == ff1T(x[0]))
		{
			Xopt.x = x[0];
			Xopt.y = x[1];
			Xopt.f_calls = i;
			Sout << hcat(x[0], Xopt.x);
			Sout << hcat(Xopt.y, i) << "\n";
			return Xopt;
		}
		else if (ff1T(x[1]) > ff1T(x[0]))
		{
			d = -d;
			x[1] = x0 + d;
			if (ff1T(x[1]) >= ff1T(x[0]))
			{
				Xopt.x = x[1];
				Xopt.y = x[0] - d;
				Xopt.f_calls = i;
				Sout << hcat(x[0], Xopt.x);
				Sout << hcat(Xopt.y, i) << "\n";
				return Xopt;
			}
		}

		while (ff1T(x[i]) > ff1T(x[i + 1]))
		{
			if (i > Nmax)
				throw "Iteration number exceeded!\n";
			i++;
			x.push_back(x[0] + pow(alpha, i) * d);
		}

		if (d > 0)
		{
			Xopt.x = x[i - 1];
			Xopt.y = x[i + 1];
			Xopt.f_calls = i;
			Sout << hcat(x[0], Xopt.x);
			Sout << hcat(Xopt.y, i) << "\n";
			return Xopt;
		}

		Xopt.x = x[i + 1];
		Xopt.y = x[i - 1];
		Xopt.f_calls = i;
		Sout << hcat(x[0], Xopt.x);
		Sout << hcat(Xopt.y, i) << "\n";
		Sout.close();
		return Xopt;

	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a,  double b, vector< double> fi,  double epsilon, matrix ud1, matrix ud2) {
	try {
		solution Xopt;
		ofstream Sout("wyniki_fib.csv");
		int k = 2;
		while (fi[k] <= ((b - a) / epsilon)) {
			k++;
		}

		 double c = b - (fi[k - 2] / fi[k - 1]) * (b - a);
		 double d = a + b - c;
		for (int i = 0; i < k - 2; i++)
		{
			if (ff(c, ud1, ud2) < ff(d, ud1, ud2))
				b = d;

			else
				a = c;

			c = b - (fi[k - i - 3] / (fi[k - i - 2])) * (b - a);
			d = a + b - c;
			Sout << hcat(i, c) << "\n";
		}
		Sout.close();

		Xopt.x = c;
		Xopt.y = Xopt.fit_fun(ff, ud1, ud2);
		Xopt.f_calls = k;
		return Xopt;
	}
	catch (string ex_info) {
		throw ("solution fib(...):\n" + ex_info);
	}
}


solution lag(matrix(*ff)(matrix, matrix, matrix),  double a,  double b,  double epsilon,  double gamma, int Nmax, matrix ud1, matrix ud2) {
	try
	{
		solution Xopt;
		int i = 0;
		ofstream Sout("wyniki_lag.csv");
		double c =  a + (b - a) / 3;
		 double L, M;
		 double d;
		 double f_a, f_b, f_c, f_d;

		do {
			f_a = ff(a, ud1, ud2)(0,0);
			f_b = ff(b, ud1, ud2)(0, 0);
			f_c = ff(c, ud1, ud2)(0, 0);

			L = f_a * (pow(b, 2)- pow(c, 2))+ f_b * (pow(c, 2)- pow(a, 2)) + f_c * (pow(a, 2)- pow(b, 2));
			M = f_a * (b - c)+ f_b * (c - a)+ f_c * (a - b);

			if (M <= 0) {
				cout << "Blad - mianownik niedodatni!\n";
				return NULL;
			}

			d = 0.5 * (L / M);
			f_d = ff(d, ud1, ud2)(0, 0);

			if (a < d && d < c) {
				if (f_d < f_c) {
					b = c;
					c = d;
				}
				else {
					a = d;
				}
			}

			else {
				if (c < d && d < b) {
					if (f_d < f_c) {
						a = c;
						c = d;
					}
					else {
						b = d;
					}
				}
			}
			Sout << hcat(i, d) << "\n";
			i++;
			if (i > Nmax)
				break;

		} while ((b - a) >= epsilon || abs(d - c) > gamma);

		Sout.close();
		Xopt.x = d;
		Xopt.y = ff(d, ud1, ud2);
		Xopt.f_calls = i;
		return Xopt;
	}

	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}

}



solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		matrix x = x0;
		matrix xB = x0;
		matrix xB_;
		int i = 0;

		while (s >= epsilon)
		{
			xB = x;
			x = HJ_trial(ff, xB, s, ud1, ud2);

			if (ff(x, ud1, ud2)(0, 0) < ff(xB, ud1, ud2)(0, 0))
			{
				while (ff(x, ud1, ud2)(0, 0) < ff(xB, ud1, ud2)(0, 0))
				{
					xB_ = xB;
					xB = x;
					x = 2 * xB - xB_;
					x = HJ_trial(ff, xB, s, ud1, ud2);
					i++;
					if (i > Nmax)
						break;
				}
				x = xB;
			}
			else
				s = alpha * s;

			if (i > Nmax)
				break;
		}

		Xopt.x = xB;
		Xopt.y = ff(Xopt.x, ud1, ud2)(0, 0);
		Xopt.f_calls = i;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

matrix  HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		matrix x = XB.x;
		matrix x_plus, x_minus;

		for (int j = 0; j < get_size(x)[0]; j++)
		{
			x_plus = x;
			x_plus(j, 0) += s;
			x_minus = x;
			x_minus(j, 0) -= s;

			if (ff(x_plus, ud1, ud2)(0, 0) < ff(x, ud1, ud2)(0, 0))
				x = x_plus;
			else
			{
				if (ff(x_minus, ud1, ud2)(0, 0) < ff(x, ud1, ud2)(0, 0))
					x = x_minus;
			}
		}

		return x;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		matrix x_b = x0;
		matrix s = s0;

		int n = get_len(x0[0]);
		double* ZERO = new double[n] {};

		matrix lambda(n, ZERO);
		
		matrix p(n, ZERO);
		matrix d = getBase(n);

	/*	int* size_xb = get_size(x_b[0]);
		cout << "Size x_b: " << size_xb[0] << ", " << size_xb[1] << endl;
		int* size_s = get_size(s[0]);
		cout << "Size s: " << size_s[0] << ", " << size_s[1] << endl;
		int* size_d = get_size(d[0]);
		cout << "Size d: " << size_d[0] << ", " << size_d[1] << endl;*/
	


		int i = 0;
		bool przerwa = true;

		do {
			for (int j = 0; j < n; j++) {
				matrix S = get_row(s, j);
				x_b = x0[i];
				// === DEBUG
				/*cout << "x_b:\n" << x_b << endl;
				cout << "S:\n" << S << endl;
				cout << "d:\n" << d << endl;*/
				// ===
				//
				//if (ff(x_b + (s[j] * d[j]), ud1, ud2) < ff(x_b, ud1, ud2)) {	// Tutaj wymiary macierzy siê nie zgadzaj¹...
				//
				//	x_b = x_b + (s[j] * d[j]);
				//	lambda[j] = lambda[j] + s[j];
				//	s[j] = s[j] * alpha;
				//}
				//else {
				//	s[j] = s[j] * (-beta);
				//	p[j] = p[j] + 1;
				//}

				if (ff(x_b + (S * d[j]), ud1, ud2) < ff(x_b, ud1, ud2)) {	// Tutaj wymiary macierzy siê nie zgadzaj¹...

					x_b = x_b + (S * d[j]);
					lambda[j] = lambda[j] + S;
					S = S * alpha;
				}
				else {
					for (int i = 0; i < n; i++) {
						S[i] = (-beta) * S[i];
					}
					//S = (-beta) * S;
					matrix p_ = get_row(p, j);
					p.set_row(p_, j);
				}
			}
			i++;

			bool zmiana_bazy = true;
			for (int j = 0; j < n; j++) {
				if (get_row(lambda, j) == 0.0 || get_row(p,j) == 0.0) {
					zmiana_bazy = false;
					break;
				}
			}

			if (zmiana_bazy) {
				matrix Q = Rosen_getQ(n, lambda, d);
				matrix vj(n, n, 0.0);
				vj.set_col(Q[0], 0);
				for (int j = 1; j < n; j++) {
					vj.set_col(Rosen_getVj(n, j, Q, d), j);
				}

				for (int j = 0; j < n; j++) {
					//d[j] = (vj[j] / norm(vj[j]));
					matrix temp = vj[j] / norm(vj[j]);
					d.set_col(temp, j);
				}

				matrix lambda_(n, ZERO);
				lambda = lambda_;
				matrix p_(n, ZERO);
				p = p_;
				s = s0;
			}

			if (i > Nmax) {
				przerwa = false;
			}
			for (int j = 0; j < n; j++) {
				double a = s(j);
				if (abs(a) > epsilon) {
					przerwa = false;
				}
			}

		} while (przerwa == true);

		Xopt.x = x_b;
		Xopt.y = ff(Xopt.x, ud1, ud2)(0, 0);
		Xopt.f_calls = i;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

matrix getBase(int n) {
	double** base = new double* [n];
	for (int i = 0; i < n; ++i) {
		base[i] = new double[n];
	}

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (j == i) {
				base[i][j] = 1.0; 
			}
			else {
				base[i][j] = 0.0; 
			}
		}
	}
	matrix baza(n, n, base);
	return baza;
}

matrix Rosen_getQ(int n, matrix lambda, matrix d) {
	double** tab = new double* [n];
	for (int i = 0; i < n; ++i) {
		tab[i] = new double[n];
	}

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (j > i) {
				tab[i][j] = 0.0;
			}
			else {
				tab[i][j] = m2d(lambda[i]);
			}
		}
	}
	matrix Q_l(n, n, 0.0);
	for (int j = 0; j < n; j++) {
		Q_l.set_col(d[j], j);
	}
	matrix Q_r(n, n, tab);
	return Q_l * Q_r;
}

matrix Rosen_getVj(int n, int j, matrix Q, matrix d) {
	double* ZERO = new double[n] {};
	matrix suma(n, ZERO);
	for (int k = 1; k <= j - 1; k++) {
		suma = suma + ((trans(Q[j]) * d[k]) * d[k]);
	}
	matrix result = Q[j] - suma;
	return result;
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
