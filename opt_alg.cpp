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

long double* expansion(matrix(*ff)(matrix, matrix, matrix), long double x0, long double d, long double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		long double* p = new long double[2]{ 0, 0 };
		int i = 0;
		vector<double> x;
		x.push_back(x0);			// x0
		x.push_back(x[0] + d);	// x1
		if (ff1T(x[1]) == ff1T(x[0]))
		{
			p[0] = x[0];
			p[1] = x[1];
			return p;
		}
		else if (ff1T(x[1]) > ff1T(x[0]))
		{
			d = -d;
			x[1] = x[0] + d;
			if (ff1T(x[1]) >= ff1T(x[0]))
			{
				p[0] = x[1];
				p[1] = x[0] - d;
				return p;
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
			p[0] = x[i - 1];
			p[1] = x[i + 1];
			return p;
		}

		p[0] = x[i + 1];
		p[1] = x[i - 1];
		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), long double a, long double b, vector<long double> fi, long double epsilon, matrix ud1, matrix ud2) {
	try {
		solution Xopt;
		//ofstream Sout("wyniki_fib.csv");
		int k = 2;
		while (fi[k] <= ((b - a) / epsilon)) {
			k++;
		}

		long double c = b - (fi[k - 2] / fi[k - 1]) * (b - a);
		long double d = a + b - c;
		for (int i = 0; i < k - 2; i++)
		{
			if (ff(c, ud1, ud2) < ff(d, ud1, ud2))
				b = d;

			else
				a = c;

			c = b - (fi[k - i - 3] / (fi[k - i - 2])) * (b - a);
			d = a + b - c;
		}

		Xopt.x = c;
		Xopt.y = Xopt.fit_fun(ff, ud1, ud2);
		Xopt.f_calls = k;
		return Xopt;
	}
	catch (string ex_info) {
		throw ("solution fib(...):\n" + ex_info);
	}
}

//solution fib(matrix(*ff)(matrix, matrix, matrix), long double A, long double B, vector<long double> fi, int Nmax, long double epsilon, matrix ud1, matrix ud2) {
//	try {
//		solution Xopt;
//		ofstream Sout("wyniki_fib.csv");
//		int k = 2;
//		while (fi[k] <= ((B - A) / epsilon)) {
//			k++;
//		}
//
//		vector<long double> a(k), b(k), c(k), d(k);
//
//		a[0] = A;
//		b[0] = B;
//		c[0] = b[0] - ((fi[k - 1]) / fi[k]) * (b[0] - a[0]);
//		d[0] = a[0] + b[0] - c[0];
//		int i = 0;
//		solution A(A), B(B), C(c[0]), D(d[0]);
//
//		for (i; i < k - 2; i++) {
//			A.fit_fun(ff, ud1, ud2);
//			B.fit_fun(ff, ud1, ud2);
//			C.fit_fun(ff, ud1, ud2);
//			D.fit_fun(ff, ud1, ud2);
//
//			if (C.y < D.y) {
//				a[i + 1] = a[i];
//				b[i + 1] = d[i];
//			}
//			else {
//				b[i + 1] = b[i];
//				a[i + 1] = c[i];
//			}
//
//			c[i + 1] = b[i + 1] - ((fi[k - i - 2]) / fi[k - i - 1]) * (b[i + 1] - a[i + 1]);
//			d[i + 1] = a[i + 1] + b[i + 1] - c[i + 1];
//			Sout << hcat(i, c[i]) << "\n";
//
//			A.x = a[i];
//			B.x = b[i];
//			C.x = c[i];
//			D.x = d[i];
//		}
//
//		Sout.close();
//		Xopt.x = c[i];
//		Xopt.y = Xopt.fit_fun(ff, ud1, ud2);
//		return Xopt;
//	}
//	catch (const std::string& ex_info) {
//		throw std::string("solution fib(...):\n" + ex_info);
//	}
//}

solution lag(matrix(*ff)(matrix, matrix, matrix), long double a, long double b, long double epsilon, long double gamma, int Nmax, matrix ud1, matrix ud2) {
	try
	{
		solution Xopt;
		int i = 0;

		long double c = a + (b - a) / 3;
		long double licznik, mianownik;
		long double d = b - (b - a) / 3;
		long double f_a, f_b, f_c, f_d;

		do {
			f_a = ff(a, ud1, ud2)(0, 0);
			f_b = ff(b, ud1, ud2)(0, 0);
			f_c = ff(c, ud1, ud2)(0, 0);

			licznik = f_a * (pow(b, 2)
				- pow(c, 2))
				+ f_b * (pow(c, 2)
					- pow(a, 2)) + f_c * (pow(a, 2)
						- pow(b, 2));
			mianownik = f_a * (b - c)
				+ f_b * (c - a)
				+ f_c * (a - b);

			if (mianownik <= 0) {
				//throw ("Blad - mianownik niedodatni!\n");
				cout << "Blad - mianownik niedodatni!\n";
				return NULL;
			}

			d = 0.5 * (licznik / mianownik);
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
			i++;
			if (i > Nmax)
				break;

		} while (b - a >= epsilon || abs(f_b - f_a) > gamma);

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

//solution lag(matrix(*ff)(matrix, matrix, matrix), long double A, long double B, long double epsilon, long double gamma, int Nmax, matrix ud1, matrix ud2) {
//	try
//	{
//		solution Xopt;
//		vector<long double> a(Nmax), b(Nmax), c(Nmax), d(Nmax);
//		a[0] = A;
//		b[0] = B;
//		c[0] = A + (B - A) / 3.0;
//		long double licznik, mianownik;
//		int i = 0;
//
//		do {
//			long double f_a = ff(a[i], ud1, ud2)(0, 0);
//			long double f_b = ff(b[i], ud1, ud2)(0, 0);
//			long double f_c = ff(c[i], ud1, ud2)(0, 0);
//			
//			licznik = f_a * (pow(b[i], 2))
//						- pow(c[i], 2) 
//						+ f_b * pow(c[i], 2)
//						- pow(a[i], 2) + f_c * (pow(a[i], 2))
//						- pow(b[i], 2);
//			mianownik = f_a * (b[i] - c[i])
//								+ f_b * (c[i] - a[i])
//								+ f_c * (a[i] - b[i]);
//
//			if (mianownik <= 0) {
//				//throw ("Blad - mianownik niedodatni!\n");
//				cout << "Blad - mianownik niedodatni!\n";
//				return NULL;
//			}
//				
//			d[i] = 0.5 * (licznik / mianownik);
//			long double f_d = ff(d[i], ud1, ud2)(0, 0);
//
//			if (a[i] < d[i] && d[i] < c[i]) {
//				if (f_d < f_c) {
//					a[i + 1] = a[i];
//					c[i + 1] = d[i];
//					b[i + 1] = c[i];
//				}
//				else {
//					a[i + 1] = d[i];
//					c[i + 1] = c[i];
//					b[i + 1] = b[i];
//				}
//			}
//
//			else {
//				if (c[i] < d[i] && d[i] < b[i]) {
//					if (f_d < f_c) {
//						a[i + 1] = c[i];
//						c[i + 1] = d[i];
//						b[i + 1] = b[i];
//					}
//					else {
//						a[i + 1] = a[i];
//						c[i + 1] = c[i];
//						b[i + 1] = d[i];
//					}
//				}
//				/*else {
//					cout << "Blad - algorytm nie jest zbiezny!\n";
//					return NULL;
//				}*/
//					//throw ("Blad - algorytm nie jest zbiezny!\n");
//					//break;
//			}
//			cout << i << endl;
//			i++;
//
//			/*if (i > Nmax) {
//				cout << "Blad - nie udalo sie osiagnac doklanosci epsilon\n";
//				return NULL;
//			}*/
//				//throw ("Blad - nie udalo sie osiagnac doklanosci epsilon\n");
//				//break;
//		} while (b[i] - a[i] >= epsilon || abs(d[i] - d[i - 1]) > gamma);
//
//		Xopt.x = d[i];
//		return Xopt;
//	}
//
//	catch (string ex_info)
//	{
//		throw ("solution lag(...):\n" + ex_info);
//	}
//
//}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
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
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
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
