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
