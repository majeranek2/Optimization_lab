#include"opt_alg.h"
double const PHI = 1.61803398;

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
			Xopt.fit_fun(ff, ud1, ud2); //wywolanie f. celu nna rzecz rozwiazania fitcals+1
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1; //dowolne jak chcemy tutaj 1-> dziala mamy rozwiazanie i finish
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;	//-> nie znaleziono za duzo iteacji 
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

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2] { 0, 0 };
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
				exit(1);
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

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, vector<int> fi, int Nmax, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		int k = 2;
		while (fi[k] <= (b - a) / epsilon) {
			k++;
		}
		k++;
		cout << k;
		vector <double> a_, b_, c_, d_;
		a_.push_back(a);
		b_.push_back(b);
		c_.push_back(b_[0] - (fi[k - 1] /(fi[k]) * (b_[0] - a_[0])));
		d_.push_back(a_[0] + b_[0] - c_[0]);

		int i = 0;
		for (i; i <= k-4; i++) {
			if (ff1T(c_[i]) < ff1T(d_[i])) {
				a_.push_back(a_[i]);
				b_.push_back(d_[i]);
			}
			else {
				b_.push_back(b_[i]);
				a_.push_back(c_[i]);
			}
			c_.push_back(b_[i + 1] - (fi[k - i - 2] / fi[k - i - 1]) * (b_[i + 1] - a_[i + 1]));
			d_.push_back(a_[i + 1] + b_[i + 1] - c_[i + 1]);
		}
		solution Xopt;
		Xopt.x = c_[k-3];
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}




solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2) {
	try {

		int i = 0;
		double a_k = a, b_k = b, c_k = a_k + (b_k - a_k) / 3.0;

		while (true) {
			double l = ff(matrix(a_k), ud1, ud2)(0, 0) * ((b_k * b_k) - (c_k * c_k)) +
				ff(matrix(b_k), ud1, ud2)(0, 0) * ((c_k * c_k) - (a_k * a_k)) +
				ff(matrix(c_k), ud1, ud2)(0, 0) * ((a_k * a_k) - (b_k * b_k));

			double m = ff(matrix(a_k), ud1, ud2)(0, 0) * (b_k - c_k) +
				ff(matrix(b_k), ud1, ud2)(0, 0) * (c_k - a_k) +
				ff(matrix(c_k), ud1, ud2)(0, 0) * (a_k - b_k);

			if (m <= 0) {
				std::cout << "Nie dzia³a" << std::endl;
				return solution();
			}

			double d_k = 0.5 * l / m;

			if (a_k < d_k && d_k < c_k) {
				if (ff(matrix(d_k), ud1, ud2)(0, 0) < ff(matrix(c_k), ud1, ud2)(0, 0)) {
					b_k = c_k;
					c_k = d_k;
				}
				else {
					a_k = d_k;
				}
			}
			else if (c_k < d_k && d_k < b_k) {
				if (ff(matrix(d_k), ud1, ud2)(0, 0) < ff(matrix(c_k), ud1, ud2)(0, 0)) {
					a_k = c_k;
					c_k = d_k;
				}
				else {
					b_k = d_k;
				}
			}
			else {
				throw "Invalid operation: d_k is out of range";
			}

			i++;
			solution::f_calls++;

			if (i > Nmax) {
				throw "Exceeded maximum number of function calls.";
			}

			if (b_k - a_k < epsilon || abs(d_k - 1) < gamma) {
				solution Xopt;
				Xopt.x = d_k;
				Xopt.y = ff(matrix(Xopt.x), ud1, ud2)(0, 0);

				return Xopt;
			}
		}
	}
	catch (const char* ex_info) {
		throw "solution lag(...):\n" + std::string(ex_info);
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
