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

//solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, vector<int> fi, int Nmax, double epsilon, matrix ud1, matrix ud2)
//{
//	try
//	{
//		
//		ofstream Sout("wyniki_fib.csv");
//		int k = 2;
//		while (fi[k] <= static_cast<int>((b - a) / epsilon))
//		{
//			k++;
//		}
//
//		std::vector<double> a_(k);
//		std::vector<double> b_(k);
//		std::vector<double> c_(k);
//		std::vector<double> d_(k);
//
//		a_[0] = a;
//		b_[0] = b;
//		c_[0] = b_[0] - (static_cast<double>(fi[k - 1]) / fi[k]) * (b_[0] - a_[0]);
//		d_[0] = a_[0] + b_[0] - c_[0];
//		int i = 0;
//		for (i; i < k - 2; i++)
//		{
//			if (ff1R(c_[i], ud1, ud2) < ff1R(d_[i], ud1, ud2))
//			{
//				a_[i + 1] = a_[i];
//				b_[i + 1] = d_[i];
//			}
//			else
//			{
//				b_[i + 1] = b_[i];
//				a_[i + 1] = c_[i];
//			}
//
//			c_[i + 1] = b_[i + 1] - (static_cast<double>(fi[k - i - 2]) / fi[k - i - 1]) * (b_[i + 1] - a_[i + 1]);
//			d_[i + 1] = a_[i + 1] + b_[i + 1] - c_[i + 1];
//			Sout << hcat(i, c_[i]) << "\n";
//			
//			}
//		
//		Sout.close();
//		solution Xopt;
//		Xopt.x = c_[i];
//		
//
//		Xopt.y = Xopt.fit_fun(ff, ud1, ud2);
//		return Xopt;
//	}
//	catch (const std::string& ex_info)
//	{
//		throw std::string("solution fib(...):\n" + ex_info);
//	}
//
//}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, vector<int> fi, int Nmax, double epsilon, matrix ud1, matrix ud2) {
	try {
		solution Xopt;
		ofstream Sout("wyniki_fib.csv");
		int k = 2;
		while (fi[k] <= static_cast<int>((b - a) / epsilon)) {
			k++;
		}

		std::vector<double> a_(k);
		std::vector<double> b_(k);
		std::vector<double> c_(k);
		std::vector<double> d_(k);

		a_[0] = a;
		b_[0] = b;
		c_[0] = b_[0] - (static_cast<double>(fi[k - 1]) / fi[k]) * (b_[0] - a_[0]);
		d_[0] = a_[0] + b_[0] - c_[0];
		int i = 0;
		solution A(a), B(b), C(c_[0]), D(d_[0]);

		for (i; i < k - 2; i++) {
			A.fit_fun(ff, ud1, ud2);
			B.fit_fun(ff, ud1, ud2);
			C.fit_fun(ff, ud1, ud2);
			D.fit_fun(ff, ud1, ud2);

			if (C.y < D.y) {
				a_[i + 1] = a_[i];
				b_[i + 1] = d_[i];
			}
			else {
				b_[i + 1] = b_[i];
				a_[i + 1] = c_[i];
			}

			c_[i + 1] = b_[i + 1] - (static_cast<double>(fi[k - i - 2]) / fi[k - i - 1]) * (b_[i + 1] - a_[i + 1]);
			d_[i + 1] = a_[i + 1] + b_[i + 1] - c_[i + 1];
			Sout << hcat(i, c_[i]) << "\n";

			A.x = a_[i];
			B.x = b_[i];
			C.x = c_[i];
			D.x = d_[i];
		}

		Sout.close();
		Xopt.x = c_[i];
		Xopt.y = Xopt.fit_fun(ff, ud1, ud2);
		return Xopt;
	}
	catch (const std::string& ex_info) {
		throw std::string("solution fib(...):\n" + ex_info);
	}
}


solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{	
		ofstream Sout("wyniki_lag.csv");
		solution Xopt;
		Xopt.ud = b - a;

		solution A(a), B(b), C, D, D_old(a);
		C.x = (a + b) / 2;

		A.fit_fun(ff, ud1, ud2);
		B.fit_fun(ff, ud1, ud2);
		C.fit_fun(ff, ud1, ud2);

		double l, m;
		int i = 0;
		while (true)
		{
			l = m2d(A.y * (pow(B.x) - pow(C.x)) + B.y * (pow(C.x) - pow(A.x)) + C.y * (pow(A.x) - pow(B.x)));
			m = m2d(A.y * (B.x - C.x) + B.y * (C.x - A.x) + C.y * (A.x - B.x));

			//cout << abs(m2d(B.x)) << " " << abs(m2d(A.x)) << "\n";
			//cout << abs(m2d(B.x) - m2d(A.x)) << "\n";

			if (m <= 0)
			{
				Xopt = D_old;
				Xopt.flag = 2;
				return Xopt;
			}

			D.x = 0.5 * l / m;
			D.fit_fun(ff, ud1, ud2);

			if (A.x <= D.x && D.x <= C.x)
			{
				if (D.y < C.y)
				{
					B = C;
					C = D;
				}
				else
				{
					A = D;
				}
			}
			else if (C.x <= D.x && D.x <= B.x)
			{
				if (D.y < C.y)
				{
					A = C;
					C = D;
				}
				else
				{
					B = D;
				}
			}
			else
			{
				Xopt = D_old;
				Xopt.flag = 2;

				//cout << abs(m2d(B.x)) << " " << abs(m2d(A.x)) << "\n";
				//cout << abs(m2d(B.x) - m2d(A.x)) << "\n";
				//cout << m2d(D.x) << "\n\n\n";

				return Xopt;
			}

			Xopt.ud.add_row((B.x - A.x)());

			if (B.x - A.x < epsilon || abs(D.x() - D_old.x()) < gamma)
			{
				Xopt = D;
				Xopt.flag = 0;

				//cout << abs(m2d(B.x)) << " " << abs(m2d(A.x)) << "\n";
				//cout << abs(m2d(B.x) - m2d(A.x)) << "\n";
				//cout << m2d(D.x) << "\n\n\n";

				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt = D;
				Xopt.flag = 1;
				break;
			}
			
			D_old = D;
			Sout << D.x() << "\n";
			i++;

			//cout << abs(m2d(B.x)) << " " << abs(m2d(A.x)) << "\n";
			//cout << abs(m2d(B.x) - m2d(A.x)) << "\n";
			//cout << m2d(D.x) << "\n\n\n";
		}
		//cout << abs(m2d(B.x)) << " " << abs(m2d(A.x)) << "\n";
		//cout << abs(m2d(B.x) - m2d(A.x)) << "\n";
		//cout << m2d(D.x) << "\n\n\n";
		Sout.close();

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}
//solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2) {
//		try {
//			std::ofstream Sout("wyniki_lag.csv");
//			int i = 0;
//			std::vector<double> a_k = { a };
//			std::vector<double> b_k = { b };
//			std::vector<double> c_k = { a + (b - a) / 3.0 };
//			double d_k = 0.0;
//			double previous_d_k = 0.0;
//
//			do {
//				double l = ff(matrix(a_k.back()), ud1, ud2)(0, 0) * (b_k.back() * b_k.back() - c_k.back() * c_k.back()) +
//					ff(matrix(b_k.back()), ud1, ud2)(0, 0) * (c_k.back() * c_k.back() - a_k.back() * a_k.back()) +
//					ff(matrix(c_k.back()), ud1, ud2)(0, 0) * (a_k.back() * a_k.back() - b_k.back() * b_k.back());
//
//				double m = ff(matrix(a_k.back()), ud1, ud2)(0, 0) * (b_k.back() - c_k.back()) +
//					ff(matrix(b_k.back()), ud1, ud2)(0, 0) * (c_k.back() - a_k.back()) +
//					ff(matrix(c_k.back()), ud1, ud2)(0, 0) * (a_k.back() - b_k.back());
//
//
//				if (m <= 0) {
//					std::cout << "Nie dzia³a" << std::endl;
//					return solution();
//				}
//				else {
//					d_k = 0.5 * l / m;
//
//					if (a_k.back() < d_k && d_k < c_k.back()) {
//						if (ff(matrix(d_k), ud1, ud2)(0, 0) < ff(matrix(c_k.back()), ud1, ud2)(0, 0)) {
//							b_k.push_back(c_k.back());
//							c_k.push_back(d_k);
//						}
//						else {
//							a_k.push_back(d_k);
//						}
//					}
//					else if (c_k.back() < d_k && d_k < b_k.back()) {
//						if (ff(matrix(d_k), ud1, ud2)(0, 0) < ff(matrix(c_k.back()), ud1, ud2)(0, 0)) {
//							a_k.push_back(c_k.back());
//							c_k.push_back(d_k);
//						}
//						else {
//							b_k.push_back(d_k);
//						}
//					}
//					else {
//						throw "Invalid operation: d_k is out of range";
//					}
//				}
//
//				Sout << hcat(i, d_k) << "\n";
//				i++;
//				solution::f_calls++;
//
//				if (i > Nmax) {
//					throw "Exceeded the maximum number of function calls.";
//				}
//
//				previous_d_k = d_k;
//
//			} while (std::abs(b_k.back() - a_k.back()) >= epsilon || std::abs(d_k - previous_d_k) > gamma);
//
//			Sout.close();
//
//			solution Xopt;
//			Xopt.x = d_k;
//			Xopt.y = ff(matrix(Xopt.x), ud1, ud2)(0, 0);
//			return Xopt;
//
//		}
//		catch (const char* ex_info) {
//			throw std::string("solution lag(...):\n") + std::string(ex_info);
//		}
//	}






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
