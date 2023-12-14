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
		
		ofstream Sout("HJ.csv");
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
				Sout << hcat(x(0,0), x(1,0));
				Sout << "\n";
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
		Sout.close();
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



solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, matrix alpha, matrix beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		matrix x_b = x0;	
		matrix s = s0;
		int n = get_len(x_b);

		double* ZERO = new double[n] {};
		matrix lambda(n, ZERO);
		matrix p(n, ZERO);
		matrix d = getBase(n);
		int i = 0;
		ofstream Sout("Rosen.csv");
		while (i <= Nmax) {
			// DEBUG
			/*cout << "\nNUMER ITERACJI: " << i << endl;
			cout << "Polozenie punktu x" << i << ": " << endl << x_b << endl;
			cout << "Macierz d: " <<endl << d << endl;
			cout << "Macierz s: " <<endl << s << endl;
			cout << "-----------------------\n";*/

			for (int j = 0; j < n; j++) {
				//cout << "Kierunek: " << j << endl;
				if (ff(x_b + (get_row(s, j) * d[j]), ud1, ud2) < ff(x_b, ud1, ud2)) {
					//cout << "WARUNEK SPE£NIONY: " << ff(x_b + (get_row(s, j) * d[j]), ud1, ud2) << " < " << ff(x_b, ud1, ud2) << endl;
					Sout << hcat(x_b(0, 0), x_b(1, 0));
					Sout << "\n";
					x_b = x_b + (get_row(s, j) * d[j]);
					add_scalar_to_matrix_row(lambda, get_row(s, j), j);
					multiply_matrix_row_by_scalar(s, alpha, j);
					
				}
				else {
					//cout << "WARUNEK NIE JEST SPE£NIONY: " << ff(x_b + (get_row(s, j) * d[j]), ud1, ud2) << " > " << ff(x_b, ud1, ud2) << endl;
					multiply_matrix_row_by_scalar(s, -beta, j);
					matrix one(1.0);
					add_scalar_to_matrix_row(p, one, j);
				}
				// DEBUG
				/*cout << endl << "Polozenie punktu x" << i << ": " << endl << x_b << endl << endl;
				cout << "Lambda: " << endl << lambda << endl << endl;
				cout << "p: " << endl << p << endl << endl;
				cout << "s: " << endl << s << endl << endl;*/
			}
			
			matrix zero(0.0);
			bool zmiana_bazy = false;
			for (int j = 0; j < n; j++) {
				if (get_row(lambda, j) != zero && get_row(p, j) != zero) {
					zmiana_bazy = true;
				}
				else {
					zmiana_bazy = false;
					break;
				}
			}

			if (zmiana_bazy) {
				matrix Q = Rosen_getQ(n, lambda, d);
				//cout << "\n\n\n TUTAJ MACIERZ Q!!\n" << Q << endl;

				for (int j = 0; j < n; j++) {
					/*cout << "SANITY CHECK DLA J = " << j << endl;
					cout << "MACIERZ Q:\n" << Q << endl;
					cout << "ROSEN GET VJ:\n" << Rosen_getVj(n, j, Q, d) << endl;
					cout << "NORM:\n" << norm(Rosen_getVj(n, j, Q, d)) << endl;
					cout << "d" << j << "\n" << (Rosen_getVj(n, j, Q, d) / norm(Rosen_getVj(n, j, Q, d))) << endl;*/
					d.set_col((Rosen_getVj(n, j, Q, d) / norm(Rosen_getVj(n, j, Q, d))), j);
				}

				for (int k = 0; k < n; k++) {
					lambda.set_row(zero, k);	// zerujemy macierze
					p.set_row(zero, k);	
				}
				s = s0;
			}

			// DEBUG
			if (zmiana_bazy) {
				//cout << "ZMIANA BAZY.\n" << "Nowa baza:\n" << d;
			}
			//cout << "===============================\n";
			//cout << "===============================\n\n";
		
			bool ifBreak = false;
			i++;
			for (int j = 0; j < n; j++) {
				double a = m2d(get_row(s, j));
				if (abs(a) < epsilon) {
					cout << "Uzyskano zbieznosc po " << i << " iteracjach!\n\n";
					ifBreak = true;
					break;
				}
			}
			if (ifBreak == true)
				break;
		}
		// tutaj: warunek na epsilon
		Sout << hcat(x_b(0, 0), x_b(1, 0));
		Sout << "\n";
		Xopt.x = x_b;
		Xopt.y = ff(Xopt.x, ud1, ud2)(0, 0);
		Xopt.f_calls = i;
		Sout.close();
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

matrix add_scalar_to_matrix_row(matrix &m, matrix scalar, int row) {
	matrix element = get_row(m, row);
	element = element + scalar;
	m.set_row(element, row);
	return m;
}

matrix multiply_matrix_row_by_scalar(matrix& m, matrix scalar, int row) {
	matrix element = get_row(m, row);
	element = scalar * element;
	m.set_row(element, row);
	return m;
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
				tab[i][j] = m2d(get_row(lambda, i));
			}
		}
	}

	matrix Q_r(n, n, tab);
	//cout << "TWORZE MACIERZ Q:\nd:\n" << d << "\nQ_r:\n" << Q_r << endl;
	return d * Q_r;
}

matrix Rosen_getVj(int n, int j, matrix Q, matrix d) {
	if (j == 0) {
		matrix result = Q[0];
		return result;
	}
	else {
		double* ZERO = new double[n] {};
		matrix suma(n, ZERO);
		for (int k = 0; k <= j-1; k++) {
			suma = suma + ((trans(Q[j]) * d[k]) * d[k]);
			//cout << "k = " << k << ", suma matrix = \n" << suma << endl;
		}
		matrix result = Q[j] - suma;
		/*cout << "TUTAJ BUG???\n\n";
		cout << "Q[j] = \n" << Q[j] << endl;
		cout << "suma matrix = \n" << suma << endl;
		cout << "result matrix = \n" << result << endl;*/

		return result;
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
		int n = get_size(x0)[0];	// n - liczba wymiarow
		cout << "n = " << n << endl;
		matrix p = x0;
		double* ZERO = new double[n+1] {};
		matrix values(n + 1, ZERO);
		int pmin;	// do oznaczenia, który wierzcho³ek jest pmin
		int pmax;
		bool condition = true;	// czy przejsc do nastepnej iteracji?
		int calls = 0;

		while (condition) {
			//condition = false;

			// Obliczenie wartoœci w ka¿dym z wierzcho³ków simpleksu
			for (int i = 0; i < n + 1; i++) {
				values.set_row(ff(p[i], ud1, ud2), i);
			}

			// Ustalenie, które wierzcho³ki to pmax, a które to pmin
			pmax = get_max(values, n);
			pmin = get_min(values, n);

			// Wyznaczanie œrodka ciê¿koœci simpleksu
			matrix sum(n, ZERO);
			for (int i = 0; i < n + 1; i++) {
				if (i != pmax)
					sum = sum + p[i];
			}
			matrix p_srodek = sum / n;
			
			// Proba odbicia
			matrix p_odb = p_srodek + alpha * (p_srodek - p[pmax]);
			cout << "P ODBITE:\n" << p_odb << endl;

			if (ff(p_odb, ud1, ud2) < ff(p[pmin], ud1, ud2)) {
				matrix p_e = p_srodek + gamma * (p_odb - p_srodek);
				if (ff(p_e, ud1, ud2) < ff(p_odb, ud1, ud2))	// ekspansja
					p[pmax] = p_e;
				else	// odbicie
					p[pmax] = p_odb;
			}

			else{
				matrix p_z = p_srodek + beta * (p[pmax] - p_srodek);
				if (ff(p_z, ud1, ud2) >= ff(p[pmax], ud1, ud2)) {
					for (int i = 0; i < n + 1; i++) {
						if (i != pmin)
							p[i] = delta * (p[i] + p[pmin]);	// redukcja
					}
				}
				else
					p[pmax] = p_z;
			}

			// warunek stopu - czy poprawny? co z tym max?
			/*for (int i = 0;  i < n + 1; i++) {
				if (norm(pmin - p[i]) >= epsilon) {
					condition = true;
					break;
				}
			}*/
			calls++;
			if (calls > Nmax)
				break;
		}
		Xopt.x = p[pmin];
		Xopt.y = ff(Xopt.x, ud1, ud2)(0, 0);
		Xopt.f_calls = calls;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

matrix get_simplex(int n) {
	matrix simplex(n, n + 1, 0.0);
	for (int i = 0; i < n; i++) {
		simplex(i, i + 1) = 1.0;
	}
	return simplex;
}

int get_max(matrix p, int n) {
	double val = p(0, 0);
	matrix pmax(val);
	int maxindex = 1;

	for (int i = 1; i < n + 1; i++) {
		if (get_row(p, i) > pmax) {
			pmax = get_row(p, i);
			maxindex = i;
		}
	}
	return maxindex;
}

int get_min(matrix p, int n) {
	double val = p(0, 0);
	matrix pmin(val);
	int minindex = 1;

	for (int i = 1; i < n + 1; i++) {
		if (get_row(p, i) < pmin) {
			pmin = get_row(p, i);
			minindex = i;
		}
	}
	return minindex;
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
