/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include"opt_alg.h"

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		lab2();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	system("pause");
	return 0;
}

/*
void lab1()
{	
	ofstream Sout("wyniki.csv");
	 double epsilon = 1e-5;
	 double gamma = 1e-9;
	 double x0= 95;
	 double d = 3.0;
	 double alpha = 1.6;
	 solution optT;
	 solution optRfib;
	 solution optRlag;
	 solution optE;
	 double* p = new  double[2]{};
	 double pom;
	 int Nmax = 1000;
	 double a, b;
	 matrix ud1, ud2;


	optE = expansion(ff1T, x0, d, alpha, Nmax, ud1, ud2);
	a = optE.x(0,0);
	b = optE.y(0,0);
	 
	 cout << "===== FUNKCJA TESTOWA =====\n";
	cout << "Rzeczywiste minimum funkcji: x = 62.7481, y = -0.9211\n";
	cout << "===========================\n\n";

	cout << "=== METODA EKSPANSJI ======\n";
	cout << "Parametry poczatkowe: x0 = " << x0 << ", d = " << d << endl;
	cout << "Uzyskany przedzial: [" << a << ", " << b << "]\n\n";


	


	cout << "=== METODA FIBONACCIEGO ===\n";
	vector< double> fi;
	fi.push_back(1);
	fi.push_back(1);
	for (int i = 2; i < 50; i++) {
		fi.push_back(fi[i - 2] + fi[i - 1]);
	}
	

	optT = fib(ff1T, a, b, fi, epsilon);
	cout << optT << endl;
	solution::clear_calls();


	cout << "=== METODA LAGRANGE'A =====\n";
	optT = lag(ff1T, a, b, epsilon, gamma, Nmax);
	if (optT.x != NULL) cout << optT << endl;
	solution::clear_calls();
	
	
	//ZADANIE 1
	for (alpha = 1.2; alpha < 2.5; alpha += 0.3) {
		
		for (int i = 0; i < 100; i++) {
			a = 0;
			b = 100;
			x0 = x0-0.2;
			optE = expansion(ff1T, x0, d, alpha, Nmax, ud1, ud2);
			a = optE.x(0, 0);
			b = optE.y(0, 0);
			Sout << hcat(a, b);
			Sout << hcat(x0, optE.f_calls);

			optT = fib(ff1T, a, b, fi, epsilon);
			Sout << hcat(optT.x, optT.y);
			Sout << hcat(optT.f_calls, optT.f_calls);
			solution::clear_calls();

			optT = lag(ff1T, a, b, epsilon, gamma, Nmax);
			Sout << hcat(optT.x, optT.y);
			Sout << optT.f_calls << "\n";
			solution::clear_calls();
		}
		Sout << "\n\n";
	}Sout.close();

	cout << "===== PROBLEM RZECZYWISTY =====\n";
	a = 1e-4;
	b = 1e-1;
	cout << "===============================\n\n";
	cout << "=== METODA FIBONACCIEGO =======\n";
	optRfib = fib(ff1R, a, b, fi, epsilon);
	cout << optRfib << endl;
	solution::clear_calls();

	cout << "=== METODA LAGRANGE'A =====\n";
	optRlag = lag(ff1R, a, b, epsilon, gamma, Nmax);
	cout << optRlag << endl;
	solution::clear_calls();


	// SYMULACJA - FIBONACCI
	ofstream FIB("symulacja_fibonacci.csv");
	matrix x(optRfib.x);
	double Va = 5;
	double Vb = 1;
	double Tb = 10;
	double tend = 1000;
	double t0 = 0;
	double dt = 1;
	matrix y0 = matrix(3, new double[3] { Va, Vb, Tb });
	matrix* y = solve_ode(df1, t0, dt, tend, y0, ud1, x);
	double max = y[1](0, 2);
	for (int t = 0; t < 1000; t++)
	{
		FIB << hcat(t, y[1](t, 0)) << hcat(y[1](t, 1), y[1](t, 2)) << "\n";
	}
	FIB.close();

	// SYMULACJA - LAGRANGE
	ofstream LAG("symulacja_lagrange.csv");
	x = matrix(optRlag.x);
	y = solve_ode(df1, 0, 1, 1000, y0, ud1, x);
	max = y[1](0, 2);
	for (int t = 0; t < 1000; t++)
	{
		LAG << hcat(t, y[1](t, 0)) << hcat(y[1](t, 1), y[1](t, 2)) << "\n";
	}
	LAG.close();
	solution::clear_calls();
	y[0].~matrix();
	y[1].~matrix();

	// SYMULACJA - LAGRANGE
	ofstream LAG("symulacja_lagrange.csv");
	x = matrix(optRlag.x);
	y = solve_ode(df1, 0, 1, 1000, y0, ud1, x);
	max = y[1](0, 2);
	for (int t = 0; t < 1000; t++)
	{
		LAG << hcat(t, y[1](t, 0)) << hcat(y[1](t, 1), y[1](t, 2)) << "\n";
	}
	LAG.close();
	solution::clear_calls();
}
*/


void lab2()
{
	srand(time(NULL));
	matrix ud1, ud2;
	solution optT;
	double epsilon = 0.01;
	double alpha = 0.5;
	int Nmax = 100;
	matrix x0(2, 1);
	matrix x(2, 1);
	x0(0, 0) = 0.5;  // Wartoœæ pocz¹tkowa dla x1
	x0(1, 0) = 0.5;
	double s = 0.5;
	//optT = HJ(ff2T, x0, s, alpha, epsilon,Nmax,ud1,ud2);
	//cout << optT << endl;
	//solution::clear_calls();

	//z wykladu dla sprawdzenia:---------------------------------------------
	/*
	x0(0, 0) = -0.5;  // Wartoœæ pocz¹tkowa dla x1
	x0(1, 0) = 1;
	s = 0.5;
	epsilon = 0.01;
	alpha = 0.5;
	optT = HJ(ff2T_2, x0, s, alpha, epsilon, Nmax, ud1, ud2);
	cout << optT << endl;
	solution::clear_calls();
	*/

	//Funkcja testowa celu zadanie:
	ofstream Sout("wyniki_lab_2.csv");
	for (double krok = 0.4; krok < 1.5; krok += 0.5) {
		for (int i = 0; i < 100; i++) {
			x0(0, 0) = rand() % 10;
			x0(1, 0) = rand() % 10;
			optT = HJ(ff2T, x0, krok, alpha, epsilon, Nmax, ud1, ud2);
			x = optT.x;
			Sout << hcat(x0(0, 0), x0(1, 0));
			Sout << hcat(x(0, 0), x(1, 0));
			Sout << optT.f_calls << "\n";
			//Zamieniæ wy¿sza linija na to przy dodaniu rosenbroka, bo /n rozpierdzieli 
			// b³edy bedzue calls dwa razy wypisany jak ostatnio
			//Sout << hcat(optT.f_calls, optT.f_calls);
			solution::clear_calls();


			//algorytm Rosenbroka
			/*optT = Rosen(ff2T, x0, krok, alpha, epsilon, Nmax, ud1, ud2);
			x = optT.x;
			Sout << hcat(x0(0, 0), x0(1, 0));
			Sout << hcat(x(0, 0), x(1, 0));
			Sout << optT.f_calls << "\n";
			solution::clear_calls();*/

		}
		Sout << "\n\n";
	}
	Sout.close();


	matrix X0(2, 2);
	// losowe liczby double od 0 do 1
	X0(0, 0) = -1.0 + 2.0 * static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
	X0(0, 1) = -1.0 + 2.0 * static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
	X0(1, 0) = -1.0 + 2.0 * static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
	X0(1, 1) = -1.0 + 2.0 * static_cast <double> (rand()) / static_cast <double> (RAND_MAX);

	alpha = 3.0;
	double beta = 0.5;
	int n = get_len(X0[0]);
	matrix S(n, n, 1.0);
	cout << "Punkty poczatkowe:\n";
	cout << "x0: " << endl << X0 << endl;

	optT = Rosen(ff2T, X0, S, alpha, beta, epsilon, Nmax);
	cout << optT << endl;
	solution::clear_calls();

	//Problem rzeczywtsty:------------------------------------
	x0(0, 0) = 5;
	x0(1, 0) = 5;
	s = 1;

	//optT = HJ(ff2R, x0, s, alpha, epsilon, Nmax, ud1, ud2);
	//cout << optT << endl;
	//solution::clear_calls();
}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
