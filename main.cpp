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
		lab1();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	system("pause");
	return 0;
}


void lab1()
{
	long double epsilon = 1e-5;
	long double gamma = 1e-7;
	long double x0 = 50;
	long double d = 5.0;
	long double alpha = 2.0;
	long double* p = new long double[2]{};
	int Nmax = 100;
	matrix ud1, ud2;

	p = expansion(ff1T, x0, d, alpha, Nmax, ud1, ud2);
	long double a = p[0];
	long double b = p[1];

	cout << "===== FUNKCJA TESTOWA =====\n";
	cout << "Rzeczywiste minimum funkcji: x = 62.7481, y = -0.9211\n";
	cout << "===========================\n\n";

	cout << "=== METODA EKSPANSJI ======\n";
	cout << "Parametry poczatkowe: x0 = " << x0 << ", d = " << d << endl;
	cout << "Uzyskany przedzial: [" << a << ", " << b << "]\n\n";


	cout << "=== METODA FIBONACCIEGO ===\n";
	vector<long double> fi;
	fi.push_back(1);
	fi.push_back(1);
	for (int i = 2; i < 50; i++) {
		fi.push_back(fi[i - 2] + fi[i - 1]);
	}

	solution optT;
	optT = fib(ff1T, a, b, fi, epsilon);
	cout << optT << endl;
	solution::clear_calls();


	cout << "=== METODA LAGRANGE'A =====\n";
	optT = lag(ff1T, a, b, epsilon, gamma, Nmax);
	if (optT.x != NULL) cout << optT << endl;
	solution::clear_calls();

	cout << "===== PROBLEM RZECZYWISTY =====\n";
	a = 0;
	b = 100;
	cout << "===============================\n\n";
	cout << "=== METODA FIBONACCIEGO =======\n";
	solution optR;
	optR = fib(ff1R, a, b, fi, epsilon);
	cout << optR << endl;
	solution::clear_calls();

	cout << "=== METODA LAGRANGE'A =====\n";
	optR = lag(ff1R, a, b, epsilon, gamma, Nmax);
	cout << optR << endl;
	solution::clear_calls();

}

void lab2()
{

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
