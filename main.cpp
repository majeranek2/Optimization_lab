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
	ofstream Sout("wyniki.csv");
	 double epsilon = 1e-5;
	 double gamma = 1e-9;
	 double x0= 95;
	 double d = 3.0;
	 double alpha = 1.6;
	 //srand(time(NULL));
	 solution optT;
	 solution optR;
	 solution optE;
	 double* p = new  double[2]{};
	 double pom;
	 int Nmax = 100;
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

	/*cout << "===== PROBLEM RZECZYWISTY =====\n";
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
	solution::clear_calls();*/

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
