/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
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

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;				//dok�adno�� oblicze�
	int Nmax = 10000;					//max liczba iteracji
	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);	//klasa matrix i wywo�anie konstruktora domy�lny z rozmiarem macierzy, lb-> dolen ograniczenie ub->g�rne ogranicznie 2 wiersze 1 kolumna wartsoc 5 i punkt a
	solution opt;	//obiekt w ktorym przechowujemy roziwaznie , przechowuje liczbe wywo�a� funkjci celu i rozwiazanie 
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a); //algorytm optymalizacji argumenty: funkcja celu, 2-> l. wymair�w, lb- lu-> ograniczenia na losowany punkt  eps-> dok�adno�c, N-> l. iteracji i a -> punkt 
	cout << opt << endl << endl;
	solution::clear_calls();

	//Wahadlo
	// wahadlo o dl l i masie m na ppocztaku jest nie ruchome kat wychylenia teta w radianach i na po� sekundy dodajemy 
	// sily nie przekroczy� tata_opt wyliczy� max monemt by nie wychyl sie za bardzo 
	//mc^2O(..) +mglsinO+bO(.)-M=0
	//szukamy x
	//O(t) i dO(t)/dt
	Nmax = 1000;	
	epsilon = 1e-2;
	lb = 0;
	ub = 5;
	double teta_opt = 1;	//
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(opt.x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
	ofstream Sout("symulacja_lab0.csv");
	Sout << hcat(Y[0], Y[1]);
	Sout.close();
	Y[0].~matrix();
	Y[1].~matrix();
}


void lab1()
{
	matrix ud1, ud2;
	int Nmax = 100;
	double x0 = 52.3;
	double d = 10;
	double alpha = 1.5;
	double* p = new double[2] { -100,100 };
	p = expansion(ff1T, x0, d, alpha, Nmax , ud1, ud2);
	cout << "(" << p[0] << ", " << p[1] << ")" << endl << endl;

	double a = p[0];
	double b = p[1];
	double epsilon = 0.0001;
	double gamma = 0.00001;

	vector<int> fi;
	fi.push_back(1);
	fi.push_back(1);
	for (int i = 2; i < 100; i++) {
		fi.push_back(fi[i - 2] + fi[i - 1]);
	}
	
	solution opt1;
	opt1 = fib(ff1T, a, b, fi, Nmax, epsilon, ud1, ud2);
	cout << "Optimal point FIB: "  <<opt1.x << endl;
	solution::clear_calls();

	solution opt2;
	opt2= lag(ff1T, a, b, epsilon, gamma, Nmax, ud1, ud2);
	cout << "Optimal point: "  << opt2.x << endl;
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
