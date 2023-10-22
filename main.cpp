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
		lab0();
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
	double epsilon = 1e-2;				//dok³adnoœæ obliczeñ
	int Nmax = 10000;					//max liczba iteracji
	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);	//klasa matrix i wywo³anie konstruktora domyœlny z rozmiarem macierzy, lb-> dolen ograniczenie ub->górne ogranicznie 2 wiersze 1 kolumna wartsoc 5 i punkt a
	solution opt;	//obiekt w ktorym przechowujemy roziwaznie , przechowuje liczbe wywo³añ funkjci celu i rozwiazanie 
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a); //algorytm optymalizacji argumenty: funkcja celu, 2-> l. wymairów, lb- lu-> ograniczenia na losowany punkt  eps-> dok³adnoœc, N-> l. iteracji i a -> punkt 
	cout << opt << endl << endl;
	solution::clear_calls();

	//Wahadlo
	// wahadlo o dl l i masie m na ppocztaku jest nie ruchome kat wychylenia teta w radianach i na po³ sekundy dodajemy 
	// sily nie przekroczyæ tata_opt wyliczyæ max monemt by nie wychyl sie za bardzo 
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
