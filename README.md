# Optymalizacja_lab
Sprawozdanie lab 2:
[https://aghedupl-my.sharepoint.com/:w:/g/personal/pawinska_student_agh_edu_pl/EUOUpvdyD0lDtlc9eS0vvqoBmTz0mMja4GqCCBLX56xw3A?e=xqubJ2](https://aghedupl-my.sharepoint.com/:w:/g/personal/pawinska_student_agh_edu_pl/Ed0YmHMnY1VJkq3HAJoUPJkBbsQhd8lRynAbp8p-Bor-Dg?e=PMGSIj)

funkcje piszemy w pliku opt_alg.cpp
wywołujwemy je w main lab_1
aktualizacje git hub!
1. ściagnij zmiany przed edycja bo zabijamy za bałagfan :)
2. zmainy wproawdzamy
3. WIDOK-> ZMIANY GIT-> zazanzcaamy plik zmieniony dajemy komneatrz do zmiany i wypychamy 


solution Xopt;
		int i = 0;
		int c = a+(b - a) / c;
		vector <double> a_, b_, c_,d_;
		a_.push_back(a);
		b_.push_back(b);
		c_.push_back(c);
		int l, m;
		l = ff1T(a_[i]) * (pow(b_[i], 2) - pow(c_[i], 2)) + ff1T(b[i]) * (pow(c_[i], 2) - pow(a_[i], 2)) + ff1T(c_[i]) * (pow(a_[i], 2) - pow(b_[i], 2));
		m = ff1T(a_[i]) * (b_[i] - c[i]) + ff1T(b_[i]) * (c_[i] - a_[i]) + ff1T(c_[i]) * (a_[i] - b_[i]);
		if (m <= 0)
			return error
			d_.push_back( 0, 5 * l / m);
		if (a_[i] < d_[i] < c_[i]) {
			if (ff1T(d_[i]) < ff1T(c_[i])) {
				a_.push_back(a_[i]);
				c_.push_back(d_[i]);
				b_.push_back(c_[i]);
			}
			else {
				a_.push_back(d_[i]);
				c_.push_back(c_[i]);
				b_.push_back(b_[i]);
			}
		}
		else{
			if (c_[i] < d_[i] < b_[i]) {
				if (ff1T(d_[i] < ff1T(c_[i])) {
					a_.push_back(c_[i]);
					c_.push_back(d_[i]);
					b_.push_back(b_[i]);
				}
				else {
					a_.push_back(a_[i]);
					c_.push_back(c_[i]);
					b_.push_back(d_[i]);
				}
			}
			else {
				return error;
			}
			i = i + 1;
			if (i > Nmax)
				return error;
			if (b_[i] - a_[i] < epsilon || abs(d_[i] - d_[i - 1]) < y)
				return d_[i];

			}
