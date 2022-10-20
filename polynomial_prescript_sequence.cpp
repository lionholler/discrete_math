#include <vector>
#include <iostream>
#include <cmath>

using namespace std;

const double EPS = 1e-9;
const int INF = 2;

long double fac(const int n) {
	long double r = 1;
	for (int i = n; i > 0; i--) {
		r *= i;
	}
	return r;
}

long double pow(const int b, const int x) {
	long double r = 1;
	for (int i = x; i > 0; i--) {
		r *= b;
	}
	return r;
}

int main() {

	/**		*		*		*		*		**

	 Prediction algorithm for the polynomial
	 prescript  of  a  polynomial  prescript 
	 sequence.  Using  a   finite difference 
	 table   and  Gauss-Jordan  method   for 
	 solving  system   of  linear equations.
	 If  a  constant term  can be found from 
	 the finite difference table a system of
	 linear equations  can be  build solving 
	 into the coefficients of  the resulting 
	 polynomial  prescript  of the sequence.

	 Coefficient a is precalculated with the
	 formula  *a* = *(cte term)* /  *power*!

	 assumptions:
	 	- the given sequence starts at A0
	 	  allowing  to  immediately  find 
	 	  coefficient d.
	 	- suffiencent  terms are given to 
	 	  reach a constant row.
	 	  min  required  sequence length:
	 	  power of polynom + 2.
	 	  longer  lengths  are encouraged 
	 	  to   avoid  false  convergence.

	 compile w/:
		g++ -std=c++11 -O3 xx.cpp -o xx

	 complexity:
	    - Gauss-Jordan:  if  m=n  O(n^3)
	      generally   O(min(n,m)  *  nm)

	 **		*		*		*		*		**/
	
	vector<double> in{-3,-3,-3,0,9,27,57};
	const vector<double> in_ = in;

	bool is_const;
	double p, CTE_T, POWER;
	for (int i = 0; i < 100; i++) {
		is_const = true;
		for (int i = 1; i < in.size(); i++) {
			in[i-1] = in[i] - in[i-1];
			if (i == 1) {
				p = in[i-1];
				continue;
			} else if (is_const && abs(p - in[i-1]) > EPS) {
				is_const = false;
			}
			p = in[i-1];
		}
		in.resize(in.size() - 1);
		if (in.size() < 2) {
			cout << "end without conversion to constant form, too few terms. Exit 0." << endl;
			exit(0);
		}
		if (is_const) {
			CTE_T = in[0];
			POWER = i+1;
			break;
		}
	}

	double a, d;
	a = CTE_T / fac(POWER);
	d = in_[0];

	cout << endl << '\t' << "Cte term: " << CTE_T << " power: " << POWER << " a: " << a << " d: " << d << endl;

	vector<vector<double>> system(max(4,(int)POWER+1), vector<double>((int)POWER+2, 1));
	for (int i = 0; i < system.size(); i++) {
		for (int j = 0; j < system[i].size(); j++) {
			if (i == 1) {
				if (j == 0) {
					system[i][j] = 1;
				} else if (j == system[i].size()-1) {
					system[i][j] = a;
				} else {
					system[i][j] = 0;
				}
				continue;
			}
			if (j == system[i].size()-1) {
				system[i][j] = in_[i];
				continue;
			}
			system[i][j] = pow(i, system[i].size()-j-2);
		}
	}

	vector<double> coeff(system[0].size() - 1, 0.);

	auto gauss = [&]() {
		int n = (int) system.size();
		int m = (int) system[0].size()-1;
		vector<int> w(m, -1);
		for (int c = 0, r = 0; c < m && r < n; ++c) {
			int select = r;
			for (int i = r; i < n; ++i) {
				if (abs(system[i][c]) > abs(system[select][c])) {
					select = i;
				}
			}
			if (abs(system[select][c]) < EPS) {
				continue;
			}
			for (int i = c; i <= m; i++) {
				swap(system[select][i], system[r][i]);
			}
			w[c] = r;
			for (int i = 0; i < n; i++) {
				if (i != r) {
					const double div = system[i][c] / system[r][c];
					for (int j = c; j <= m; j++) {
						system[i][j] -= system[r][j] * div;
					}
				}
			}
			++r;
		}
		for (int i = 0; i < m; i++) {
			if (w[i] != -1) {
				coeff[i] = system[w[i]][m] / system[w[i]][i];
			}
		}
		for (int i = 0; i < n; i++) {
			double sum = 0;
			for (int j = 0; j < m; j++) {
				sum += coeff[j] * system[i][j];
			}
			if (abs(sum - system[i][m]) > EPS) {
				return 0;
			}
		}
		for (int i = 0; i < m; ++i) {
			if (w[i] == -1) {
				return INF;
			}
		}
		return 1;
	};

	int ret = gauss();
	if (ret == 0) {
		cout << "equation system has no solution, exit 0" << endl;
		exit(0);
	} else if (ret == INF) {
		cout << "equation system has infinite many solutions, exit 0" << endl;
		exit(0);
	}

	cout << endl << '\t' << "An = ";
	for (int i = 0; i < coeff.size(); i++) {
		if (abs(0 - coeff[i]) < EPS) {
			continue;
		}
		cout << coeff[i] * (coeff[i] < 0 ? -1 : 1);
		
		if (i < coeff.size() - 1) {
			if (coeff.size()-i-1 == 1) {
				cout << " * n";
			} else {
				cout << " * (n ^ " << coeff.size()-i-1 << ")";
			}
			if (coeff[i+1] < 0 && abs(0 - coeff[i+1]) > EPS) {
				cout << " - ";
			} else if (abs(0 - coeff[i+1]) > EPS) {
				cout << " + ";
			}
		}
	}
	cout << endl << endl;;

	return 0;
}
