#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include "QuadMatrix.h"
#include "Matrix.h"
#include "Norma.h"

//template<class T> T permEps() { return T(); }
//template<> float permEps<float>()   { return (float)1e-6; }
//template<> double permEps<double>() { return 1e-15; }

template<class T>
pair<vector<T>, QuadMatrix<T>> gaussLinSolve(QuadMatrix<T> A, vector<T> b) {
	size_t n = A.order();

	QuadMatrix<T> C(n);
	vector<T> y(n);

	for (int i = 0; i < n; ++i) {
		int m = i;
		for (int k = i; k < n; ++k) {
			if (fabs(A(k, i)) > fabs(A(m, i))) {
				m = k;
			};
		}
		for (int k = 0; k < n; ++k) {
			swap(A(m, k), A(i, k));
		}
		swap(b[i], b[m]);

		if (fabs(A(i, i)) <= 1e-7) {
			cout << "Нулевой элемент на главной диагонали!";
			system("pause");
			exit(1);
		}
		
		y[i] = b[i] / A(i, i);

		C(i, i) = 1;
		
		for (int j = i + 1; j < n; ++j) {
			C(i, j) = A(i, j) / A(i, i);
		}
		for (int pi = i + 1; pi < n; ++pi) {
			b[pi] -= A(pi, i) * y[i];
			for (int pj = i + 1; pj < n; ++pj) {
				A(pi, pj) -= A(pi, i) * C(i, pj);
			}
		}	
	}
	/*cout << endl;
	cout << "Прямой метод Гаусса: ";
	for (int i = 0; i < n; ++i) {
		cout << std::endl;
		for (int j = 0; j < n; ++j) {
			cout << C(i, j) << ' ';
		}
	}*/
	//C.print();
	//cout << endl;

	/*ofstream ansFile;
	ansFile.open("AnswerFile.txt", ios_base::out);
	ansFile << endl;
	ansFile << "Прямой метод Гаусса: ";
	for (int i = 0; i < n; ++i) {
		ansFile << endl;
		for (int j = 0; j < n; ++j) {
			ansFile << C(i, j) << ' ';
		}
	}
	ansFile << endl;
	ansFile.close();*/
	

	return { upperTriagLinSolve(C, y), C };
}

template<class T>
vector<T> upperTriagLinSolve(QuadMatrix<T> A, vector<T> b) {
	size_t n = A.order();
	vector<T> x(n);
	
	for (int i = n - 1; i >= 0; i--) {
		if (fabs(A(i, i)) <= 1e-7) {
			cout << "Нулевой элемент на главной диагонали!";
			system("pause");
			exit(1);
		}
		x[i] = b[i];
		for (int j = i + 1; j < n; j++) {
			x[i] -= A(i, j) * x[j];
		}

		x[i] /= A(i, i);
	}
	return x;
}

template<class T>
pair<QuadMatrix<T>, QuadMatrix<T>> qrDecomposition(QuadMatrix<T> A) {
	size_t n = A.order();
	vector<T> x(n);
	QuadMatrix<T> Q(n);
	T c, s;

	for (int i = 0; i < n; ++i) {
		Q(i, i) = 1;
	}

	
	for (int i = 0; i < n; ++i) {
			int m = i;
			for (int k = i; k < n; ++k) {
				if (fabs(A(k, i)) > fabs(A(m, i))) {
					m = k;
				};
			}
			for (int k = 0; k < n; ++k) {
				swap(A(m, k), A(i, k));
				swap(Q(k, m), Q(k, i));
			}
			
			if (fabs(A(i, i)) <= 1e-7) {
				cout << "Нулевой элемент на главной диагонали!";
				system("pause");
				exit(1);
			}

		for (int j = i + 1; j < n; ++j) {
			c = (A(i, i)) / sqrt(A(i, i) * A(i, i) + A(j, i) * A(j, i));
			s = (A(j, i)) / sqrt(A(i, i) * A(i, i) + A(j, i) * A(j, i));

			for (int k = 0; k < n; ++k) {
				T aa = A(i, k);
				T ab = A(j, k);

				A(i, k) = c * aa + s * ab;
				A(j, k) = c * ab - s * aa;

				T qa = Q(k, i);
				T qb = Q(k, j);

				Q(k, i) = c * qa + s * qb; 
				Q(k, j) = c * qb - s * qa;
			}
			A(j, i) = 0;
		}
	}
	return { Q, A };
}

template<class T>
vector<T> qrLinSolve(QuadMatrix<T> Q, QuadMatrix<T> A, vector<T> b) {
	size_t n = A.order();

	vector<T> br = gaussLinSolve(Q, b).first;
	vector<T> x = upperTriagLinSolve(A, br);

	return x;
}

template<class T> 
T cond(QuadMatrix<T> A, T (&f)(QuadMatrix<T> A)) {
	QuadMatrix inversA = A.inv();
	return f(inversA) * f(A);
}

template<class T>
vector<T> sum(const vector<T>& a, const vector<T>& b) {
	int n = min(a.size(), b.size());
	vector<T> res(n);

	for (int i = 0; i < n; ++i) {
		res[i] = a[i] + b[i];
	}

	return res;
}

template<class T>
vector<T> diff(const vector<T>& a, const vector<T>& b) {
	int n = min(a.size(), b.size());
	vector<T> res(n);

	for (int i = 0; i < n; ++i) {
		res[i] = a[i] - b[i];
	}

	return res;
}

template<class T>
T condEstimate(QuadMatrix<T> matrix, T(&f)(vector<T>)) {
	int n = matrix.order();
	vector<T> kx(n, 1);
	auto kb = mul<T>(matrix, kx);

	cout << endl << "kb: ";
	for (auto bi : kb) {
		cout << bi << ' ';
	}
	cout << endl;

	vector<T> db(n, 0.1);
	auto px = gaussLinSolve(matrix, sum(kb, db)).first;
	auto dx = diff(px, kx);

	return f(dx) * f(kb) / f(db) / f(kx);
}

template<class T> 
T normDiscrepancyVectorGauss(QuadMatrix<T> A, vector<T> b, T(&f)(vector<T>)) {
	vector<T> x = gaussLinSolve(A, b).first;
	vector<T> b1 = mul(A, x);
	return f(diff(b,b1));
}

template<class T>
T normDiscrepancyVectorQR(QuadMatrix<T> A, vector<T> b, T(&f)(vector<T>)) {
	QuadMatrix<T> Q = qrDecomposition(A).first;
	QuadMatrix<T> R = qrDecomposition(A).second;
	vector<T> x = qrLinSolve(Q, R, b);
	vector<T> b1 = mul(A, x);
	return f(diff(b, b1));
}