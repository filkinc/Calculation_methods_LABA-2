#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include "QuadMatrix.h"

//Кубическая норма для вектора и матрицы
template <class T>
T norm_inf(vector<T> x) {
	size_t n = x.size();
	T xMax = 0;

	for (int i = 0; i < n; ++i) {
		if (fabs(x[i]) > xMax) {
			xMax = fabs(x[i]);
		}
	}
	return xMax;
}
template <class T>
T norm_inf(QuadMatrix<T> A) {
	size_t n = A.order();
	T  aSum, aMax = 0;

	for (int i = 0; i < n; ++i) {
		aSum = 0;
		for (int j = 0; j < n; j++) {
			aSum += fabs(A(i, j));
		}
		if (aSum > aMax) {
			aMax = aSum;
		}
	}
	return aMax;
}

//Октаэдральная норма для вектора и матрицы
template<class T>
T norm_1(vector<T> x) {
	size_t n = x.size();
	T xSum = 0;
	for (int i = 0; i < n; ++i) {
		xSum += fabs(x[i]);
	}
	return xSum;
}
template<class T>
T norm_1(QuadMatrix<T> A) {
	size_t n = A.order();
	T aSum, aMax = 0;
	for (int j = 0; j < n; ++j) {
		aSum = 0;
		for (int i = 0; i < n; ++i) {
			aSum += fabs(A(i, j));
		}
		if (aSum > aMax) {
			aMax = aSum;
		}
	}
	return aMax;
}

//Шаровая норма для вектора и матрицы(норма Фробениуса)
//template<class T>
//T norm_2(vector<T> x) {
//	size_t n = x.size();
//	T xSum = 0;
//	for (int i = 0; i < n; ++i) {
//		xSum += x[i] * x[i];
//	}
//	return sqrt(xSum);
//}
//Это не точно
//template<class T>
//T norm_F(QuadMatrix<T> A) {
//	size_t n = A.order();
//	T aSum = 0;
//	for (int i = 0; i < n; ++i) {
//		for (int j = 0; j < n; ++j) {
//			aSum += A(i, j) * A(i, j);
//		}
//	}
//	return sqrt(aSum);
//}