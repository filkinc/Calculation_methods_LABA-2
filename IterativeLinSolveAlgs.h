#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include "QuadMatrix.h"
#include "Matrix.h"
#include "Norma.h"
#include "LinSolveAlgs.h"

using namespace std;

enum IterMethodResCodes {IMRC_OK, IMRC_NORM_GE_1, IMRC_ITER_OVERFLOW, IMRC_BAD_MATRIX};

template<typename T>
struct iterativeMethodResult {
    int iterationCount;
    vector<T> sol;
    QuadMatrix<T> C;
    IterMethodResCodes resCode;
};

template<typename T>
struct IterLinSolveOptions {
    vectorNorm<T> vNorm;
    matrixNorm<T> mNorm;
    T eps;
};


template<typename T>
iterativeMethodResult<T> stationaryIterativeMethod(const QuadMatrix<T>& A, const vector<T>& b,
                                                   const QuadMatrix<T>& C, const vector<T>& y,
                                                   IterLinSolveOptions<T> opts)
{
    int n = C.order();

    if (opts.mNorm(C) >= 1) {
        cerr << "Норма матрицы не меньше 1. Метод может не сходится!" << endl << endl;
        return {0, {}, {{}}, IMRC_NORM_GE_1};
    }

    int itCount = 0;

    vector<T> prevX(n, 0);
    vector<T> curX = prevX;

    for (itCount = 1; itCount <= 1000 && opts.vNorm(diff(mul(A, curX), b)) > opts.eps; ++itCount) {
        prevX = curX;
        curX = sum(mul(C, prevX), y);
    }

    if (itCount >= 1000) {
        cerr << "Количетсво итераций превысило допустимый порог!" << endl << endl;
        return {itCount, curX, C, IMRC_ITER_OVERFLOW};
    }

    return {itCount, curX, C, IMRC_OK};
}

template<typename T>
iterativeMethodResult<T> simpleIterationMethod(const QuadMatrix<T>& A, const vector<T>& b, T tau, IterLinSolveOptions<T> opts) {
    int n = A.order();

    QuadMatrix<T> C = -tau * A;
    vector<T> y(n);

    for (int i = 0; i < n; ++i) {
        C(i, i) += 1;
        y[i] = tau * b[i];
    }

    return stationaryIterativeMethod(A, b, C, y, opts);
};

template<typename T>
T findBestIterParam(const QuadMatrix<T>& A, const vector<T>& b, IterLinSolveOptions<T> opts) {
    T bestTau = 0;
    int bestTauIterCount = 1e6;

    for (T p = 0.005; p <= 0.1; p += 0.005) {
        T tau = p;
        iterativeMethodResult<T> res = simpleIterationMethod<T>(A, b, tau, opts);

        if (res.resCode == IMRC_OK && res.iterationCount <= bestTauIterCount) {
            bestTau = tau;
            bestTauIterCount = res.iterationCount;
        }
    }

    if (bestTau == 0) {
        cerr << "Лучший итерационный параметр не был найден!" << endl;
    }

    return bestTau;
}

template<typename T>
iterativeMethodResult<T> jacobiMethod(const QuadMatrix<T>& A, const vector<T>& b, IterLinSolveOptions<T> opts) {
    int n = A.order();

    QuadMatrix<T> C(n);
    vector<T> y(n);

    for (int i = 0; i < n; ++i) {
        if (A(i, i) < opts.eps) {
            cerr << "Диагональные элементы матрицы должны быть ненулевыми!" << endl;
            return {0, {}, {{}}, IMRC_BAD_MATRIX}; // диаг. элементы должны быть по модулю > 0.
        }

        y[i] = b[i] / A(i, i);

        for (int j = 0; j < n; ++j) {
            if (i == j) {
                C(i, j) = 0;
                continue;
            }

            C(i, j) = -A(i, j) / A(i, i);
        }
    }

    return stationaryIterativeMethod(A, b, C, y, opts);
}

