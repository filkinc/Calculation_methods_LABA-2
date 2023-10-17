#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include <numeric>
#include "QuadMatrix.h"
#include "Matrix.h"
#include "Norma.h"
#include "LinSolveAlgs.h"
#include "Diag3Matrix.h"

using namespace std;

enum IterMethodResCodes {IMRC_OK, IMRC_NORM_GE_1, IMRC_ITER_OVERFLOW, IMRC_BAD_MATRIX};

template<typename T>
struct iterativeMethodResult {
    int iterationCount;
    vector<T> sol;
    QuadMatrix<T> C;
    vector<vector<T>> xs;
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
        cerr << "Ќорма матрицы C должна быть меньше 1!" << endl << endl;
        return {0, {}, C, {}, IMRC_NORM_GE_1};
    }

    int itCount = 0;

    vector<T> prevX(n, 0);
    vector<T> curX = prevX;
    vector<vector<T>> xs = {curX};

    for (itCount = 1; itCount <= 1000 && opts.vNorm(diff(mul(A, curX), b)) > opts.eps; ++itCount) {
        prevX = curX;
        curX = sum(mul(C, prevX), y);

        xs.push_back(curX);
    }

    if (itCount >= 1000) {
        cerr << " оличество итераций привысело заданный порог" << endl << endl;
        return {itCount, curX, C, xs, IMRC_ITER_OVERFLOW};
    }

    return {itCount, curX, C, xs, IMRC_OK};
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
        cerr << "Ќи дл€ одног из рассмотренных итерационных параметров метод не сходитс€!" << endl;
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
            cerr << "ћатрица система должна содержать ненулевые диагональные элементы!" << endl;
            return {0, {}, {{}}, {}, IMRC_BAD_MATRIX}; // ????. ???????? ?????? ???? ?? ?????? > 0.
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

template<typename T>
iterativeMethodResult<T> relaxationMethod(const QuadMatrix<T>& A, const vector<T>& b, T w, IterLinSolveOptions<T> opts, bool useStationaryMethod = false) {
    int n = A.order();

    QuadMatrix<T> C(n);
    QuadMatrix<T> C1(n), C2(n);

    for (int i = 0; i < n; ++i) {
        if (A(i, i) < opts.eps) {
            cerr << "ћатрица система должна содержать ненулевые диагональные элементы!" << endl;
            return {0, {}, {{}}, {}, IMRC_BAD_MATRIX}; // ????. ???????? ?????? ???? ?? ?????? > 0.
        }
    }

    for (int i = 0; i < n; ++i) {
        C1(i, i) = 1;
        for (int j = 0; j < i; ++j) {
            C1(i, j) = w * A(i, j) / A(i, i);
        }
    }

    for (int i = 0; i < n; ++i) {
        C2(i, i) = 1 - w;
        for (int j = i + 1; j < n; ++j) {
            C2(i, j) = -w * A(i, j) / A(i, i);
        }
    }

    C = C1.inv() * C2;

    if (opts.mNorm(C) >= 1) {
        cerr << "Ќорма матрицы C должна быть меньше 1!" << endl << endl;
        return {0, {}, C, {}, IMRC_NORM_GE_1};
    }

    if (useStationaryMethod) {
        vector<T> y(n);
        for (int i = 0; i < n; ++i) {
            y[i] = w * b[i] / A(i, i);
        }

        y = mul(C1.inv(), y);

        return stationaryIterativeMethod(A, b, C, y, opts);
    }

    int itCount = 0;

    vector<T> prevX(n, 0);
    vector<T> curX = prevX;
    vector<vector<T>> xs = {curX};

    for (itCount = 1; itCount <= 1000 && opts.vNorm(diff(mul(A, curX), b)) > opts.eps; ++itCount) {
        prevX = curX;

        for (int i = 0; i < n; ++i) {
            T prevXsum = 0;
            T curXsum = 0;

            for (int j = 0; j < i; ++j) {
                curXsum += curX[j] * A(i, j) / A(i, i);
            }

            for (int j = i + 1; j < n; ++j) {
                prevXsum += prevX[j] * A(i, j) / A(i, i);
            }

            curX[i] = -w * curXsum + (1 - w) * prevX[i] - w * prevXsum + w * b[i] / A(i, i);
        }

        xs.push_back(curX);
    }

    if (itCount >= 1000) {
        cerr << " оличество итераций привысело заданный порог" << endl << endl;
        return {itCount, curX, C, xs, IMRC_ITER_OVERFLOW};
    }

    return {itCount, curX, C, xs, IMRC_OK};
}

template<typename T>
iterativeMethodResult<T> relaxation3diagMethod(const vector<T>& a,
                                               const vector<T>& b,
                                               const vector<T>& c,
                                               const vector<T>& d,
                                               T w,
                                               IterLinSolveOptions<T> opts)
{
    int n = a.size();

    for (int i = 0; i < n; ++i) {
        if (b[i] < opts.eps) {
            cerr << "ћатрица система должна содержать ненулевые диагональные элементы!" << endl;
            return {0, {}, {{}}, {}, IMRC_BAD_MATRIX}; // ????. ???????? ?????? ???? ?? ?????? > 0.
        }
    }

    int itCount = 0;

    vector<T> prevX(n, 0);
    vector<T> curX = prevX;
    vector<vector<T>> xs = {curX};



    for (itCount = 1; itCount <= 1000 && opts.vNorm(diff(diag3MatrixMul(a, b, c, curX), b)) > opts.eps; ++itCount) {
        prevX = curX;

        for (int i = 0; i < n; ++i) {
            T prevXsum = i == n - 1 ? 0 : prevX[i + 1] * c[i] / b[i];
            T curXsum = i == 0 ? 0 : curX[i - 1] * a[i] / b[i];

            curX[i] = -w * curXsum + (1 - w) * prevX[i] - w * prevXsum + w * d[i] / b[i];
        }

        xs.push_back(curX);
    }

    if (itCount >= 1000) {
        cerr << " оличество итераций привысело заданный порог" << endl << endl;
        return {itCount, curX, {{}}, xs, IMRC_ITER_OVERFLOW};
    }

    return {itCount, curX, {{}}, xs, IMRC_OK};
}

template<typename T>
iterativeMethodResult<T> seidelMethod(const QuadMatrix<T>& A, const vector<T>& b, IterLinSolveOptions<T> opts, bool useStationaryMethod = false) {
    return relaxationMethod<T>(A, b, 1, opts, useStationaryMethod);
}

template<typename T>
iterativeMethodResult<T> seidel3diagMethod(const vector<T>& a,
                                          const vector<T>& b,
                                          const vector<T>& c,
                                          const vector<T>& d,
                                          IterLinSolveOptions<T> opts) {
    return relaxation3diagMethod<T>(a, b, c, d, 1, opts);
}
