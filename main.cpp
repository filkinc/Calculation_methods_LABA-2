#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include <utility>
#include "QuadMatrix.h"
#include "Matrix.h"
#include "Norma.h"
#include "IterativeLinSolveAlgs.h"

using namespace std;

#define zebelTrajectoryFileName "zebelTrajectory.txt"
#define variantNumber 4

template<typename T>
void getVariantDiag3Matrix(vector<T>& a, vector<T>& b, vector<T>& c, vector<T>& d) {
    int n = 200 + variantNumber;
    a = vector<T>(n, 1);
    b = vector<T>(n, 4);
    c = vector<T>(n, 1);
    d = vector<T>(n);
    a[0] = 0;
    c[n - 1] = 0;

    d[0] = 6;
    d[n - 1] = 9 - 3 * (n % 2);

    for (int i = 1; i < n - 1; i++) {
        d[i] = 10 - 2 * ((i + 1) % 2);
    }
}

int main() {
    IterLinSolveOptions<double> opts;
    opts.eps = 1e-7;
    opts.vNorm = norm_inf;
    opts.mNorm = norm_inf;

    QuadMatrix<double> A({
        {15., 2, -3, 7},
        {-5, 11, 2, -3},
        {0, -1, 7, 4},
        {12, 0, -6, 20}
    });

    vector<double> b = {53, -90, 107, 68};

    auto res = simpleIterationMethod<double>(A, b, 0.05, opts);

    res.C.print();

    for (auto& xi : res.sol) {
        cout << xi << ' ';
    }
    cout << endl << endl;


    auto res1 = jacobiMethod<double>(A, b, opts);

    res1.C.print();

    for (auto& xi : res1.sol) {
        cout << xi << ' ';
    }
    cout << endl << endl;

    // auto bestTau = findBestIterParam<double>(A, b, opts);

    // cout << bestTau << endl << endl;

    auto resZ = seidelMethod(A, b, opts);

    resZ.C.print();

    for (auto& xi : resZ.sol) {
        cout << xi << ' ';
    }
    cout << endl << endl;

    ofstream zebelTrajectoryOutput((string)"C:\\Users\\Alex\\Documents\\Qt\\BMSTU\\CalcMethods Semester 5\\CalcMethods_Lab_2\\" + zebelTrajectoryFileName);

    for (auto& x : resZ.xs) {
        for (auto& xi : x) {
            zebelTrajectoryOutput << xi << " ";
        }
        zebelTrajectoryOutput << endl;
    }

    zebelTrajectoryOutput.close();

    auto resR = relaxationMethod(A, b, 1.5, opts);

    resR.C.print();

    for (auto& xi : resR.sol) {
        cout << xi << ' ';
    }
    cout << endl << endl;

    vector<double> a1, b1, c1, d1;
    getVariantDiag3Matrix<double>(a1, b1, c1, d1);

    auto res3diagZ = seidel3diagMethod(a1, b1, c1, d1, opts);

    for (auto& xi : res3diagZ.sol) {
        cout << xi << ' ';
    }
    cout << endl << endl;

    return 0;
}
