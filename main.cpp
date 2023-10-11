#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include "QuadMatrix.h"
#include "Matrix.h"
#include "Norma.h"
#include "IterativeLinSolveAlgs.h"

using namespace std;

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

    auto bestTau = findBestIterParam<double>(A, b, opts);

    cout << bestTau << endl;

    return 0;
}
