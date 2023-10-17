#pragma once

#include <vector>
#include <utility>

using namespace std;

template<typename T>
vector<T> diag3MatrixMul(const vector<T>& a, const vector<T>& b, const vector<T>& c, const vector<T>& v) {
    int n = v.size();

    vector<T> res(n);

    if (n == 1) {
        return {b[0] * v[0]};
    }

    res[0] = b[0] * v[0] + c[0] * v[1];
    res[n - 1] = a[n - 1] * v[n - 2] + b[n - 1] * v[n - 1];

    for (int i = 1; i < n - 1; i++) {
        res[i] = a[i] * v[i - 1] + b[i] * v[i] + c[i] * v[i + 1];
    }

    return res;
}
