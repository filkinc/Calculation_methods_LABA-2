#pragma once

#include<vector>

using namespace std;

template<class T> 
class Matrix
{
protected:
	vector<vector<T>> elems;
public:
	T& operator() (size_t i, size_t j) {
		return elems[i][j];
	}

    T operator() (size_t i, size_t j) const {
        return elems[i][j];
    }

    size_t rowCount() const {
		return elems.size();
	}
    size_t jumCount() const {
		return elems[0].size();
	}

    friend vector<T> operator* (const Matrix<T>& A, vector<T> b);
	friend vector<T> mul (const Matrix<T>& A, vector<T> b);

	Matrix(size_t n, size_t m) {
		elems = vector<vector<T>>(n, vector<T>(m));
    }

    Matrix(vector<vector<T>> data) {
        elems = vector<vector<T>>(data.size());

        for (int i = 0; i < data.size(); ++i) {
            elems[i] = data[i];
        }
    }
};

template<class T>
vector<T> operator* (const Matrix<T>& A, vector<T> b) {
	size_t n = A.rowCount(), m = A.jumCount();
	vector<T> c(m);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			c[i] = A(i, j) * b[j];
		}
	}
	return c;
}

template<class T>
vector<T> mul (const Matrix<T>& A, vector<T> b) {
	size_t n = A.rowCount(), m = A.jumCount();
	vector<T> c(m);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			c[i] = A(i, j) * b[j];
		}
	}
	return c;
}
