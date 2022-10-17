#pragma once

#include <vector>
#include <iostream>
#include <thread>

using namespace std;

typedef long long int64;

struct MatrixSize
{
	size_t rows;
	size_t cols;

	MatrixSize(size_t rows, size_t cols) : rows(rows), cols(cols)
	{ }
};

template <typename Ty>
class Matrix
{
private:
	size_t _n;
	size_t _m;
	vector<vector<Ty>> _matrix;
public:
	Matrix(size_t n, size_t m)
	{
		_n = n;
		_m = m;
		_matrix.assign(_n, vector<Ty>(_m));
	}

	Matrix(vector<vector<Ty>> &matrix)
	{
		_n = matrix.size();
		if (_n < 1)
		{
			throw new invalid_argument("Invalid matrix size.");
		}
		
		_m = matrix.front().size();
		if (_m < 1)
		{
			throw new invalid_argument("Invalid matrix size.");
		}

		_matrix = vector<vector<Ty>>(matrix);
	}

	MatrixSize size()
	{
		return MatrixSize(_n, _m);
	}

	Ty at(size_t i, size_t j)
	{
		return _matrix[i][j];
	}

	void setElement(size_t i, size_t j, Ty val)
	{
		_matrix[i][j] = val;
	}

	Ty det()
	{
		if (_n != _m)
		{
			throw new invalid_argument("Matrix is not square.");
		}
		if (_n == 1)
			return _matrix[0][0];
		if (_n == 2)
			return _matrix[0][0] * _matrix[1][1] - _matrix[0][1] * _matrix[1][0];
		Ty result = 0;
		for (int i = 0; i < _n; i++)
		{
			result += _matrix[i][0] * cofactor(i, 0);
		}

		return result;
	}

	Matrix* parallelInverse()
	{
		Ty determinant = this->det();
		if (determinant == 0)
		{
			throw new invalid_argument("There is no matrix inverse. Det = 0.");
		}
		Matrix* cMatrix = parallelCofactorMatrix();
		Matrix* tCofactorMatrix = cMatrix->T();
		Matrix* result = tCofactorMatrix->dot(1.0 / determinant);
		delete cMatrix;
		delete tCofactorMatrix;
		return result;
	}

	Matrix* inverse()
	{
		Ty determinant = this->det();
		if (determinant == 0)
		{
			throw new invalid_argument("There is no matrix inverse. Det = 0.");
		}
		Matrix* cMatrix = cofactorMatrix();
		Matrix* tCofactorMatrix = cMatrix->T();
		Matrix* result = tCofactorMatrix->dot(1.0 / determinant);
		delete cMatrix;
		delete tCofactorMatrix;
		return result;
	}

	Ty cofactor(size_t row, size_t col)
	{
		vector<vector<Ty>> minorVectors(_n - 1, vector<Ty>(_m - 1));
		for (int i = 0, i1 = 0; i < _n; i++)
		{
			if (i == row)
				continue;
			for (int j = 0, j1 = 0; j < _m; j++)
			{
				if (j == col)
					continue;
				minorVectors[i1][j1] = _matrix[i][j];
				j1++;
			}
			i1++;
		}

		Matrix minorMatrix(minorVectors);
		Ty minorMatrixDet = minorMatrix.det();
		return ((row + col) % 2 == 0) ? minorMatrixDet : -minorMatrixDet;
	}

	Matrix* parallelCofactorMatrix()
	{
		Matrix<Ty>* resultMatrix = new Matrix(_n, _m);
		vector<thread> threads;
		for (int i = 0; i < _n; i++)
		{
			for (int j = 0; j < _m; j++)
			{
				threads.push_back(thread(&Matrix<Ty>::calculateCofactor, resultMatrix, this, i, j));
			}
		}

		for (auto& th : threads)
		{
			th.join();
		}
		return resultMatrix;
	}

	Matrix* cofactorMatrix()
	{
		Matrix* cofactorMatrix = new Matrix(_n, _m);
		for (int i = 0; i < _n; i++)
		{
			for (int j = 0; j < _m; j++)
			{
				cofactorMatrix->calculateCofactor(this, i, j);
			}
		}
		return cofactorMatrix;
	}

	Matrix* T()
	{
		vector<vector<Ty>> vectors(_n, vector<Ty>(_m));
		for (int i = 0; i < _n; i++)
		{
			for (int j = 0; j < _m; j++)
			{
				vectors[i][j] = _matrix[j][i];
			}
		}
		Matrix* tMatrix = new Matrix(vectors);
		return tMatrix;
	}

	void calculateCofactor(Matrix* matrix, size_t i, size_t j)
	{
		this->setElement(i, j, matrix->cofactor(i, j));
	}

	Matrix<double>* dot(double factor)
	{
		vector<vector<double>> product(_n, vector<double>(_m));
		for (int i = 0; i < _n; i++)
		{
			for (int j = 0; j < _m; j++)
			{
				product[i][j] = double(factor * _matrix[i][j]);
			}
		}
		Matrix* result = new Matrix(product);
		return result;
	}
};

