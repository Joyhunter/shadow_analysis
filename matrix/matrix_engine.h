#pragma once

#include "matrix.h"

namespace Utility
{

	class MatrixFunction
	{
	public:
		void solveSparseLinearSystem(
			int rows, int cols, 
			vector<int>& A_I, vector<int>& A_J, vector<double>& A_Val, 
			vector<double>& b, 
			vector<double>& x)
		{
			sparse::matrix A;
			A.create(rows, cols);

			doFv(i, A_I){
				A.add(A_I[i], A_J[i], A_Val[i]);
			}

			double* b_ = new double[b.size()];
			doFv(i, b) b_[i] = b[i];
			double* x_ = new double[x.size()];

			A.solve(x_, b_);
			
			doFv(i, x) x[i] = x_[i];

			delete [] b_;
			delete [] x_;
		}
	};

}