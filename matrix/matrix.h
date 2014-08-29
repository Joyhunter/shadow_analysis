#pragma once

#ifdef KUN_SPARSE_MATRIX_DLL_EXPORT

	#define KUN_SPARSE_MATRIX_DLL_API __declspec(dllexport)

#else

	#define KUN_SPARSE_MATRIX_DLL_API __declspec(dllimport)
	
	#ifdef _DEBUG
		#pragma comment(lib,"sparsematd.lib")
	#else
		#pragma comment(lib,"sparsemat.lib")
	#endif

#endif


namespace sparse{

	class KUN_SPARSE_MATRIX_DLL_API matrix
	{
	public:
		matrix(void);
		~matrix(void);
		
		// calling this function (create) would destroy previous data
		bool create(int nRows,int nCols);
		void destroy();

		// user must ensure no duplicated elements are inserted
		// or this method would return false.
		bool add(int rowIdx,int colIdx, double val);

		int getNumofRows() const ;
		int getNumofColumns() const;
		int getNumOfNonZeros() const;

		/*

		it is user's responsibility to choose which method to use

		*/

		// this method solves Ax = b
		// for square matrix, it would use LU decomposition;
		// for rectangular matrix, it would use QR decomposition (solving least square form);
		// users must ensure enough memory of x and b are allocated
		// length of x should be equals to getNumOfColumns();
		// length of b should be equals to getNumOfRows();
		bool solve(double* x, const double* b) const;

		// this method solves A_{t}Ax = A_{t}b
		// In theory, it should get the same result as solving Ax = b;
		// Advantage: solving A_{t}Ax = A_{t}b may be faster than solving Ax = b, for rectangular matrix A
		// Disadvantage: but it is more numerical instable than solving Ax = b.
		// users must ensure enough memory of x and b are allocated
		// length of x should be equals to getNumOfColumns();
		// length of b should be equals to getNumOfRows();
		bool solveSym(double* x, const double* b) const;
		
		// calculate b = Ax
		// users must ensure enough memory of x and b are allocated
		// length of x should be equals to getNumOfColumns();
		// length of b should be equals to getNumOfRows();
		bool multiply(const double* x, double* b) const;

	protected:
		int m_nCols;
		int m_nRows;
		
		void * _p;
		
	};

	class KUN_SPARSE_MATRIX_DLL_API symmatrix
	{
	public:
		symmatrix();
		~symmatrix();

		bool create(int nRows);
		void destroy();

		//user must ensure no duplicated elements are inserted
		//e.g. : if user insert elements at (1,2) and at (2,1), it would be considered as dulipcated
		bool add(int rowIdx,int colIdx, double val);

		int getNumofRows() const;

		//only lower triangular parts are counted in.
		int getNumOfNonZeros() const;
		
		// this method solves Ax = b
		// Since A is symmetric,
		// It would first try cholesky decomposition (CHOLMOD) to solve;
		// If failed, it would convert A to a normal unsymmetric matrix and use LU decomposition to solve.
		// users must ensure enough memory of x and b are allocated
		// length of x should be equals to getNumOfColumns();
		// length of b should be equals to getNumOfRows();
		bool solve(double*x, const double*b) const;
		
		// calculate b = Ax
		// users must ensure enough memory of x and b are allocated
		// length of x should be equals to getNumOfColumns();
		// length of b should be equals to getNumOfRows();
		bool multiply(const double*x, double*b) const;

		bool convert2unsymmatrix(matrix& outMat) const;

	protected:
		int m_nRows;

		// only lower triangular parts are stored (e.g. row index > column index)
		void * _p;
				
	};
}

