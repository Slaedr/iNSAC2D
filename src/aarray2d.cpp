/**
 * @file aarray2d.cpp
 * @brief Some method implementations for the 2d arrays class.
 * 
 * Part of FVENS.
 * @author Aditya Kashi
 */

#include "aarray2d.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>

namespace amat {

template <typename T>
Array2d<T>::Array2d() : nrows{0}, ncols{0}, size{0}, elems{nullptr}
{ }

// Full-arg constructor
template <typename T>
Array2d<T>::Array2d(a_int nr, a_int nc)
{
	assert(nr > 0);
	assert(nc > 0);
	nrows = nr; ncols = nc;
	size = nrows*ncols;
	elems = new T[nrows*ncols];
}

template <typename T>
Array2d<T>::Array2d(const Array2d<T>& other)
{
	nrows = other.nrows;
	ncols = other.ncols;
	size = nrows*ncols;
	elems = new T[nrows*ncols];
	for(a_int i = 0; i < nrows*ncols; i++)
	{
		elems[i] = other.elems[i];
	}
}

template <typename T>
Array2d<T>::~Array2d()
{
	delete [] elems;
}

template <typename T>
Array2d<T>& Array2d<T>::operator=(const Array2d<T>& rhs)
{
#ifdef DEBUG
	if(this==&rhs) return *this;		// check for self-assignment
#endif
	nrows = rhs.nrows;
	ncols = rhs.ncols;
	size = nrows*ncols;
	delete [] elems;
	elems = new T[nrows*ncols];
	for(a_int i = 0; i < nrows*ncols; i++)
	{
		elems[i] = rhs.elems[i];
	}
	return *this;
}

/// Separate setup function in case no-arg constructor has to be used
/** \deprecated Please use resize() instead.
 */
template <typename T>
void Array2d<T>::setup(const a_int nr, const a_int nc)
{
	assert(nr > 0);
	assert(nc > 0);
	nrows = nr; ncols = nc;
	size = nrows*ncols;
	delete [] elems;
	elems = new T[nrows*ncols];
}

/// Sets a new size for the array, deletes the contents and allocates new memory
template <typename T>
void Array2d<T>::resize(const a_int nr, const a_int nc)
{
	assert(nr > 0);
	assert(nc > 0);
	nrows = nr; ncols = nc;
	size = nrows*ncols;
	delete [] elems;
	elems = new T[nrows*ncols];
}

/// Setup without deleting earlier allocation: use in case of Array2d<t>* (pointer to Array2d<t>)
template <typename T>
void Array2d<T>::setupraw(const a_int nr, const a_int nc)
{
	assert(nr > 0);
	assert(nc > 0);
	nrows = nr; ncols = nc;
	size = nrows*ncols;
	delete [] elems;
	elems = new T[nrows*ncols];
}

/// Fill the matrix with zeros.
template <typename T>
void Array2d<T>::zeros()
{
	for(a_int i = 0; i < size; i++)
		elems[i] = (T)(0.0);
}

template <typename T>
void Array2d<T>::ones()
{
	for(a_int i = 0; i < size; i++)
		elems[i] = 1;
}

template <typename T>
void Array2d<T>::identity()
{
	T one = (T)(1);
	T zero = (T)(0);
	for(a_int i = 0; i < nrows; i++)
		for(a_int j = 0; j < ncols; j++)
			if(i==j) operator()(i,j) = one;
			else operator()(i,j) = zero;
}

template <typename T>
void Array2d<T>::setdata(const T* A, a_int sz)
{
#ifdef DEBUG
	if(sz != size)
	{
		std::cout << "\nError in setdata: argument size does not match matrix size";
		return;
	}
#endif
	for(a_int i = 0; i < nrows; i++)
		for(a_int j = 0; j < ncols; j++)
			elems[i*ncols+j] = A[i*ncols+j];
}

template <typename T>
void Array2d<T>::mprint() const
{
	std::cout << "\n";
	for(a_int i = 0; i < nrows; i++)
	{
		for(a_int j = 0; j < ncols; j++)
			std::cout << std::setw(WIDTH) << std::setprecision(WIDTH/2+1) << elems[i*ncols+j];
		std::cout << std::endl;
	}
}

template <typename T>
void Array2d<T>::fprint(std::ofstream& outfile) const
{
	//outfile << '\n';
	outfile << std::setprecision(MATRIX_DOUBLE_PRECISION);
	for(a_int i = 0; i < nrows; i++)
	{
		for(a_int j = 0; j < ncols; j++)
			outfile << " " << elems[i*ncols+j];
		outfile << '\n';
	}
}

template <typename T>
void Array2d<T>::fread(std::ifstream& infile)
{
	infile >> nrows; infile >> ncols;
	size = nrows*ncols;
	delete [] elems;
	elems = new T[nrows*ncols];
	for(a_int i = 0; i < nrows; i++)
		for(a_int j = 0; j < ncols; j++)
			infile >> elems[i*ncols + j];
}

template <typename T>
T Array2d<T>::average() const
{
	T avg = 0;
	for(a_int i = 0; i < size; i++)
		avg += elems[i];
	avg = avg/size;
	return avg;
}

template <typename T>
T Array2d<T>::l2norm() const
{
	T tot = 0;
	for(a_int i = 0; i < size; i++)
	{
		tot += elems[i]*elems[i];
	}
	tot = std::sqrt(tot);
	return tot;
}

/*
/// Returns sum of products of respective elements of flattened arrays 
/// containing matrix elements of this and A
template <typename T>
T Array2d<T>::dot_product(const Array2d<T>& A)
{
const T* elemsA = A.elems;
T ans = 0;
for(a_int i = 0; i < size; i++)
{
ans += elems[i]*elemsA[i];
}
return ans;
}

/// Computes 1-norm (max column-sum norm) of the matrix
template <typename T>
T Array2d<T>::matrixNorm_1() const
{
T max = 0, sum;
a_int i,j;
for(j = 0; j < ncols; j++)
{
sum = 0;
for(i = 0; i < nrows; i++)
{
sum += (T)( fabs(get(i,j)) );
}
if(max < sum) max = sum;
}
return max;
}*/

template class Array2d<a_real>;
template class Array2d<a_int>;

}
