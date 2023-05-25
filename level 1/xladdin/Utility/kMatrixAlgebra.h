#pragma once

//	includes
#include "kMatrix.h"
#include "kInlines.h"

//	class
namespace kMatrixAlgebra
{
	//	mat mult res = A * B
	template <class U, class V, class W>
	void	mmult(
		const kMatrix<U>&	a,
		const kMatrix<V>&	b,
		kMatrix<W>&			ab)
	{
		//	dims
		int m1 = a.rows();
		int m2 = a.cols();
		int m3 = b.cols();

		ab.resize(m1, m3, 0.0);

		//	calc
		for(int i=0;i<m1;++i)
		{
			for(int k=0;k<m2;++k)
			{
				for(int j=0;j<m3;++j)
					ab(i,j) += a(i,k)*b(k,j);
			}
		}

		//	done
		return;
	}

	//	tridag: solves A u = r when A is tridag
	template <class V>
	void	tridag(
		const kMatrix<V>&	    A,		//	n x 3
		const kVector<V>&	    r,
		kVector<V>&				u,
		kVector<V>&				gam)
	{
		// todo: implement
	}

	//	band diagonal matrix vector multiplication
	template <class V>
	void	banmul(
		const kMatrix<V>&	A,		//	n x (m1 + 1 + m2)
		int					m1,
		int					m2,
		const kVector<V>&	b,
		kVector<V>&			x)
	{
		// todo: implement
	}

	//	transposing a symmetric ban matrix
	template <class V>
	void	transpose(
		const int			mm,
		kMatrix<V>&			A)
	{
		int n = A.rows()-1;
		int i, j, k, l;
		int jl;
		for(i=0;i<=n;++i)
		{
			jl = max(0, i - mm);
			for(j=jl;j<i;++j)
			{
				k = j - i + mm;
				l = i - j + mm;
				kInlines::swap(A(i, k), A(j, l));
			}
		}

		//	done
		return;
	}

	//	transposing a ban matrix
	template <class V>
	void	transpose(
		const kMatrix<V>&	A,
		const int			m1,
		const int			m2,
		kMatrix<V>&			At)
	{
		//	resize
		At.resize(A.rows(), A.cols());

		//	loop
		int n = A.rows() - 1;
		int i, j, k, l, jl, ju;
		for (i = 0; i <= n; ++i)
		{
			jl = max(0, i - m1);
			ju = min(i + m2, n);
			for(j=jl;j<=ju;++j)
			{
				k = j - i + m1;
				l = i - j + m2;
				At(j, l) = A(i, k);
			}
		}

		//	done
		return;
	}

}
