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
		const kMatrixView<V>	A,		//	n x 3
		const kVectorView<V>	r,
		kVectorView<V>			u,
		kVectorView<V>			gam)
	{
		//	helps
		V bet;
		int j;

		//	dim
		int n = A.rows();

		//	go
		u(0) = r(0)/(bet=A(0,1));
		for(j=1;j<n;++j)
		{
			gam(j) = A(j-1,2)/bet;
			bet    = A(j,1)-A(j,0)*gam(j);
			u(j)   = (r(j)-A(j,0)*u(j-1))/bet;
		}
		for(j=n-2;j>=0;--j)
		{
			u(j) -= gam(j+1)*u(j+1);
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
		//	dim
		int n = A.rows();

		//	check dim
		if(u.size()<n) u.resize(n);
		if(gam.size()<n) gam.resize(n);

		tridag(A(), r(), u(), gam());
		//	done
		return;
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
		int n = A.rows()-1;
		x.resize(n+1,m1+1+m2);
		V xi;
		for(int i = 0;i<=n;++i)
		{
			int jl = max<int>(0, i - m1);
			int ju = min<int>(i + m2, n);
			xi = 0.0;
			for(int j = jl;j<=ju;++j)
			{
				int k = j - i + m1;
				xi += A(i,k)*b(j);
			}
			x(i) = xi;
		}

		//	done
		return;
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

	//	solve: A u = r where A is left diag
	template <class V>
	void leftdag(
		const kMatrix<V>&	A,	//	nx2
		const kVector<V>&	r,
		kVector<V>&			u)
	{
		//	resize
		u.resize(A.rows());

		//	tjek
		if(!A.rows()) return;

		//	dims
		int n = A.rows() - 1;

		//	helps
		int i;

		//	loop
		u(0) = r(0)/A(0, 1);
		for(i=1;i<=n;++i)
		{
			u(i) = (r(i) - A(i, 0) * u(i-1))/A(i,1);
		}

		//	done
		return;
	}

	//	solve: A u = r where A is right diag
	template <class V>
	void rightdag(
		const kMatrix<V>&	A,	//	nx2
		const kVector<V>&	r,
		kVector<V>&			u)
	{
		//	resize
		u.resize(A.rows());

		//	tjek
		if (!A.rows()) return;

		//	dims
		int n = A.rows() - 1;

		//	helps
		int i;

		//	loop
		u(n) = r(n) / A(n, 0);
		for(i=n-1;i>=0;--i)
		{
			u(i) = (r(i) - A(i,1)*u(i+1))/A(i,0);
		}

		//	done
		return;
	}


}
