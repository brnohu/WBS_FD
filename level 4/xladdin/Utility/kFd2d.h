#pragma once

//	desc:	2d finite difference solution for pdes of the form
//
//		0 = dv/dt + Ax v + Ay v + Axy v
//
//		Ax = -r/2 + mux d/dx + 1/2 varx d^2/dx^2
//
//		Ay = -r/2 + muy d/dy + 1/2 vary d^2/dy^2 
// 
//		Axy = cov d^2/(dx dy)
//
//	using the theta scheme
//
//		[1 - dt Ax] v(t+dt/2) = [ 1 + 1/2 dt Axy] v(t+dt)
//
//		[1 - dt Ay] v(t) = [ 1 + 1/2 dt Axy] v(t+dt/2)
//

//	includes
#include "kFiniteDifference.h"
#include "kMatrixAlgebra.h"

//	class
template <class V>
class kFd2d
{
public:

	//	init
	void	init(
		int					numV,
		const kVector<V>&	x,
		const kVector<V>&	y);

	//	init from sigma rho
	void	initFromSigmaRho(
		const kMatrix<V>&	sigmaX,
		const kMatrix<V>&	sigmaY,
		const kMatrix<V>&	rhoXY);

	//	calc Ax
	void	calcAx(
		int					j,
		V					one,
		V					dt,
		int					wind,
		bool				tr,
		kMatrix<V>&			Ax) const;
	
	//	calc Ay
	void	calcAy(
		int					i,
		V					one,
		V					dt,
		int					wind,
		bool				tr,
		kMatrix<V>&			Ay) const;

	//	calc Axy
	void	calcAxy(
		int					i,
		int					j,
		V					one,
		V					dt,
		bool				tr,
		kMatrix<V>&			Axy) const;

	//	res = Axy b
	V	multAxy(
		int					i,
		int					j,
		const kMatrix<V>&	Axy,
		const kMatrix<V>&	b) const;

	//	roll bwd
	void	rollBwd(
		V						dt,
		int						windX,
		int						windY,
		kVector<kMatrix<V> >&	res);

	//	roll fwd
	void	rollFwd(
		V						dt,
		int						windX,
		int						windY,
		kVector<kMatrix<V> >&	res);

	//	axis'es
	kVector<V>	myX, myY;
	
	//	r, mux, varx, muy, vary, covxy
	kMatrix<V>	myR, myMuX, myVarX, myMuY, myVarY, myCovXY;

	//	result
	kVector<kMatrix<V> >	myRes;

	//	private parts
private:

	//	x operators
	kMatrix<V>	myDxd, myDxu, myDx, myDxx;

	//	y operators
	kMatrix<V>	myDyd, myDyu, myDy, myDyy;

	//	matrixes
	kMatrix<V>	myAx, myAy, myAxy;
	kVector<V>	myWs;
	
	//	value helpers
	kVector<kMatrix<V> >	myVx, myVy, myVxy;
	


};

//	init
template <class V>
void	
kFd2d<V>::init(
	int					numV,
	const kVector<V>&	x,
	const kVector<V>&	y)
{
	//	set x
	myX = x;
	myY = y;

	//	dims
	int numX = myX.size();
	int numY = myY.size();

	//	resize
	myR.resize(numX, numY, 0.0);
	myMuX.resize(numX, numY, 0.0);
	myVarX.resize(numX, numY, 0.0);
	myMuY.resize(numX, numY, 0.0);
	myVarY.resize(numX, numY, 0.0);
	myCovXY.resize(numX, numY, 0.0);

	//	calc x operators
	kFiniteDifference::dx(-1, myX, myDxd);
	kFiniteDifference::dx( 0, myX, myDx);
	kFiniteDifference::dx( 1, myX, myDxu);
	kFiniteDifference::dxx(   myX, myDxx);

	//	calc y operators
	kFiniteDifference::dx(-1, myY, myDyd);
	kFiniteDifference::dx( 0, myY, myDy);
	kFiniteDifference::dx( 1, myY, myDyu);
	kFiniteDifference::dxx(   myY, myDyy);

	//	resize results
	myRes.resize(numV);
	myVx.resize(numV);
	myVy.resize(numV);
	myVxy.resize(numV);
	for(int h = 0; h < numV; ++h)
	{
		myRes(h).resize(numX, numY);
		myVx(h).resize(numY, numX);
		myVy(h).resize(numX, numY);
		myVxy(h).resize(numY, numX);
	}

	//	ws
	myWs.resize(max(numX, numY));

	//	done
	return;
}

//	init from sigma rho
template <class V>
void	
kFd2d<V>::initFromSigmaRho(
	const kMatrix<V>&	sigmaX,
	const kMatrix<V>&	sigmaY,
	const kMatrix<V>&	corrXY)
{
	//	dims
	int numX = myX.size()-1;
	int numY = myY.size()-1;

	//	helps
	int i, j;
	for(i=0;i<=numX;++i)
	{
		for(j=0;j<=numY;++j)
		{
			myCovXY(i,j) = i==0 || i==numX || j==0 || j==numY ? 0.0 : sigmaX(i,j)*sigmaY(i,j)*corrXY(i,j);
			myVarX(i,j)  = sigmaX(i,j)*sigmaX(i,j);
			myVarY(i,j)  = sigmaY(i,j)*sigmaY(i,j);
		}
	}
	
	//	done
	return;
}

//	calc Ax
template <class V>
void	
kFd2d<V>::calcAx(
	int			j,
	V			one,
	V			dt,
	int			wind,
	bool		tr,
	kMatrix<V>&	Ax) const
{
	//	dims
	int numX = myX.size();
	int numC = myDx.cols();
	int mm   = numC/2;

	//	resize
	Ax.resize(numX, numC);

	//	operator
	const kMatrix<V>* Dx = &myDx;
	if	   (wind< 0)  Dx = &myDxd;
	else if(wind==1)  Dx = &myDxu;

	//	helps
	V mu, var;
	int i,k;

	//	loop
	for(i=0;i<numX;++i)
	{
		mu  = myMuX(i, j);
		var = myVarX(i, j);

		if(wind==2)
		{
			if(mu<0.0) Dx = &myDxd;
			else	   Dx = &myDxu;
		}
		for(k=0;k<numC;++k)
		{
			Ax(i, k) = dt*(mu*(*Dx)(i,k) + 0.5*var*myDxx(i,k));
		}
		Ax(i,mm) += 1.0 - dt*0.5*myR(i,j);
	}

	//	tr
	if(tr) kMatrixAlgebra::transpose(mm, Ax);

	//	done
	return;
}

//	calc Ay
template <class V>
void	
kFd2d<V>::calcAy(
	int			i,
	V			one,
	V			dt,
	int			wind,
	bool		tr,
	kMatrix<V>& Ay) const
{
	//	dims
	int numY = myY.size();
	int numC = myDy.cols();
	int nn   = numC / 2;

	//	resize
	Ay.resize(numY, numC);

	//	operator
	const kMatrix<V>* Dy = &myDy;
	if     (wind< 0)  Dy = &myDyd;
	else if(wind==1)  Dy = &myDyu;

	//	helps
	V mu, var;
	int j, k;

	//	loop
	for(j=0;j<numY;++j)
	{
		mu  = myMuY(i, j);
		var = myVarY(i, j);
		if(wind==2)
		{
			if(mu<0.0) Dy = &myDyd;
			else	   Dy = &myDyu;
		}
		for(k=0;k<numC;++k)
		{
			Ay(j,k) = dt*(mu*(*Dy)(j,k) + 0.5*var*myDxx(j,k));
		}
		Ay(j,nn) += 1.0 - dt * 0.5 * myR(i,j);
	}

	//	tr
	if(tr) kMatrixAlgebra::transpose(nn, Ay);

	//	done
	return;
}

//	calc Axy
template <class V>
void	
kFd2d<V>::calcAxy(
	int			i,
	int			j,
	V			one,
	V			dt,
	bool		tr,
	kMatrix<V>& Axy) const
{
	//	dims
	int numX  = myX.size()-1;
	int numY  = myY.size()-1;
	int numCx = myDx.cols();
	int numCy = myDy.cols();
	int mm    = numCx/2;
	int nn    = numCy/2;

	//	resize
	Axy.resize(numCx, numCy);

	//	loop bounds
	int ml = max(0, i - mm);
	int mu = min(i + mm, numX);
	int nl = max(0, j - nn);
	int nu = min(j + nn, numY);

	//	helps
	int k, l, m, n;

	//	fill
	if(!tr)
	{
		for(m=ml;m<=mu;++m)
		{
			k = m - i + mm;
			for(n=nl;n<=nu;++n)
			{
				l = n - j + nn;
				Axy(k,l) = dt*myCovXY(i,j)*myDx(i,k)*myDy(j,l);
			}
		}
		
	}
	else
	{
		for(m=ml;m<=mu;++m)
		{
			k = m-i+mm;
			for(n=nl;n<=nu;++n)
			{
				l = n-j+nn;
				Axy(k,l) = dt*myCovXY(m,n)*myDx(m,k)*myDy(n,l);
			}
		}
	}
	Axy(mm, nn) += one;

	//	done
	return;
}

//	mult Axy
template <class V>
V
kFd2d<V>::multAxy(
	int					i,
	int					j,
	const kMatrix<V>&	Axy,
	const kMatrix<V>&	b) const
{
	//	dims
	int numX = myX.size()-1;
	int numY = myY.size()-1;
	int mm   = myDx.cols()/2;
	int nn   = myDy.cols()/2;

	//	dims
	int ml = max(0, i - mm);
	int mu = min(i + mm, numX);
	int nl = max(0, j - nn);
	int nu = min(j + nn, numY);

	//	helps
	int k, l, m, n;

	//	loop
	V res = 0.0;
	for(m=ml;m<=mu;++m)
	{
		k = m - i + mm;
		for(n=nl;n<=nu;++n)
		{
			l = n - j + nn;
			res += Axy(k, l) * b(m, n);
		}
	}

	//	done
	return res;
}


//	roll bwd
template <class V>
void	
kFd2d<V>::rollBwd(
	V						dt,
	int						windX,
	int						windY,
	kVector<kMatrix<V> >&	res)
{
	//	dims 
	int numX = myX.size();
	int numY = myY.size();
	int numV = res.size();

	//	helps
	int h, i, j;

	//	explicit xy
	for(i=0;i<numX;++i)
	{
		for(j=0;j<numY;++j)
		{
			calcAxy(i, j, 1.0, 0.5*dt, false, myAxy);
			for(h=0;h<numV;++h)
			{
				myVxy(h)(j, i) = multAxy(i, j, myAxy, res(h));
			}
		}
	}

	//	implicit x
	for(j=0;j<numY;++j)
	{
		calcAx(j, 1.0, -dt, windX, false, myAx);
		for(h=0;h<numV;++h)
		{
			kMatrixAlgebra::tridag(myAx(), myVxy(h)(j), myVx(h)(j), myWs());
		}
	}

	//	transpose
	for(h=0;h<numV;++h)
	{
		kMatrixAlgebra::transpose(myVx(h), res(h));
	}

	//	explicit xy
	for(i=0;i<numX;++i)
	{
		for(j=0;j<numY;++j)
		{
			calcAxy(i, j, 1.0, 0.5*dt, false, myAxy);
			for(h=0;h<numV;++h)
			{
				myVy(h)(i,j) = multAxy(i, j, myAxy, res(h));
			}
		}
	}

	//	implicit y 
	for(i=0;i<numX;++i)
	{
		calcAy(i,1.0,-dt,windY,false,myAy);
		for(h=0;h<numV;++h)
		{
			kMatrixAlgebra::tridag(myAy(), myVy(h)(i), res(h)(i), myWs());
		}
	}

	//	done
	return;
}

//	roll fwd
template <class V>
void	
kFd2d<V>::rollFwd(
	V						dt,
	int						windX,
	int						windY,
	kVector<kMatrix<V> >&	res)
{
	//	dims
	int numX = myX.size();
	int numY = myY.size();
	int numV = res.size();

	//	helps
	int h, i, j;

	//	implicit y
	for(i=0;i<numX;++i)
	{
		calcAy(i, 1.0, -dt, windY, true, myAy);
		for(h=0;h<numV;++h)
		{
			kMatrixAlgebra::tridag(myAy(), res(h)(i), myVy(h)(i), myWs());
		}
	}

	//	explicit xy
	for(i=0;i<numX;++i)
	{
		for(j=0;j<numY;++j)
		{
			calcAxy(i, j, 1.0, 0.5*dt, true, myAxy);
			for(h=0;h<numV;++h)
			{
				myVxy(h)(j,i) = multAxy(i, j, myAxy, myVy(h));
			}
		}
	}

	//	implicit x
	for(j=0;j<numY;++j)
	{
		calcAx(j, 1.0, -dt, windX, true, myAx);
		for(h=0;h<numV;++h)
		{
			kMatrixAlgebra::tridag(myAx(), myVxy(h)(j), myVx(h)(j), myWs());
		}
	}

	//	transpose
	for(h=0;h<numV;++h)
	{
		kMatrixAlgebra::transpose(myVx(h), myVy(h));
	}

	//	explicit xy
	for(i=0;i<numX;++i)
	{
		for(j=0;j<numY;++j)
		{
			calcAxy(i, j, 1.0, 0.5*dt, true, myAxy);
			for(h=0;h<numV;++h)
			{
				res(h)(i,j) = multAxy(i, j, myAxy, myVy(h));
			}
		}
	}

	//	done
	return;
}

