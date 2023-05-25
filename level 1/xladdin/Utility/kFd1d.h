#pragma once

//	desc:	1d finite difference solution for pdes of the form
//
//		0 = dV/dt + A V
//
//		A = -r + mu d/dx + 1/2 var d^2/dx^2
//
//	using the theta scheme
//
//		[1-theta dt A] V(t) = [1 + (1-theta) dt A] V(t+dt)
//

//	includes
#include "kFiniteDifference.h"
#include "kMatrixAlgebra.h"

//	class declaration
template <class V>
class kFd1d
{
public:

	//	init 
	void	init(
		int						numV,
		const kVector<V>&		x,
		bool					log);

	//	operator
	void	calcAx(
		V						one,
		V						dtTheta,
		int						wind,
		bool					tr,
		kMatrix<V>&				A) const;

	//	roll bwd
	void	rollBwd(
		V						dt,
		bool					update,
		V						theta,
		int						wind,
		kVector<kVector<V>>&	res);

	//	roll fwd
	void	rollFwd(
		V						dt,
		bool					update,
		V						theta,
		int						wind,
		kVector<kVector<V>>&	res);

	//	x, r, mu, var
	kVector<V>	myX, myR, myMu, myVar;

	//	vector of results
	kVector<kVector<V>> myRes;

	//	private parts
private:

	//	diff operators
	kMatrix<V>	myDxd, myDxu, myDx, myDxx;

	//	operator matrix
	kMatrix<V>	myAe, myAi;

	//	helper
	kVector<V>	myVs, myWs;

	
};

//	init
template <class V>
void
kFd1d<V>::init(
	int					numV,
	const kVector<V>&	x,
	bool				log)
{
	myX = x;
	myRes.resize(numV);
	for(int k=0;k<numV;++k) myRes[k].resize(myX.size());

	//	resize params
	myR.resize(myX.size(), 0.0);
	myMu.resize(myX.size(), 0.0);
	myVar.resize(myX.size(), 0.0);

	kFiniteDifference::dx(-1,myX,myDxd);
	kFiniteDifference::dx( 0,myX,myDx);
	kFiniteDifference::dx( 1,myX,myDxu);
	kFiniteDifference::dxx(  myX,myDxx);

	if (myX.empty()) return;

	int numC = myDx.cols();

	//	log transform case
	if(log)
	{
		int n = myX.size()-1;
		for(int i = 1;i<n;++i)
		{
			for(int j=0;j<numC;++j)
				myDxx(i,j) -= myDx(i,j);
		}
	}

	myAe.resize(myX.size(),numC);
	myAi.resize(myX.size(),numC);
	myVs.resize(myX.size());
	myWs.resize(myX.size());

	//	done
	return;
}

//	construct operator
template <class V>
void
kFd1d<V>::calcAx(
	V				one,
	V				dtTheta,
	int				wind,
	bool			tr,
	kMatrix<V>&		A) const
{
	// todo: implement
}

//	roll bwd
template <class V>
void
kFd1d<V>::rollBwd(
	V						dt,
	bool					update,
	V						theta,
	int						wind,
	kVector<kVector<V>>&	res)
{
	//	helps
	int k;

	//	dims
	int n = myX.size();
	int mm = 1;
	int numV = (int)res.size();

	//	explicit
	if (theta != 1.0)
	{
		if (update) calcAx(1.0, dt * (1.0 - theta), wind, false, myAe);
		for (k = 0; k < numV; ++k)
		{
			myVs = res[k];
			kMatrixAlgebra::banmul(myAe, mm, mm, myVs, res[k]);
		}
	}

	//	implicit
	if (theta != 0.0)
	{
		if (update) calcAx(1.0, -dt * theta, wind, false, myAi);
		for (k = 0; k < numV; ++k)
		{
			myVs = res[k];
			kMatrixAlgebra::tridag(myAi, myVs, res[k], myWs);
		}
	}

	//	done
	return;
}

//	roll fwd
template <class V>
void
kFd1d<V>::rollFwd(
	V						dt,
	bool					update,
	V						theta,
	int						wind,
	kVector<kVector<V>>&	res)
{
	// todo: implement
}
