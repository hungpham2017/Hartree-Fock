#include <fstream>
#include <vector>
#include <iterator>
#include <iostream>
#include <lawrap/blas.h>
#include <lawrap/lapack.h>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>

//Copy a vector X to a vector Y:
std::vector<double> copyvec(std::vector<double> X)
{
	return X;
}

//Print a vector X
void printvec(std::vector<double> X)
{
	for (int i = 0; i < X.size(); i++)
	{
		std::cout << X[i] << "\n";
	}
}

//Add two vectors -> X + Y
std::vector<double> addvec(std::vector<double> X, std::vector<double> Y)
{
	LAWrap::axpy(X.size(), 1.0, X.data(), 1, Y.data(), 1);
	return Y;
}

//Add two vectors with coefficients -> aX + bY
std::vector<double> addvec2(double a, std::vector<double> X, double b, std::vector<double> Y)
{
	int size = X.size();
	std::vector<double> temp(size);
	for (int i = 0; i < size; i++)
	{
		temp[i] = a*X[i] + b*Y[i];
	}
	 
	return temp;	
}
	
//Subtract two vectors -> X - Y
std::vector<double> subvec(std::vector<double> X, std::vector<double> Y)
{
	LAWrap::axpy(X.size(), -1.0, Y.data(), 1, X.data(), 1);
	return Y;
}

//Element-wise sum of a 4-tensor X
double sumten(int nao, std::vector<std::vector<std::vector<std::vector<double>>>> X)
{
	double sum = 0;
	for (int p = 0; p < nao; p++)
	{
		for (int q = 0; q < nao; q++)
		{
			for (int r = 0; r < nao; r++)
			{
				for (int s = 0; s < nao; s++)
				{
					sum = sum + X[p][q][r][s];
				}
			}
		}
	}
	return sum;
}

//Dot product of two square vectors with the same size -> X*Y
std::vector<double> dotvec(std::vector<double> X, std::vector<double> Y)
{
	int size = sqrt(X.size());
	std::vector<double> product(X.size());
	LAWrap::gemm('N', 'N', size, size, size, 1.0, X.data(), size, Y.data(), size, 0.0, product.data(), size);
	return product;
}

//READ a vector and turns it into a matrix
std::vector<std::vector<double>> readmax(std::string file, int nao)
{
	std::ifstream Cdata(file);
	std::vector<double> Cvec;	
		std::copy(
		std::istream_iterator<double>(Cdata), 
		std::istream_iterator<double>(), 
		std::back_inserter(Cvec));
	std::vector<std::vector<double>> Cmax;
	Cmax.resize(nao, std::vector<double> (nao));
	for (int i = 0; i < nao; i++)
	{
		Cmax[i].resize ( nao );
		for (int j = 0; j < nao; j++)
		{
			Cmax[i][j] = Cvec[i*nao + j];
		}
	}
	return Cmax;
}

//READ a vector and turns it into a 4-rank tensor
std::vector<std::vector<std::vector<std::vector<double>>>> readten(std::string file, int nao)
{
	std::ifstream Cdata(file);
	std::vector<double> Cvec;	
		std::copy(
		std::istream_iterator<double>(Cdata), 
		std::istream_iterator<double>(), 
		std::back_inserter(Cvec));
	std::vector<std::vector<std::vector<std::vector<double>>>> Cten;
	Cten.resize(nao, std::vector<std::vector<std::vector<double>>> (nao, std::vector<std::vector<double>> (nao, std::vector<double> (nao))));
	for (int i = 0; i < nao; i++)
	{
		for (int j = 0; j < nao; j++)
		{
			for (int k = 0; k < nao; k++)
			{
				for (int l = 0; l < nao; l++)
				{
					Cten[i][j][k][l] = Cvec[i*pow(nao,3) + j*pow(nao,2) + k*nao + l];
				}
			}	
		}
	}
	return Cten;
}

//READ a vector
std::vector<double> readvec(std::string file)
{
	std::ifstream data(file);
	std::vector<double> vec;	
		std::copy(
		std::istream_iterator<double>(data), 
		std::istream_iterator<double>(), 
		std::back_inserter(vec));
	return vec;
}

//Diagonalize a Fock in the orthogonal basis
std::vector<std::vector<double>> diag(std::vector<double> F, std::vector<double> UsU)
{
	//Compute the UsU.T @ F @ UsU
	int nao = sqrt(F.size());
	std::vector<double> temp1(nao);		  // to store  eigen values
	std::vector<double> temp2(nao*nao);	  // to store  @ UsU
	std::vector<double> FT(nao*nao);
	std::vector<double> eigval(nao*nao, 0);	
	LAWrap::gemm('N', 'N', nao, nao, nao, 1.0, F.data(), nao, UsU.data(), nao, 0.0, temp2.data(), nao);
	LAWrap::gemm('T', 'N', nao, nao, nao, 1.0, UsU.data(), nao, temp2.data(), nao, 0.0, FT.data(), nao);
	
	//Diagonalize FT (Fock in the orthogonal basis)
	LAWrap::heev('V', 'U',  nao, FT.data(), nao, temp1.data());
	
	std::vector<double> C(nao*nao);	
	//Retrun C~ to C, C here is not in the same order in the C.ravel() in python
	LAWrap::gemm('N', 'N', nao, nao, nao, 1.0, UsU.data(), nao, FT.data(), nao, 0.0, C.data(), nao);	
	
	
	for (int i = 0; i < nao; i++)
	{
		eigval[i*nao + i] = temp1[i];
	}
	std::vector<std::vector<double>> result;
	result.resize(2, std::vector<double> (nao*nao));
	
	for (int i = 0; i < nao*nao; i++)
	{
			result[0][i] = eigval[i];
	}

	for (int i = 0; i < nao*nao; i++)
	{
			result[1][i] = C[i];
	}
	
	return result;	
}

//Compute VHF
std::vector<double> VHFcal(std::vector<std::vector<std::vector<std::vector<double>>>> g, std::vector<double> D)
{
	//Compute VHF = 2*J - K
	int nao = sqrt(D.size());
	std::vector<double> VHF(nao*nao, 0);
	double test = 0;
	for (int p = 0; p < nao; p++)
	{
		for (int q = 0; q < nao; q++)
		{
			for (int r = 0; r < nao; r++)
			{
				for (int s = 0; s < nao; s++)
				{
					VHF[p*nao + q] = VHF[p*nao + q] + D[r*nao +s]*( g[p][q][r][s] - 0.5 * g[p][r][q][s]);
				}
			}
		}
	}
	
	return VHF;
}

//Construct D from C*C.T
std::vector<double>  RDM1cal(int nocc, std::vector<double> C)
{
	
	int nao = sqrt(C.size());
	std::vector<double> D(C.size());
	LAWrap::gemm('N', 'T', nao, nao, nocc, 2.0, C.data(), nao, C.data(), nao, 0.0, D.data(), nao);
	return D;
}

//Compute E_elec from H, F, D
double E_elec(std::vector<double> H, std::vector<double> F, std::vector<double> D)
{
	std::vector<double> HplusF(H.size());
	double E = 0;
	HplusF = addvec(H, F);
	for (int i = 0; i < H.size(); i++)
	{
		E = E + 0.5 * HplusF[i] * D[i];
	}
}

int main()
{
	//SCF parameters
	int nao = 7;
	int nocc = 5;
	int ncyc = 30;	
	double e_conv = 1.e-10;
	double d_conv = 1.e-10;
	double damp_value = 0.20;
	int damp_start = 5;

	//READ g vector	and turn it to 4-rank tensor
	std::vector<std::vector<std::vector<std::vector<double>>>> g;
	g = readten("g.data", nao);
	
	//READ S vector	
	std::vector<double> Svec;
	Svec = readvec("S.data");
	
	//READ H vector	
	std::vector<double> Hvec;
	Hvec = readvec("H.data");

	//READ C vector	
	std::vector<double> Cvec;
	Cvec = readvec("eivec.data");
	
	//Calculate tranformation matrix s = s^-1/2, UsU = S^-1/2 = Us^-1/2U*
	std::vector<double> W(nao);
	std::vector<double> U(nao*nao, 0.0);
	std::vector<double> s(nao*nao);	
	std::vector<double> Us(nao*nao);
	std::vector<double> UsU(nao*nao);	

	U = copyvec(Svec);
	LAWrap::heev('V', 'U',  nao, U.data(), nao, W.data()); //After this step, U is the eigenvector matrix
	for (int i = 0; i < nao; i++)
	{
		s[i*nao + i] = 1.0/sqrt(W[i]);
	}
	LAWrap::gemm('N', 'N', nao, nao, nao, 1.0, U.data(), nao, s.data(), nao, 1.0, Us.data(), nao);
	LAWrap::gemm('N', 'T', nao, nao, nao, 1.0, Us.data(), nao, U.data(), nao, 1.0, UsU.data(), nao);
	

	//Diagonalize the core H to get initial C, then D
	std::vector<std::vector<double>> result;
	result.resize(2, std::vector<double> (nao*nao));	
	result = diag(Hvec, UsU);  //result[0][...] for eigenvalues, result[1][...] for eigenvectors
	
	std::vector<double> C(nao*nao);	
	std::vector<double> D(nao*nao);
	C = copyvec(result[1]);
	D = Dcal(nocc, C);	
	
	std::vector<double> Fold(nao*nao, 0);
	std::vector<double> VHF(nao*nao);
	std::vector<double> Fock(nao*nao);
	std::vector<double> Fock_old(nao*nao);
	std::vector<double> Fock_new(nao*nao);
	
	double Eold = 0.0;
	double Ediff = 0.0;
	double Etotal = 0.0;	
	double Eelec = 0.0;
	double Enuc = 0.0;	
	
	//SCF ITERATION
	std::cout << "START - SCF:" << "\n";	
	for (int i = 0; i < ncyc; i++)
	{

		VHF = VHFcal(g, D);
		Fock_new = addvec(Hvec, VHF);
		
		//Damping
		if (i >= damp_start)
		{ 	
			Fock = addvec2(damp_value, Fock_old, 1.0 - damp_value, Fock_new); 
		}
		else
		{Fock = Fock_new;}
		
		Fock_old = Fock_new;
		
		Eelec = E_elec(Hvec, Fock, D);
		Etotal = Eelec + Enuc;
		Ediff = Etotal - Eold;
		Eold = Etotal;

		if (i >= damp_start)
		{ 	
			std::cout << "Damping SCF cycle: " << i << " E_total " << Etotal << " E_diff " <<  Ediff << "\n";
		}
		else
		{		
			std::cout << "SCF cycle: " << i << " E_total " << Etotal << " E_diff " <<  Ediff << "\n";
		}
		
		if ( std::abs(Ediff) < e_conv)
		{
			break;
		}
		
		diag(Fock, UsU);
		std::vector<std::vector<double>> diag_result;
		diag_result.resize(2, std::vector<double> (nao*nao));
		diag_result = diag(Fock, UsU);
		C = copyvec(diag_result[1]);
		D = RDM1cal(nocc, C);			
	}
	
	
	//printvec(D);
	//std::cout << g[0][2][6][0] << "\n";
	


	
	

	return 0;
}
