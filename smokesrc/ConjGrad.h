// ==========================================================================
// Copyright (C) 2007 Ben Sunshine-Hill
// ==========================================================================

#ifndef CONJGRAD_H
#define CONJGRAD_H

#pragma warning(disable: 4244 4267 4996)

#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/matrix.hpp"
#include <boost/numeric/ublas/io.hpp>

#include <utility>


/** Use the conjugate gradient method to solve a symmetric, positive semidefinite system of linear equations.
  *
  * This is an iterative method which attempts to converge on a solution.
  * The method will exit when either the number of iterations reaches max_iter,
  * or the 2-norm of the residual goes below tol.
  *
  * @param A         The system of linear equations to solve. This must be symmetric and positive definite.
  * @param b         The constant terms of the system.
  * @param x         Where the output is stored. This must be a vector of the appropriate size. It is overwritten by this method.
  * @param max_iter  The maximum number of iterations the solver will make.
  * @param tol       The maximum tolerable error, below which the solver will consider the solution to have converged.
  */
template<typename Matrix, typename Vector>
bool cg_solve(const Matrix &A, const Vector &b, Vector &x, int max_iter, double tol)
{
	std::fill(x.begin(), x.end(), 0);

	Vector r = b;
	Vector p = b;
	Vector r_old;
	Vector temp;

	// CG loop
	for(int niter = 0; niter < max_iter; niter++)
	{
		temp = prod(A,p);
		double alpha = inner_prod(r,r)/inner_prod(p,temp);
		x += (p*alpha);
		r_old = r;
		r -= (temp*alpha);
		double residn = norm_2(r);
		if(residn < tol){
            std::cout << "numiters: "<< niter <<std::endl;
			return true;
		}
		double beta = inner_prod(r,r)/inner_prod(r_old,r_old);      
		p = r + p*beta;
	}

   std::cout << "WARNING: cg_solve did not converge" << std::endl;
	return false;
}

int getIndex(int i, int j, int k)
{
   extern int theDim[3]; // Naughty globals...

   if (i < 0 || i > theDim[0]-1) return -1;
   if (j < 0 || j > theDim[1]-1) return -1;
   if (k < 0 || k > theDim[2]-1) return -1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];
   return col+row+stack;
}

void getCell(int index, int& i, int& j, int& k)
{
   extern int theDim[3]; // Naughty globals...
   j = (int) index/(theDim[0]*theDim[2]);           // stack
   k = (int) (index - j*theDim[0]*theDim[2])/theDim[0];       // row
   i = index - j*theDim[0]*theDim[2] - k*theDim[0]; // col
}

template<typename Matrix, typename Vector>
void applyPreconditioner(const Matrix &A, const Vector &precon, const Vector &r, Vector& z)
{
    unsigned int numCells = r.size();
    ublas::vector<double> q(numCells);
    std::fill(q.begin(), q.end(), 0);

    // First solve Lq = r
    for (unsigned int index = 0; index < numCells; index++)
    {
        //std::cout << "Initializing row " << row << "/" << numCells << std::endl;
        int i, j, k; getCell(index, i, j, k); // Each row corresponds to a cell

        int neighbori = getIndex(i-1,j,k);
        int neighborj = getIndex(i,j-1,k);
        int neighbork = getIndex(i,j,k-1);
        double termi = neighbori != -1? A(index,neighbori) * precon(neighbori) * q(neighbori) : 0;
        double termj = neighborj != -1? A(index,neighborj) * precon(neighborj) * q(neighborj) : 0;
        double termk = neighbork != -1? A(index,neighbork) * precon(neighbork) * q(neighbork) : 0;
        
        double t = r(index) - termi - termj - termk;
        q(index) = t * precon(index);
    }

    //std::cout << "q: " << q << std::endl;

    // Next solve transpose(L)z = q
    for (int index = numCells-1; index >= 0; index--)
    {
        //std::cout << "Initializing row " << row << "/" << numCells << std::endl;
        int i, j, k; getCell(index, i, j, k); // Each row corresponds to a cell

        int neighbori = getIndex(i+1,j,k);
        int neighborj = getIndex(i,j+1,k);
        int neighbork = getIndex(i,j,k+1);
        double termi = neighbori != -1? A(index,neighbori) * precon(index) * z(neighbori) : 0;
        double termj = neighborj != -1? A(index,neighborj) * precon(index) * z(neighborj) : 0;
        double termk = neighbork != -1? A(index,neighbork) * precon(index) * z(neighbork) : 0;
        
        double t = q(index) - termi - termj - termk;
        z(index) = t * precon(index);
    }
    //std::cout << "z: " << z << std::endl;
}

template<typename Matrix, typename Vector>
bool cg_psolve(const Matrix &A, const Vector &precon, const Vector &b, Vector &p, int max_iter, double tol)
{
   std::fill(p.begin(), p.end(), 0);
   Vector r = b;
   Vector z = b;
   Vector s = b;
   //std::cout << "r: " << r << std::endl;
   applyPreconditioner(A, precon, r, s); 

   double resign;
   double sigma = inner_prod(s, r);

	for(int niter = 0; niter < max_iter; niter++)
	{
      z = prod(A, s);
      double alpha = sigma/inner_prod(s, z);
      p += alpha*s;
      r -= alpha*z;
      resign = norm_2(r);
      if(resign < tol) 
      { 
        //std::cout << "numiters: "<< niter <<std::endl;
        return true; 
      }
      
      applyPreconditioner(A, precon, r, z);
      double sigma_new = inner_prod(z, r);
      double beta = sigma_new/sigma;
      s = z + beta*s;
      sigma = sigma_new;
   }

   std::cout << "WARNING: cg_psolve did not converge: " << resign << std::endl;
   return false;
}

#endif