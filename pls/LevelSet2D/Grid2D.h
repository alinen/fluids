/**************************************************************************
	ORIGINAL AUTHOR: 
		Emud Mokhberi (emud@ucla.edu)
	MODIFIED BY:
	
	CONTRIBUTORS:

-----------------------------------------------
	
 ***************************************************************
 ******General License Agreement and Lack of Warranty ***********
 ****************************************************************

 This software is distributed for noncommercial use in the hope that it will 
 be useful but WITHOUT ANY WARRANTY. The author(s) do not accept responsibility
 to anyone for the consequences of using it or for whether it serves any 
 particular purpose or works at all. No guarantee is made about the software 
 or its performance.

 You are allowed to modify the source code, add your name to the
 appropriate list above and distribute the code as long as 
 this license agreement is distributed with the code and is included at 
 the top of all header (.h) files.

 Commercial use is strictly prohibited.
***************************************************************************/


/*
	Grid2D : A class for representing and working with a 2D grid of points

	The grid is created with a 1 cell buffer on each side. However, note that these
	buffer cells are not hidden from the user. That means that indexing the grid
	with for instance (0,0) will return a buffer cell and not the first non-buffer cell.
	
	The grid can be initialized based on an input array that does not contain the 1 cell
	buffer. This can be done when the grid is constructed or by calling the "set" function.
	Note that the 1 cell boundary has to be set after this is done to ensure proper boundary
	conditions.
	
	Boundary functions are provided to set the boundary cells to certain values depending
	on the use of the grid. The current boundary functions are:
	Dirichlet - boundary cells are set to 0
	Neumann - boundary cells are set to the value of the adjacent cell
	Velocity U - Left and Right boundarys are set to the negative value of their neighbors
				 Top and Bottom are set as Neumann
	Velocity V - Left and Right boundarys are set as Neumann
				 Top and Bottom are set to the negative value of their neighbors
	SignedDistance - This is currently set to 3 times the cell size with the only purpose
					 being to make sure that the boundary is always considered outside of
					 of the implicit surface

	Created by Emud Mokhberi: UCLA : 09/04/04
*/

#ifndef GRID2D_H
#define GRID2D_H
#include "main.h"


class Grid2D
{
private:
	inline int GI(int i, int j) const { return i+(Nx+2)*j; }

	int Nx, Ny, size;
	Double *grid;

public:
	Grid2D(int nx, int ny) : Nx(nx), Ny(ny) 
		{ size = (nx+2) * (ny+2); grid = new Double[size]; fill(grid, grid+size, 0.); }
	Grid2D(const Grid2D &gi) 
		{ gi.GetSize(Nx, Ny, size); if(!size) grid = NULL;
		  else { grid = new Double[size]; FOR_GRID2D grid[i] = gi[i];} }
	Grid2D(int nx, int ny, const Double val[]) : Nx(nx), Ny(ny)
		{ size = (nx+2) * (ny+2); grid = new Double[size]; 
		  for(int j=1; j<=ny; j++) for(int i=1; i<=nx; i++) grid[GI(i,j)] = val[(j-1)*nx + (i-1)]; }
	~Grid2D() { if(grid != NULL) delete [] grid; }
			
	inline Double& operator[] (int index) { return grid[index]; }
	inline const Double& operator[] (int index) const { return grid[index]; }
	inline Double& operator() (int i, int j) { return grid[GI(i,j)]; }
	inline const Double& operator() (int i, int j) const { return grid[GI(i,j)]; }
	
	operator Double*() { return &grid[0]; }
	operator const Double*() { return &grid[0]; }

	inline Grid2D& operator=(const Grid2D &gi)  { FOR_GRID2D grid[i] = gi[i]; return *this; }
	inline Grid2D& operator+=(const Grid2D &gi) { FOR_GRID2D grid[i] += gi[i]; return *this;}
	inline Grid2D& operator-=(const Grid2D &gi) { FOR_GRID2D grid[i] -= gi[i]; return *this;}
	inline Grid2D& operator*=(Double c)		 { FOR_GRID2D grid[i] *= c; return *this; }
	inline Grid2D& operator/=(Double c) 
		{ Double cInv = 1. / c; FOR_GRID2D grid[i] *= cInv; return *this; }

	inline void set(const Double val[])
		{ for(int j=1; j<=Ny; j++) for(int i=1; i<=Nx; i++) grid[GI(i,j)] = val[(j-1)*Nx + (i-1)]; }
	inline Double dot(const Grid2D& gi) const
		{ Double ret = 0; FOR_GRID2D ret += grid[i] * gi[i]; return ret; }
	Double SquaredLength() const { return (*this).dot(*this); }
	Double Length() const { return sqrt(SquaredLength()); }
	void Normalize() { (*this) /= Length();}
	inline void GetSize(int &nx, int &ny, int &s) const { nx = Nx; ny = Ny; s = size; }
	inline void Clear() { fill(grid, grid+size, 0.); }

	//b = -1 -> Dirichlet
	//b = 0  -> Neumann
	//b = 1  -> U Velocities
	//b = 2  -> V Velocities
	inline void SetBoundary(int b)
	{
		if(b == 0)		SetBoundaryNeumann();
		else if(b == 1)	SetBoundaryU();
		else if(b == 2)	SetBoundaryV();
		else 			SetBoundaryDirichlet();
	}
	inline void SetBoundaryDirichlet();
	inline void SetBoundaryNeumann();
	inline void SetBoundaryU();
	inline void SetBoundaryV();
	inline void SetBoundaryCorner();	
	inline void SetBoundarySignedDist();
};

inline void Grid2D::SetBoundaryDirichlet()
{
	for(int i=1 ; i<=Ny; i++ ) grid[GI(0   ,i)] = grid[GI(Nx+1,i)] = 0;
	for(int i=1 ; i<=Nx; i++ ) grid[GI(i,0   )] = grid[GI(i,Ny+1)] = 0;
	grid[GI(0   ,0   )] = grid[GI(0   ,Ny+1)] = 0;
	grid[GI(Nx+1,0   )] = grid[GI(Nx+1,Ny+1)] = 0;
}

inline void Grid2D::SetBoundaryNeumann()
{
	for(int i=1 ; i<=Ny ; i++ ) {
		grid[GI(0   ,i)] = grid[GI(1 , i)];
		grid[GI(Nx+1,i)] = grid[GI(Nx, i)];
	}
	for(int i=1 ; i<=Nx; i++ ) { 
		grid[GI(i,0   )] = grid[GI(i , 1)];
		grid[GI(i,Ny+1)] = grid[GI(i ,Ny)];
	}
	SetBoundaryCorner();
}

inline void Grid2D::SetBoundaryU()
{
	for(int i=1 ; i<=Ny ; i++ ) {
		grid[GI(0   ,i)] = -grid[GI( 1, i)];
		grid[GI(Nx+1,i)] = -grid[GI(Nx, i)];
	}
	for(int i=1 ; i<=Nx; i++ ) { 
		grid[GI(i,0   )] = grid[GI(i , 1)];
		grid[GI(i,Ny+1)] = grid[GI(i ,Ny)];
	}
	SetBoundaryCorner();
}

inline void Grid2D::SetBoundaryV()
{
	for(int i=1 ; i<=Ny ; i++ ) {
		grid[GI(0   ,i)] = grid[GI(1 , i)];
		grid[GI(Nx+1,i)] = grid[GI(Nx, i)];
	}
	for(int i=1 ; i<=Nx; i++ ) { 
		grid[GI(i,0   )] = -grid[GI( i, 1)];
		grid[GI(i,Ny+1)] = -grid[GI( i,Ny)];
	}
	SetBoundaryCorner();
}

inline void Grid2D::SetBoundaryCorner()
{
	grid[GI(0   ,0   )] = 0.5*(grid[GI(1 ,0   )]+grid[GI(0   ,1 )]);
	grid[GI(0   ,Ny+1)] = 0.5*(grid[GI(1 ,Ny+1)]+grid[GI(0   ,Ny)]);
	grid[GI(Nx+1,0   )] = 0.5*(grid[GI(Nx,0   )]+grid[GI(Nx+1,1 )]);
	grid[GI(Nx+1,Ny+1)] = 0.5*(grid[GI(Nx,Ny+1)]+grid[GI(Nx+1,Ny)]);
}

inline void Grid2D::SetBoundarySignedDist()
{
    for(int i=1 ; i<=Ny ; i++ ) {
		grid[GI(0   ,i)] = 3*HH;
		grid[GI(Nx+1,i)] = 3*HH;
	}
	for(int i=1 ; i<=Nx; i++ ) { 
		grid[GI(i,0   )] = 3*HH;
		grid[GI(i,Ny+1)] = 3*HH;
    }
	SetBoundaryCorner();
}

inline ostream& operator<<(ostream& out, const Grid2D& g)
{	
	int nx, ny, size; 
	g.GetSize(nx, ny, size);
	size = 0;
	for(int j=0; j<ny+2; j++) 
	{ 
		for(int i=0; i<nx+2; i++)
		{
			out << g[size] << " ";
			size++;
		}
		out << endl;
	}		
	return out;
}

inline istream& operator>>(istream& in, Grid2D& g)
{	
	int nx, ny, size; 
	g.GetSize(nx, ny, size);
	for(int i=0; i<size; i++) in >> g[i]; 
	return in; 
}

#endif