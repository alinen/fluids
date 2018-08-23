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
	Grid : A class for representing and working with a 3D grid of points

	The grid is created with a 1 cell buffer on each side. However, note that these
	buffer cells are not hidden from the user. That means that indexing the grid
	with for instance (0,0,0) will return a buffer cell and not the first non-buffer cell.
	
	The grid can be initialized based on an input array that does not contain the 1 cell
	buffer. This can be done when the grid is constructed or by calling the "set" function.
	Note that the 1 cell boundary has to be set after this is done to ensure proper boundary
	conditions.
	
	
	Boundary functions are provided to set the boundary cells to certain values depending
	on the use of the grid. The current boundary functions are:
	Dirichlet - boundary cells are set to 0
	Neumann - boundary cells are set to the value of the adjacent cell
	Velocity U - x-axis boundaries are set to the negative value of their neighbors
				 y-axis and z-axis boundaries are set as Neumann
	Velocity V - y-axis boundaries are set to the negative value of their neighbors
				 x-axis and z-axis boundaries are set as Neumann
	Velocity W - z-axis boundaries are set to the negative value of their neighbors
				 x-axis and y-axis boundaries are set as Neumann	
	SignedDistance - This is currently set to 3 times the cell size with the only purpose
					 being to make sure that the boundary is always considered outside of
					 of the implicit surface

	Created by Emud Mokhberi: UCLA : 09/04/04
*/

#ifndef GRID_H
#define GRID_H
#include "main.h"

class Grid
{
private:
	inline int GI(int i, int j, int k) const { return i + dj*j + dk*k; }

	int Nx, Ny, Nz, size, dj, dk;
	Double *grid;

public:
	Grid(int nx, int ny, int nz) : Nx(nx), Ny(ny), Nz(nz), dj(Nx+2), dk((Nx+2)*(Ny+2))
		{ size = (nx+2) * (ny+2) * (nz+2); grid = new Double[size]; fill(grid, grid+size, 0.); }
    Grid(int nx, int ny, int nz, Float c) : Nx(nx), Ny(ny), Nz(nz), dj(Nx+2), dk((Nx+2)*(Ny+2))
		{ size = (nx+2) * (ny+2) * (nz+2); grid = new Double[size]; fill(grid, grid+size, c); }
	Grid(const Grid &gi) 
		{ gi.GetSize(Nx, Ny, Nz, size); dj = Nx+2; dk = (Nx+2) * (Ny+2);
          grid = new Double[size]; FOR_GRID grid[i] = gi[i]; }
	Grid(int nx, int ny, int nz, const Double val[]) : Nx(nx), Ny(ny), Nz(nz)
		{ size = (nx+2) * (ny+2) * (nz+2); grid = new Double[size]; 
		  for (int k=1; k<=nz; k++) for(int j=1; j<=ny; j++) for(int i=1; i<=nx; i++) 
		  grid[GI(i,j,k)] = val[(k-1)*ny*nx + (j-1)*nx + (i-1)]; }
	~Grid() { if(grid != NULL) delete [] grid; }
			
	inline Double& operator[] (int index) { return grid[index]; }
	inline const Double& operator[] (int index) const { return grid[index]; }
	inline Double& operator() (int i, int j, int k) { return grid[GI(i,j,k)]; }
	inline const Double& operator() (int i, int j, int k) const { return grid[GI(i,j,k)]; }
	
	operator Double*() { return &grid[0]; }
	operator const Double*() { return &grid[0]; }

	inline Grid& operator=(const Grid &gi)  { FOR_GRID grid[i] = gi[i]; return *this; }
	inline Grid& operator+=(const Grid &gi) { FOR_GRID grid[i] += gi[i]; return *this;}
	inline Grid& operator-=(const Grid &gi) { FOR_GRID grid[i] -= gi[i]; return *this;}
	inline Grid& operator*=(Double c)		 { FOR_GRID grid[i] *= c; return *this; }
	inline Grid& operator/=(Double c) 
		{ Double cInv = 1. / c; FOR_GRID grid[i] *= cInv; return *this; }

    inline int GetNx() const { return Nx; }
    inline int GetNy() const { return Ny; }
    inline int GetNz() const { return Nz; }

	inline void set(const Double val[])
		{ for (int k=1; k<=Nz; k++) for(int j=1; j<=Ny; j++) for(int i=1; i<=Nx; i++) 
		  grid[GI(i,j,k)] = val[(k-1)*Ny*Nx + (j-1)*Nx + (i-1)]; }
	inline Double dot(const Grid& gi) const
		{ Double ret = 0; FOR_GRID ret += grid[i] * gi[i]; return ret; }
	inline Double normSqrd() const { return (*this).dot(*this); }
	inline Double norm() const { return sqrt(normSqrd()); }
	void Normalize() { (*this) /= norm();}
    inline void Clear() { fill(grid, grid+size, 0.); }
	inline void GetSize(int &nx, int &ny, int &nz, int &s) const 
        { nx = Nx; ny = Ny; nz = Nz; s = size; }

	//b = -1 -> Dirichlet
	//b = 0  -> Neumann
	//b = 1  -> U Velocities
	//b = 2  -> V Velocities
	//b = 3  -> W Velocities
	inline void SetBoundary(int b)
	{
		if(b == 0)		SetBoundaryNeumann();
		else if(b == 1)	SetBoundaryU();
		else if(b == 2)	SetBoundaryV();
		else if(b == 3) SetBoundaryW();
		else 			SetBoundaryDirichlet();
	}
	inline void SetBoundaryDirichlet();
	inline void SetBoundaryNeumann();
	inline void SetBoundaryU();
	inline void SetBoundaryV();
    inline void SetBoundaryW();
	inline void SetBoundaryCorner();	
    inline void SetBoundarySignedDist();
};

inline void Grid::SetBoundaryDirichlet()
{
    FOR_GRIDZY grid[GI(0,j,k)] = grid[GI(Nx+1,j,k)] = 0; END_FOR_TWO
    FOR_GRIDZX grid[GI(i,0,k)] = grid[GI(i,Ny+1,k)] = 0; END_FOR_TWO
    FOR_GRIDYX grid[GI(i,j,0)] = grid[GI(i,j,Nz+1)] = 0; END_FOR_TWO
    for(int k=1 ; k<=Nz ; k++ ) {
		grid[GI(0,0   ,k)] = grid[GI(Nx+1,0   ,k)] = 0;
        grid[GI(0,Ny+1,k)] = grid[GI(Nx+1,Ny+1,k)] = 0;
	}
	for(int j=1 ; j<=Ny ; j++ ) {
		grid[GI(0,j,0   )] = grid[GI(Nx+1,j,0   )] = 0;
        grid[GI(0,j,Nz+1)] = grid[GI(Nx+1,j,Nz+1)] = 0;
	}
	for(int i=1 ; i<=Nx; i++ ) { 
		grid[GI(i,0,0   )] = grid[GI(i,Ny+1,0   )] = 0;
        grid[GI(i,0,Nz+1)] = grid[GI(i,Ny+1,Nz+1)] = 0;
	}
    grid[GI(0   ,0   ,0   )] = grid[GI(Nx+1,0   ,0   )] = 0;
	grid[GI(0   ,Ny+1,0   )] = grid[GI(0   ,0   ,Nz+1)] = 0;
    grid[GI(Nx+1,Ny+1,0   )] = grid[GI(Nx+1,0   ,Nz+1)] = 0;
    grid[GI(0   ,Ny+1,Nz+1)] = grid[GI(Nx+1,Ny+1,Nz+1)] = 0;
}

inline void Grid::SetBoundaryNeumann()
{
    FOR_GRIDZY
        grid[GI(0   ,j,k)] = grid[GI(1 ,j,k)];
        grid[GI(Nx+1,j,k)] = grid[GI(Nx,j,k)];
    END_FOR_TWO
    FOR_GRIDZX
        grid[GI(i,0   ,k)] = grid[GI(i,1 ,k)];
        grid[GI(i,Ny+1,k)] = grid[GI(i,Ny,k)];
    END_FOR_TWO
    FOR_GRIDYX
        grid[GI(i,j,0   )] = grid[GI(i,j,1 )];
        grid[GI(i,j,Nz+1)] = grid[GI(i,j,Nz)];
    END_FOR_TWO
	SetBoundaryCorner();
}

inline void Grid::SetBoundaryU()
{
    FOR_GRIDZY
        grid[GI(0   ,j,k)] = -grid[GI(1 ,j,k)];
        grid[GI(Nx+1,j,k)] = -grid[GI(Nx,j,k)];
    END_FOR_TWO
    FOR_GRIDZX
        grid[GI(i,0   ,k)] = grid[GI(i,1 ,k)];
        grid[GI(i,Ny+1,k)] = grid[GI(i,Ny,k)];
    END_FOR_TWO
    FOR_GRIDYX
        grid[GI(i,j,0   )] = grid[GI(i,j,1 )];
        grid[GI(i,j,Nz+1)] = grid[GI(i,j,Nz)];
    END_FOR_TWO
	SetBoundaryCorner();
}

inline void Grid::SetBoundaryV()
{
    FOR_GRIDZY
        grid[GI(0   ,j,k)] = grid[GI(1 ,j,k)];
        grid[GI(Nx+1,j,k)] = grid[GI(Nx,j,k)];
    END_FOR_TWO
    FOR_GRIDZX
        grid[GI(i,0   ,k)] = -grid[GI(i,1 ,k)];
        grid[GI(i,Ny+1,k)] = -grid[GI(i,Ny,k)];
    END_FOR_TWO
    FOR_GRIDYX
        grid[GI(i,j,0   )] = grid[GI(i,j,1 )];
        grid[GI(i,j,Nz+1)] = grid[GI(i,j,Nz)];
    END_FOR_TWO
	SetBoundaryCorner();
}

inline void Grid::SetBoundaryW()
{
    FOR_GRIDZY
        grid[GI(0   ,j,k)] = grid[GI(1 ,j,k)];
        grid[GI(Nx+1,j,k)] = grid[GI(Nx,j,k)];
    END_FOR_TWO
    FOR_GRIDZX
        grid[GI(i,0   ,k)] = grid[GI(i,1 ,k)];
        grid[GI(i,Ny+1,k)] = grid[GI(i,Ny,k)];
    END_FOR_TWO
    FOR_GRIDYX
        grid[GI(i,j,0   )] = -grid[GI(i,j,1 )];
        grid[GI(i,j,Nz+1)] = -grid[GI(i,j,Nz)];
    END_FOR_TWO
	SetBoundaryCorner();
}

inline void Grid::SetBoundarySignedDist()
{
    static Double phi = 3. * HH;
    FOR_GRIDZY grid[GI(0,j,k)] = grid[GI(Nx+1,j,k)] = phi; END_FOR_TWO
    FOR_GRIDZX grid[GI(i,0,k)] = grid[GI(i,Ny+1,k)] = phi; END_FOR_TWO
    FOR_GRIDYX grid[GI(i,j,0)] = grid[GI(i,j,Nz+1)] = phi; END_FOR_TWO
    for(int k=1 ; k<=Nz ; k++ ) {
		grid[GI(0,0   ,k)] = grid[GI(Nx+1,0   ,k)] = phi;
        grid[GI(0,Ny+1,k)] = grid[GI(Nx+1,Ny+1,k)] = phi;
	}
	for(int j=1 ; j<=Ny ; j++ ) {
		grid[GI(0,j,0   )] = grid[GI(Nx+1,j,0   )] = phi;
        grid[GI(0,j,Nz+1)] = grid[GI(Nx+1,j,Nz+1)] = phi;
	}
	for(int i=1 ; i<=Nx; i++ ) { 
		grid[GI(i,0,0   )] = grid[GI(i,Ny+1,0   )] = phi;
        grid[GI(i,0,Nz+1)] = grid[GI(i,Ny+1,Nz+1)] = phi;
	}
    grid[GI(0   ,0   ,0   )] = grid[GI(Nx+1,0   ,0   )] = phi;
	grid[GI(0   ,Ny+1,0   )] = grid[GI(0   ,0   ,Nz+1)] = phi;
    grid[GI(Nx+1,Ny+1,0   )] = grid[GI(Nx+1,0   ,Nz+1)] = phi;
    grid[GI(0   ,Ny+1,Nz+1)] = grid[GI(Nx+1,Ny+1,Nz+1)] = phi;
}

inline void Grid::SetBoundaryCorner()
{
    for(int k=1 ; k<=Nz ; k++ ) {
		grid[GI(0   ,0   ,k)] = 0.5 * (grid[GI(1   ,0 ,k)] + grid[GI(0   ,1   ,k)]);
        grid[GI(Nx+1,0   ,k)] = 0.5 * (grid[GI(Nx  ,0 ,k)] + grid[GI(Nx+1,1   ,k)]);
        grid[GI(0   ,Ny+1,k)] = 0.5 * (grid[GI(0   ,Ny,k)] + grid[GI(1   ,Ny+1,k)]);
		grid[GI(Nx+1,Ny+1,k)] = 0.5 * (grid[GI(Nx+1,Ny,k)] + grid[GI(Nx  ,Ny+1,k)]);
	}
	for(int j=1 ; j<=Ny ; j++ ) {
		grid[GI(0   ,j,0   )] = 0.5 * (grid[GI(1   ,j,0 )] + grid[GI(0   ,j,1   )]);
        grid[GI(Nx+1,j,0   )] = 0.5 * (grid[GI(Nx  ,j,0 )] + grid[GI(Nx+1,j,1   )]);
        grid[GI(0   ,j,Nz+1)] = 0.5 * (grid[GI(0   ,j,Nz)] + grid[GI(1   ,j,Nz+1)]);
		grid[GI(Nx+1,j,Nz+1)] = 0.5 * (grid[GI(Nx+1,j,Nz)] + grid[GI(Nx  ,j,Nz+1)]);
	}
	for(int i=1 ; i<=Nx; i++ ) { 
		grid[GI(i,0   ,0   )] = 0.5 * (grid[GI(i,1   ,0 )] + grid[GI(i,0   ,1   )]);
        grid[GI(i,Ny+1,0   )] = 0.5 * (grid[GI(i,Ny  ,0 )] + grid[GI(i,Ny+1,1   )]);
        grid[GI(i,0   ,Nz+1)] = 0.5 * (grid[GI(i,0   ,Nz)] + grid[GI(i,1   ,Nz+1)]);
		grid[GI(i,Ny+1,Nz+1)] = 0.5 * (grid[GI(i,Ny+1,Nz)] + grid[GI(i,Ny  ,Nz+1)]);
	}
    grid[GI(0   ,0   ,0   )] = ONE_THIRD * (grid[GI(1   ,0   ,0   )] + 
                                            grid[GI(0   ,1   ,0   )] + 
                                            grid[GI(0   ,0   ,1   )] );
	grid[GI(Nx+1,0   ,0   )] = ONE_THIRD * (grid[GI(Nx  ,0   ,0   )] +
                                            grid[GI(Nx+1,1   ,0   )] +
                                            grid[GI(Nx+1,0   ,1   )] );
	grid[GI(0   ,Ny+1,0   )] = ONE_THIRD * (grid[GI(1   ,Ny+1,0   )] +
                                            grid[GI(0   ,Ny  ,0   )] +
                                            grid[GI(0   ,Ny+1,1   )] );
	grid[GI(0   ,0   ,Nz+1)] = ONE_THIRD * (grid[GI(1   ,0   ,Nz+1)] +
                                            grid[GI(0   ,1   ,Nz+1)] +
                                            grid[GI(0   ,0   ,Nz  )] );
    grid[GI(Nx+1,Ny+1,0   )] = ONE_THIRD * (grid[GI(Nx  ,Ny+1,0   )] +
                                            grid[GI(Nx+1,Ny  ,0   )] +
                                            grid[GI(Nx+1,Ny+1,1   )] );
    grid[GI(Nx+1,0   ,Nz+1)] = ONE_THIRD * (grid[GI(Nx  ,0   ,Nz+1)] +
                                            grid[GI(Nx+1,1   ,Nz+1)] +
                                            grid[GI(Nx+1,0   ,Nz  )] );
    grid[GI(0   ,Ny+1,Nz+1)] = ONE_THIRD * (grid[GI(1   ,Ny+1,Nz+1)] +
                                            grid[GI(0   ,Ny  ,Nz+1)] +
                                            grid[GI(0   ,Ny+1,Nz  )] );
    grid[GI(Nx+1,Ny+1,Nz+1)] = ONE_THIRD * (grid[GI(Nx  ,Ny+1,Nz+1)] +
                                            grid[GI(Nx+1,Ny  ,Nz+1)] +
                                            grid[GI(Nx+1,Ny+1,Nz  )] );
}

#endif