
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
	LevelSet: A class for representing and working with a 3D LevelSets
	Inputs: Grid size and cell size. For the sake of simplicity, cells are of uniform size
	
	Four different Grids are stored by the class at one time. 
	gridPhi  - contians the current levelset at any one time
	gridTemp - used for updating gridPhi
	gridPos  - used for error correction with escaped positive particles
	gridNeg  - used for error correction with escaped negative particles
	
	Finally a fastMarching grid is used for reinitializing the signed distance function

	Public Functions:
	Initialize		- This function has to be called before the simulation begins and is
					  used to initialize the level set grid values
	Update			- This function takes as input a velocity grid and timestep and updates
					  the levelset using a fast first order accurate semi-lagrangian update
	Fix				- Takes a particleSet as input and performs error correction on levelSet
	ReInitialize	- Reinitializes the grid to a signed distance grid using the fast
					  first order accurate fast marching method
	LinearSample	- Takes as input a Float position within the grid and uses the four 
					  surrounding cells for linearly interpolating the value of the 
					  LevelSet at that point
	CubicSample		- Same as LinearSample but uses Cubic interpolation. Although this is
					  more accurate, it is considerable more expensive
	eval			- A function used by Marching Cubes for visualization
			
	Private Functions:
	FixPos			- Takes as input a cell the contains error and the positive particle 
					  used to correct the value of the cell
	FixNeg			- Same as FixPos but used for negative escaped particles
	normal			- Calculates the normal at the given point. This is basically the 
					  normalized gradient
	gradient		- Calculated the gradient at the given point. if a velocity is given
					  it is used in the calculation of the gradient
	SemiLagrangiaStep   - Performs the first order accurate semi lagrangian step

	Created by Emud Mokhberi: UCLA : 09/04/04
*/

#ifndef LEVELSET_H
#define LEVELSET_H

#include "Grid.h"
#include "FastMarch.h"
#include "main.h"
#include "impSurface.h"
#include "vector.h"
	
class LevelSet: public ImpSurface {
public:
	LevelSet(int nx,int ny, int nz, Double hi) 
        : Nx(nx), Ny(ny), Nz(nz), h(hi), hInv(1./hi), size((nx+2)*(ny+2)*(nz+2)), 
        gridPhi(nx,ny,nz), gridTmp(nx,ny,nz), gridPos(nx,ny,nz), gridNeg(nx,ny,nz) {}

    inline Double& operator[] (int index) { return gridPhi[index]; }
	inline const Double& operator[] (int index) const { return gridPhi[index]; }
	inline Double& operator() (int i, int j, int k) { return gridPhi(i,j,k); }
	inline const Double& operator() (int i, int j, int k) const { return gridPhi(i,j,k); }

	void Initialize(const Grid &init) { gridPhi = init; }
	void Update(const Velocity &grid, const Double &dt);
	void Fix(const ParticleSet &particleSet);
	void ReInitialize(FastMarch &gridFM);
	void gradient(const Vector &pos, Vector &g) const;
	void gradient(const Vector &pos, const Vector &u, Vector &g);
	
	inline Double LinearSample(const Vector &pos) const;
	Double CubicSample(const Vector &pos) const;

    virtual Double	eval	(const Point3d& location)
	{
		Vector pos((location[0] + 1) * (Nx-1) * 0.5 + 1, 
				   (location[1] + 1) * (Ny-1) * 0.5 + 1,
				   (location[2] + 1) * (Nz-1) * 0.5 + 1);
		return LinearSample(pos);
	}

private:
	inline void FixPos(const Particle &particle, int i, int j, int k);
	inline void FixNeg(const Particle &particle, int i, int j, int k);
	void normal(const Vector &pos, Vector &n) const;
	void SemiLagrangianStep(int x,int y, int z, const Velocity& grid, const Double &dt);

	//grid size in each dimension
	int Nx, Ny, Nz, size;
	Double h, hInv;

	Grid gridPhi;
	Grid gridTmp;
	Grid gridPos;
	Grid gridNeg;
};


#endif