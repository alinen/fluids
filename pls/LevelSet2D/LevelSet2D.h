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
	LevelSet2D : A class for representing and working with a 2D LevelSets
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

#ifndef LEVELSET2D_H
#define LEVELSET2D_H

#include "Grid2D.h"
#include "fastmarch2d.h"
#include "main.h"
	
class LevelSet2D {
public:
	LevelSet2D(int nx,int ny, Double hi) 
        : Nx(nx), Ny(ny), h(hi), hInv(1./hi), size((nx+2)*(ny+2)), gridPhi(nx,ny), 
        gridTmp(nx,ny), gridPos(nx,ny), gridNeg(nx,ny), gridFM(nx,ny,hi) {}

    inline Double& operator[] (int index) { return gridPhi[index]; }
	inline const Double& operator[] (int index) const { return gridPhi[index]; }
	inline Double& operator() (int i, int j) { return gridPhi(i,j); }
	inline const Double& operator() (int i, int j) const { return gridPhi(i,j); }

	void Initialize(const Grid2D &init) { gridPhi = init; }
	void Update(const Velocity2D &grid, const Double &dt);
	void Fix(const ParticleSet2D &particleSet);
	void ReInitialize();
	
	inline Double LinearSample(const Vec2D &p) const;
	Double CubicSample(const Vec2D &pos) const;

private:
	inline void FixPos(const Particle2D &particle, int i, int j);
	inline void FixNeg(const Particle2D &particle, int i, int j);
	void normal(const Vec2D &pos, Vec2D &n) const;
	void gradient(const Vec2D &pos, Vec2D &g) const;
	void gradient(const Vec2D &pos, const Vec2D &u, Vec2D &g);
	void SemiLagrangianStep(int x,int y, const Velocity2D& grid, const Double &dt);

	//grid size in each dimension
	int Nx, Ny, size;
	Double h, hInv;

	Grid2D gridPhi;
	Grid2D gridTmp;
	Grid2D gridPos;
	Grid2D gridNeg;
	FastMarch2D gridFM;
};


#endif