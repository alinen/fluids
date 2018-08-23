
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
	Container2D : A class for representing and working with the levelset library
				  as a whole. 
			  
	The class holds a LevelSet, a ParticleSet, a VelocityGrid, and a timestep
	
	Functions:
	GetVelocity		- returns the velocity at the given point
	Update			- This is the actual simulator. The steps are pretty selfexplanatory
	Clear			- resets the grid to its original form
	
	MakeSphere:
	This function will initialize the values in the grid "init" to create an implicit surface 
	representing Zalesak’s disk. A function like this is necessary to create the initial grid 
	values representing the implicit surface. In the current version of this file, it is 
	assumed that this function creates a signed distance function around the implicit function.
	However, if this is not the case, the level set function "Reinitialize" has to be called
	immediately after the "Initialize" function is called. The "Reinitialize" function
	creates a signed distance function based on the interface location.
	
	Created by Emud Mokhberi: UCLA : 09/04/04
*/

#ifndef CONTAINER2D_H
#define CONTAINER2D_H

#include "main.h"
#include "LevelSet2D.h"
#include "ParticleSet2D.h"
#include "Velocity2D.h"

void MakeSphere(Grid2D &init, Double h, const Vec2D &pos, Double radius);

class Container2D
{
public:
    Container2D(int nx, int ny, Double h, Double dti) 
        : lset(nx,ny,h), pset(nx,ny,h), grid(nx,ny), init(nx,ny), dt(dti), Nx(nx), Ny(ny) 
	{ MakeSphere(init, h, (Vec2D(Nx,Ny) * Vec2D(0.5, 0.75))+Vec2D(1,1), .15 * Ny); Clear(); }
	
    inline Double& operator[] (int index) { return lset[index]; }
	inline const Double& operator[] (int index) const { return lset[index]; }
	inline Double& operator() (int i, int j) { return lset(i,j); }
	inline const Double& operator() (int i, int j) const { return lset(i,j); }

    void GetVelocity(int i, int j, Double &u, Double &v)
    {
        Vec2D uv(u,v);
        grid.GetVelocity(Vec2D(i,j),uv);
        u = uv[0];
        v = uv[1];
    }

	void Update()
	{
		//THE VELOCITY GRID WOULD BE UPDATED HERE. THIS EXAMPLE IS USING A CONTANT VELOCITY
		//GRID WHICH CREATES A RIGID BODY ROTATION OF THE IMPLICIT SURFACE
		lset.Update(grid,dt);
		pset.Update(grid,dt);
		lset.Fix(pset);
		lset.ReInitialize();
		lset.Fix(pset);
		//See Section 3.4 in the paper for activating the lines below
		//pset.Resample(lset);

		//count++;
		//if(count % 20 == 0) pset.Reseed(lset);
	}

    void Clear()
    {
		lset.Initialize(init);
		pset.Reseed(lset);
		count = 0;
    }

	LevelSet2D lset;
	ParticleSet2D pset;
	Velocity2D grid;
    int Nx, Ny;
    Double dt;
    int count;
	Grid2D init;
};

void MakeSphere(Grid2D &init, Double h, const Vec2D &pos, Double radius)
{
	FOR_ALL_LS2D
		Double val1 = ((pos - Vec2D(i,j)).Length() - radius) * h; 
        
        //Double val3 = (pos - Vec2D(i,j)).Length() - (radius - (radius*0.2));
        //val1 = max(-val3, val1);

		Double leftWall = (pos[0] - radius/6);// * h;
		Double rightWall = (pos[0] + radius/6);// * h;
		Double topWall = (pos[1] + radius/1.5);// * h;
		Double bottomWall = (pos[1] - radius);// * h;

		Double val2 = INFINITY;
		
		if(bottomWall <= j && j <= topWall && leftWall <= i && i <= rightWall)
		{
			Double top = topWall - j;
			Double bottom = j - bottomWall;
			Double left = i - leftWall;
			Double right = rightWall - i;
			if(abs(val2) > abs(top))	val2 = top;
			if(abs(val2) > abs(bottom)) val2 = bottom;
			if(abs(val2) > abs(left))	val2 = left;
			if(abs(val2) > abs(right))	val2 = right;
		}
		else if(leftWall <= i && i <= rightWall)
		{
			if(abs(val2) > abs(j-bottomWall))	val2 = j-bottomWall;
			if(abs(val2) > abs(j-topWall))		val2 = topWall-j;
		}
		else if(bottomWall <= j && j <= topWall)
		{
			if(abs(val2) > abs(i-leftWall))		val2 = i-leftWall;
			if(abs(val2) > abs(i-rightWall))	val2 = rightWall-i;
		}
		else
		{
			Double ul = sqrt((i-leftWall)*(i-leftWall)+(j-topWall)*(j-topWall));
			Double ur = sqrt((i-rightWall)*(i-rightWall)+(j-topWall)*(j-topWall));
			Double bl = sqrt((i-leftWall)*(i-leftWall)+(j-bottomWall)*(j-bottomWall));
			Double br = sqrt((i-rightWall)*(i-rightWall)+(j-bottomWall)*(j-bottomWall));
			val2 = -min(min(min(ul,ur),bl),br);
		}
        val2 *= h;

		//val = 10*cos(i/Double(size[0])*10) + 10*sin(j/Double(size[1])*10) ;
		//gridPhi(i,j) = val1;
		init(i,j) = max(val1,val2);
	END_FOR_TWO
}

#endif