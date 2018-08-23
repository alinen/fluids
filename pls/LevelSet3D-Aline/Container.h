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
	Container: A class for representing and working with the levelset library as a whole. 
			  
	The class holds a LevelSet, a ParticleSet, a VelocityGrid, and a timestep
	
	Functions:
	GetVelocity		- returns the velocity at the given point
	Update			- This is the actual simulator. The steps are pretty selfexplanatory
	Clear			- resets the grid to its original form
	
	MakeSphere:
	This function will initialize the values in the grid "init" to create an implicit surface 
	representing ZalesakÕs sphere. A function like this is necessary to create the initial grid 
	values representing the implicit surface. In the current version of this file, it is 
	assumed that this function creates a signed distance function around the implicit function.
	However, if this is not the case, the level set function "Reinitialize" has to be called
	immediately after the "Initialize" function is called. The "Reinitialize" function
	creates a signed distance function based on the interface location.
	
	Created by Emud Mokhberi: UCLA : 09/04/04
*/

#ifndef CONTAINER_H
#define CONTAINER_H

#include "main.h"
#include "LevelSet.h"
#include "ParticleSet.h"
#include "Velocity.h"

void MakeSphere(Grid &init, Double h, const Vector &pos, Double radius);

class Container
{
public:
    Container(int nx, int ny, int nz, Double h) 
        : fm(nx,ny,nz,h), lset(nx,ny,nz,h), pset(nx,ny,nz,h), 
          grid(nx,ny), init(nx,ny,nz), dt(DT), Nx(nx), Ny(ny) 
	{ MakeSphere(init, h, (Vector(Nx,Ny,Nz) * Vector(0.5, 0.75, 0.5)) + Vector(1,1,1), .15 * Ny );
	  Clear(); }
	
    inline Double& operator[] (int index) { return lset[index]; }
	inline const Double& operator[] (int index) const { return lset[index]; }
	inline Double& operator() (int i, int j,int k) { return lset(i,j,k); }
	inline const Double& operator() (int i, int j, int k) const { return lset(i,j,k); }

    void GetVelocity(int i, int j, int k, Double &u, Double &v, Double &w)
    {
        Vector uv(u,v,w);
        grid.GetVelocity(Vector(i,j,k),uv);
        u = uv[0];
        v = uv[1];
        w = uv[2];
    }

	void Update()
	{
		//THE VELOCITY GRID WOULD BE UPDATED HERE. THIS EXAMPLE IS USING A CONTANT VELOCITY
		//GRID WHICH CREATES A RIGID BODY ROTATION OF THE IMPLICIT SURFACE
		lset.Update(grid,dt);
		pset.Update(grid,dt);
		lset.Fix(pset);
		lset.ReInitialize(fm);
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

	LevelSet lset;
    FastMarch fm;
	ParticleSet pset;
	Velocity grid;
    int Nx, Ny, Nz;
    Double dt;
    int count;
	Grid init;
};

// ASN MOD - Make a sphere in the center of the grid
void MakeSphere(Grid &init, Double h, const Vector &pos, Double radius)
{
	FOR_ALL_LS
        Double cellsize = h;
        Vector center(NX * 0.5 * cellsize, NY * 0.5 * cellsize, NZ * 0.5 * cellsize);
		Vector pos(i * cellsize, j * cellsize, k * cellsize);

        Double distance = sqrt((pos[0] - center[0])*(pos[0] - center[0]) + 
                 (pos[1] - center[1])*(pos[1] - center[1]) + 
                 (pos[2] - center[2])*(pos[2] - center[2]));

        if (distance > radius)
		    init(i,j,k) = distance - radius;
        else if (fabs(distance - radius) < 0.0001)
		    init(i,j,k) = distance - radius;
        else if (distance < radius)
		    init(i,j,k) = -(radius - distance);

    END_FOR_THREE
}

/*
void MakeSphere(Grid &init, Double h, const Vector &pos, Double radius)
{
	FOR_ALL_LS
		Double val1 = ((pos - Vector(i,j,k)).Length() - radius) * h; 
        
        //Double val3 = (pos - Vector(i,j)).Length() - (radius - (radius*0.2));
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
		//gridPhi(i,j,k) = val1;
		init(i,j,k) = max(val1,val2);
	END_FOR_THREE
}*/

#endif