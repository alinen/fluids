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
	ParticleSet2D : A class for representing the set of particles used in error correcting the Level Set
	Inputs: Size of grid and cell
			
	This class only stores a list with pointers to all the particles within the grid. Note that the
	particles are not stored per cell, but for the entire grid. 

	Functions:
	Update		- takes as input a velocity grid and a timestep. Calls the update function for
				  each partricle and removes it if the particle has exited the grid.
	Resample	- Updates the radius for each particle. Only use this function is necessary
	Reseed		- Deletes all particles and creates new ones. Only use this function when
			      absolutely necessary
				  
	Created by Emud Mokhberi: UCLA : 09/04/04
*/

#ifndef PARTICLESET2D_H
#define PARTICLESET2D_H

#include "main.h"
#include "Vec2D.h"

#include "LevelSet2D.h"
#include "Particle2D.h"
#include "Velocity2D.h"

class ParticleSet2D
{
private:
	int Nx, Ny;
	list<Particle2D*> particles;
    Double hInv;
public:
	ParticleSet2D(int nx,int ny, Double hi) : Nx(nx), Ny(ny), hInv(1./hi) {}
	~ParticleSet2D() 
	{   for(Iterator it = particles.begin(); it != particles.end(); it++) delete *it; 
        particles.clear(); }

	typedef list<Particle2D*>::iterator Iterator;
	typedef list<Particle2D*>::const_iterator cIterator;
	cIterator begin() const { return particles.begin(); }
	cIterator end()   const { return particles.end(); }

	void Update(const Velocity2D& grid, const Double &dt)
	{
        static Particle2D* p;
        static Vec2D pos;
		for(Iterator it = particles.begin(); it != particles.end();)
		{
			p = *it;
            p->Update(grid, dt, hInv);
			p->GetPosition(pos);

			if( pos[0] < 0 || pos[0] > Nx+1 || 
				pos[1] < 0 || pos[1] > Ny+1 ) { it = particles.erase(it); delete p; }
			else it++;
		}
	}
	void Resample(const LevelSet2D& levelSet) // updates particle radii;
	{
		static Vec2D pos;
        Particle2D *p;
		for(Iterator it = particles.begin(); it != particles.end();) 
		{
			p = *it;
			p->GetPosition(pos);
			if(p->SetRadius(levelSet.LinearSample(pos), hInv)) it++;
			else { it = particles.erase(it); delete p; }
		}
	}
	void Reseed(const LevelSet2D& levelSet)	 // deletes particles and creates new ones
	{
		for(Iterator it = particles.begin(); it != particles.end(); it++) delete *it;
		particles.clear();
		static bool reseed;

		FOR_LS2D
			reseed = false;
			for(int dx=0; dx < 2; dx++)
				for(int dy=0; dy < 2; dy++)
					if(abs(levelSet(i+dx,j+dy)) < RESEED_THRESHOLD)
						reseed = true;
			if(reseed) {
				for(int x=0; x < PARTICLES_PER_NODE; x++) {
					Vec2D pos(Double(i) + RandomFloat(), 
							  Double(j) + RandomFloat());
					particles.push_back(new Particle2D(pos, levelSet.LinearSample(pos), hInv));
				}
			}
        END_FOR_TWO
	}
};

#endif