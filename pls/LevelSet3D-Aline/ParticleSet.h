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
	ParticleSet: A class for representing the set of particles used in error correcting the Level Set
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

#ifndef PARTICLESET_H
#define PARTICLESET_H

#include "main.h"
#include "Vector.h"

#include "LevelSet.h"
#include "Particle.h"
#include "Velocity.h"

class ParticleSet
{
private:
	int Nx, Ny, Nz;
	list<Particle*> particles;
    Double h, hInv;
public:
	ParticleSet(int nx,int ny, int nz, Double hi) : Nx(nx), Ny(ny), Nz(nz), h(hi), hInv(1./hi) {}
	~ParticleSet() 
	{   for(Iterator it = particles.begin(); it != particles.end(); it++) delete *it; 
        particles.clear(); }

	typedef list<Particle*>::iterator Iterator;
	typedef list<Particle*>::const_iterator cIterator;
	cIterator begin() const { return particles.begin(); }
	cIterator end()   const { return particles.end(); }

	void Update(const Velocity& grid, const Double &dt)
	{
        static Particle* p;
        static Vector pos;
		for(Iterator it = particles.begin(); it != particles.end();)
		{
			p = *it;
            p->Update(grid, dt, hInv);
			p->GetPosition(pos);

			if( pos[0] < 0 || pos[0] > Nx+1 || 
                pos[1] < 0 || pos[1] > Ny+1 ||
                pos[2] < 0 || pos[2] > Nz+1) { it = particles.erase(it); delete p; }
			else it++;
		}
	}
	void Resample(const LevelSet& levelSet) // updates particle radii;
	{
		static Vector pos;
        Particle *p;
		for(Iterator it = particles.begin(); it != particles.end();) 
		{
			p = *it;
			p->GetPosition(pos);
			if(p->SetRadius(levelSet.SAMPLEPHI(pos), hInv)) it++;
			else { it = particles.erase(it); delete p; }
		}
	}
	void Reseed(const LevelSet& levelSet)	 // deletes particles and creates new ones
	{
        static Double phi, ppn; 
        static bool reseed, reseed2;
		for(Iterator it = particles.begin(); it != particles.end(); it++) delete *it;
		particles.clear();

		FOR_LS
			reseed = false;
            reseed2 = false; 
            for(int dx=0; dx < 2; dx++) {
                for(int dy=0; dy < 2; dy++) {
                    for(int dz=0; dz < 2; dz++) {
                    phi = abs(levelSet(i+dx,j+dy,k+dz));
					if(phi < RESEED_THRESHOLD) reseed = true;
                    if(phi < h) reseed2 = true;
                    }
                }
            }
            if(reseed2) ppn = PARTICLES_PER_INTERFACE_NODE;
            else        ppn = PARTICLES_PER_NODE;
			if(reseed) {
				for(int x=0; x < ppn; x++) {
					Vector pos(Double(i) + RandomFloat(), 
							   Double(j) + RandomFloat(),
                               Double(k) + RandomFloat());
					particles.push_back(new Particle(pos, levelSet.SAMPLEPHI(pos), hInv));
				}
			}
               
        END_FOR_THREE
	}
};

#endif