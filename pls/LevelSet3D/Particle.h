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
	Particle: A class for representing individual particles used in error correcting the Level Set
	Inputs: Position of particle within grid, Signed Distance from the interface at the center of 
		    the particle, and the inverse of the grid cell size
			
	Each particle stores its position in the grid, its sign (whether its an exterior (positive) or
	interior(negative) particle) and its radius. The radius of the particle is set so that the
	particle is touching the interface. This is adjusted based on a minimum and maximum radius.

	Functions:
	phi			- takes as input a position and returns the signed distance of the position
				  with respect to the particle, treating the particle boundary as the interface.
				  This function is used when performing error correction on a cell using this 
				  particle.
	Sign		- returns the sign of the particle, negative for interior particles and 
				  positive for exterior particles
	Radius		- returns the radius of the particle
	GetPosition - returns the position of the particle
	SetRadius   - takes as input an updated signed distance value of the particles position within
				  the grid and resets the radius of the particle based on this new value
	Update		- takes as input a velocity grid and a timestep and updates the particle
				  using standard lagrangian method with a second order accurate runga kutta
			      time integration

	Created by Emud Mokhberi: UCLA : 09/04/04
*/

#ifndef PARTICLE_H
#define PARTICLE_H

#include "main.h"
#include "Vector.h"
#include "Velocity.h"

//Radius and position are based on index of grid. They are adjusted to reflect 
//cell size h when calculating phi.
class Particle
{
public:
	Particle(const Vector &pos, const Double &phi, const Double &hInv)
	{
		sign = phi < 0 ? -1 : 1;
		radius = sign * phi * hInv;
		if(radius >RADIUS_MAX) radius =RADIUS_MAX;
		else if(radius < RADIUS_MIN) radius = RADIUS_MIN;
		position = pos;
	}

	inline Double phi(const Vector &point, const Double &h) const 
		{ return sign * (radius - (point - position).Length()) * h; }
	inline int Sign() const { return sign; }
	Double Radius() const { return radius; }
	inline void GetPosition(Vector &pos) const { pos = position; }
	inline bool SetRadius(const Double &phi, const Double &hInv)
	{
		if((phi * sign < 0.) && (abs(phi) > PARTICLE_DELETE)) return false;
		radius = abs(phi) * hInv;
		if(radius >RADIUS_MAX) radius =RADIUS_MAX;
		else if(radius < RADIUS_MIN) radius = RADIUS_MIN;
		return true;
	}
    void Update(const Velocity &grid, const Double &dt, const Double &hInv)
    {
        //RK2 update
        static Vector u, p2;
        grid.GetVelocity(position, u);
        p2 = position + u * dt * hInv;
        grid.GetVelocity(p2, u);
        p2 += u * dt * hInv;
        position = (position + p2) * 0.5;
    }

private:
	Vector position;
	int sign;
	Double radius;
};

#endif
