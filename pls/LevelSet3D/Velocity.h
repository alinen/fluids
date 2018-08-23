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
	Velocity2: A class that generates a rigib body rotation within the grid
			   It is used for testing the functionality of the library

	Created by Emud Mokhberi: UCLA : 09/04/04
*/

#ifndef VELOCITY_H
#define VELOCITY_H

#include "main.h"
#include "Vector.h"

class Velocity
{
public:
	Velocity(int x,int y)
	{
		xs = Double(x+2) * 0.5;
		ys = Double(y+2) * 0.5;
		c = M_PI / 314;
	}
	inline void GetVelocity(const Vector &pos, Vector &u) const
	{
		//vortex around 0,0,1
		u = Vector( c * (ys - pos[1]), c * (pos[0] - xs), 0.);
	}
	Double xs,ys;
	Double c;
};

#endif