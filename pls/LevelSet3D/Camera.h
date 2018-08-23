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

#ifndef CAMERA_H
#define CAMERA_H

#include "main.h"
#include "Vector.h"

class Camera 
{
public:
	Camera(Float x,Float y,Float z);
	void MoveForward(Float dist);
	void MoveBackward(Float dist);
	void MoveLeft(Float dist);
	void MoveRight(Float dist);
	void MoveUp(Float dist);
	void MoveDown(Float dist);

	void RotateY(Float degrees);
	void RotateX(Float degrees);

    Vector GetDir();
    Vector GetPos();

	void Draw();
private:
	void CalcDir();
	Float theta,phi;
	Vector pos;
	Vector dir,up,right;
	bool dirty;

};


#endif