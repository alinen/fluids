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
	Vec2D : A class for representing and working with a 2D Vector
	
	Created by Emud Mokhberi: UCLA : 09/04/04
*/

#ifndef VECT2D_H
#define VECT2D_H
#include "main.h"

class Vec2D
{
public:
	Vec2D() { v[0] = 0; v[1] = 0; } 
	Vec2D(Float c) { v[0] = c; v[1] = c; }
	Vec2D(Float x,Float y) { v[0]=x; v[1]=y; }
	Vec2D(const Vec2D& vi) { v[0]=vi[0]; v[1]=vi[1]; }
	Vec2D(Float vec[]) { v[0]=vec[0]; v[1]=vec[1]; }
	
	void Set(Float x=0,Float y=0,Float z=0) { v[0]=x; v[1]=y; }
	
	Float& operator[] (int index) { return v[index]; }
	const Float& operator[] (int index) const { return v[index]; }
	
	operator Float*() { return &v[0]; }
	operator const Float*() { return &v[0]; }

	Float& x() { return v[0]; }
	Float& y() { return v[1]; }

	Vec2D& operator=(const Vec2D &vi) 
		{ v[0] = vi[0]; v[1] = vi[1]; return *this;}
	Vec2D& operator+=(const Vec2D &vi) 
		{ v[0] +=vi[0]; v[1] +=vi[1]; return *this;}
	Vec2D& operator-=(const Vec2D &vi) 
		{ v[0] -=vi[0]; v[1] -=vi[1]; return *this;}
	Vec2D& operator*=(const Vec2D &vi)
		{ v[0] *=vi[0]; v[1] *=vi[1]; return *this;}
	Vec2D& operator+=(Float c)
		{ v[0] += c; v[1] += c; return *this; }
	Vec2D& operator-=(Float c)
		{ v[0] -= c; v[1] -= c; return *this; }
	Vec2D& operator*=(Float c) 
		{ v[0] *= c; v[1] *= c; return *this; }
	Vec2D& operator/=(Float c) 
		{ Float cInv  = 1. / c; v[0]*=cInv; v[1]*=cInv; return *this; }

	Vec2D operator+(const Vec2D &vi) const { return Vec2D(*this)+=vi; }
	Vec2D operator-(const Vec2D &vi) const { return Vec2D(*this)-=vi; }
	Vec2D operator*(const Vec2D &vi) const { return Vec2D(*this)*=vi; }
	Vec2D operator+(Float c) const { return Vec2D(*this)+=c; }
	Vec2D operator-(Float c) const { return Vec2D(*this)-=c; }
	Vec2D operator*(Float c) const { return Vec2D(*this)*=c; }
	Vec2D operator/(Float c) const { return Vec2D(*this)/=c; }
	Vec2D operator-() const { return Vec2D(-v[0],-v[1]); }

	Float Length() const { return sqrt(SquaredLength()); }
    Float SquaredLength() const { return v[0]*v[0] + v[1]*v[1]; }
	void Normalize() { (*this) /= Length();}
	Vec2D Hat() const { return (*this)/Length(); }
	Vec2D Cross(const Vec2D& vi) const; 
	Float Dot(const Vec2D& vi) const;

private:
	Float v[3];
};

inline Vec2D operator*(Float c, const Vec2D &vi) {
	return Vec2D(c*vi[0], c*vi[1]);
}

inline Vec2D operator/(Float c, const Vec2D &vi) {
	Float cInv = 1.0/c;
	return Vec2D(cInv*vi[0], cInv*vi[1]);
}

inline Float Vec2D::Dot(const Vec2D& vi) const 
{ 
	return v[0]*vi[0] + v[1]*vi[1]; 
}

inline Float dot(const Vec2D &v1, const Vec2D &v2) {
    return v1[0] * v2[0] + v1[1] * v2[1];
}

inline ostream& operator<<(ostream& out,const Vec2D& v)
{
    return out<<v[0]<<" "<<v[1];
}

inline istream& operator>>(istream& in,Vec2D& v)
{
    return in>>v[0]>>v[1];
}

#endif