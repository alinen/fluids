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


#ifndef VECTOR_H
#define VECTOR_H
#include "main.h"

class Vector
{
public:
	Vector() { v[0] = 0; v[1] = 0; v[2] = 0; } 
	Vector(Float c) { v[0] = c; v[1] = c; v[2] = c; }
	Vector(Float x,Float y, Float z) { v[0]=x; v[1]=y; v[2]=z; }
	Vector(const Vector& vi) { v[0]=vi[0]; v[1]=vi[1]; v[2]=vi[2]; }
	Vector(Float vec[]) { v[0]=vec[0]; v[1]=vec[1]; v[2]=vec[2]; }
	
	inline void Set(Float x=0,Float y=0,Float z=0) { v[0]=x; v[1]=y; v[2]=z; }
	
	inline Float& operator[] (int index) { return v[index]; }
	inline const Float& operator[] (int index) const { return v[index]; }
	
	inline operator Float*() { return &v[0]; }
	inline operator const Float*() { return &v[0]; }

	inline Float& x() { return v[0]; }
	inline Float& y() { return v[1]; }
    inline Float& z() { return v[2]; }

	inline Vector& operator=(const Vector &vi) 
		{ v[0] = vi[0]; v[1] = vi[1]; v[2] = vi[2]; return *this;}
	inline Vector& operator+=(const Vector &vi) 
		{ v[0] +=vi[0]; v[1] +=vi[1]; v[2] += vi[2]; return *this;}
	inline Vector& operator-=(const Vector &vi) 
		{ v[0] -=vi[0]; v[1] -=vi[1]; v[2] -= vi[2]; return *this;}
	inline Vector& operator*=(const Vector &vi)
		{ v[0] *=vi[0]; v[1] *=vi[1]; v[2] *= vi[2]; return *this;}
	inline Vector& operator+=(Float c)
		{ v[0] += c; v[1] += c; v[2] += c; return *this; }
	inline Vector& operator-=(Float c)
		{ v[0] -= c; v[1] -= c; v[2] -= c; return *this; }
	inline Vector& operator*=(Float c) 
		{ v[0] *= c; v[1] *= c; v[2] *= c; return *this; }
	inline Vector& operator/=(Float c) 
		{ Float cInv  = 1. / c; v[0]*=cInv; v[1]*=cInv; v[2]*=cInv; return *this; }

	inline Vector operator+(const Vector &vi) const { return Vector(*this)+=vi; }
	inline Vector operator-(const Vector &vi) const { return Vector(*this)-=vi; }
	inline Vector operator*(const Vector &vi) const { return Vector(*this)*=vi; }
	inline Vector operator+(Float c) const { return Vector(*this)+=c; }
	inline Vector operator-(Float c) const { return Vector(*this)-=c; }
	inline Vector operator*(Float c) const { return Vector(*this)*=c; }
	inline Vector operator/(Float c) const { return Vector(*this)/=c; }
	inline Vector operator-() const { return Vector(-v[0],-v[1],-v[2]); }

	inline Float Length() const { return sqrt(SquaredLength()); }
    inline Float SquaredLength() const { return v[0]*v[0] + v[1]*v[1] + v[2]*v[2]; }
	inline void Normalize() { (*this) /= Length();}
	inline Vector Hat() const { return (*this)/Length(); }
    Vector Cross(const Vector& vi) const; 
	inline Float Dot(const Vector& vi) const;

private:
	Float v[3];
};

inline Vector operator*(Float c, const Vector &vi) {
	return Vector(c*vi[0], c*vi[1], c*vi[2]);
}

inline Vector operator/(Float c, const Vector &vi) {
	Float cInv = 1.0/c;
	return Vector(cInv*vi[0], cInv*vi[1], cInv*vi[2]);
}

inline Float Vector::Dot(const Vector& vi) const 
{ 
	return v[0]*vi[0] + v[1]*vi[1] + v[2]*vi[2]; 
}

inline Vector Vector::Cross(const Vector& vi) const 
{ 
	return Vector( v[1]*vi[2] - v[2]*vi[1],
				   v[2]*vi[0] - v[0]*vi[2],
				   v[0]*vi[1] - v[1]*vi[0]); 
}

inline Float dot(const Vector &v1, const Vector &v2) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

inline Vector cross(const Vector &v1, const Vector &v2) {
    return Vector( (v1[1]*v2[2] - v1[2]*v2[1]),
                   (v1[2]*v2[0] - v1[0]*v2[2]),
                   (v1[0]*v2[1] - v1[1]*v2[0]));
}

inline ostream& operator<<(ostream& out,const Vector& v)
{
    return out<<v[0]<<" "<<v[1]<<" "<<v[2];
}

inline istream& operator>>(istream& in,Vector& v)
{
    return in>>v[0]>>v[1]>>v[2];
}

#endif