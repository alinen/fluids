#ifndef GridData_H_
#define GridData_H_

#pragma warning(disable: 4244 4267 4996)
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
using namespace boost::numeric;
#include "vec.h"

// Implements arrays based on tips from Adam Bargteil
class GridData
{
public:
   GridData();
   GridData(const GridData& orig);
   virtual ~GridData();
   virtual GridData& operator=(const GridData& orig);
   virtual double& operator()(int i, int j, int k);
   virtual const double& operator()(int i, int j, int k) const;
   virtual double interpolate(const vec3& pt);
   virtual void initialize(double dfltValue = 0.0);
   ublas::vector<double>& data();
   virtual void getCell(const vec3& pt, int& i, int& j, int& k);

protected:

   virtual vec3 worldToSelf(const vec3& pt) const;
   double mDfltValue;
   vec3 mMax;
   ublas::vector<double> mData;
};

class GridDataX : public GridData
{
public:
   GridDataX();
   virtual ~GridDataX();
   virtual void initialize(double dfltValue = 0.0);
   virtual double& operator()(int i, int j, int k);
   virtual const double& operator()(int i, int j, int k) const;
   virtual vec3 worldToSelf(const vec3& pt) const;
};

class GridDataY : public GridData
{
public:
   GridDataY();
   virtual ~GridDataY();
   virtual void initialize(double dfltValue = 0.0);
   virtual double& operator()(int i, int j, int k);
   virtual const double& operator()(int i, int j, int k) const;
   virtual vec3 worldToSelf(const vec3& pt) const;
};

class GridDataZ : public GridData
{
public:
   GridDataZ();
   virtual ~GridDataZ();
   virtual void initialize(double dfltValue = 0.0);
   virtual double& operator()(int i, int j, int k);
   virtual const double& operator()(int i, int j, int k) const;
   virtual vec3 worldToSelf(const vec3& pt) const;
};

class CubicGridData : public GridData
{
public:
   CubicGridData();
   CubicGridData(const CubicGridData& orig);
   virtual ~CubicGridData();
   virtual double interpolate(const vec3& pt);

protected:

   double cubic(double q1, double q2, double q3, double q4, double t);
   double interpX(int i, int j, int k, double fracty, double fractx);
   double interpY(int i, int j, int k, double fracty);
};

class ClampedGridData : public GridData
{
public:
	ClampedGridData();
	virtual ~ClampedGridData();
	virtual double& operator()(int i, int j, int k);
   virtual const double& operator()(int i, int j, int k) const;
	//virtual double interpolate(const vec3& pt);
};


#endif