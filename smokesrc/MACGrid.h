// ==========================================================================
// Copyright (C) 2009 Aline Normoyle
// ==========================================================================
#ifndef MACGrid_H_
#define MACGrid_H_

#pragma warning(disable: 4244 4267 4996)
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
using namespace boost::numeric;

#include <windows.h>
#include "GL/gl.h"
#include "vec.h"
#include "GridData.h"

enum RenderMode { NONE, DENSITY, TEMPERATURE, PRESSURE };

class Camera;
class MACGrid
{
   friend MACGrid;
public:
   MACGrid();
   ~MACGrid();
   MACGrid(const MACGrid& orig);
   MACGrid& operator=(const MACGrid& orig);

   void reset();

   void draw(const Camera& c);
   void updateSources();
   void advectVelocity(double dt);
   void addExternalForces(double dt);
   void project(double dt);
   void advectTemperature(double dt);
   void advectDensity(double dt);

protected:

   // Setup
   void initialize();
    void createSolids();

   // Rendering
   struct Cube { vec3 pos; vec4 color; double dist; };
   void drawWireGrid();
   void drawSolids(const Camera& c);
   void drawSmokeCubes(const Camera& c);
   void drawSmoke(const Camera& c);
   void drawCube(const MACGrid::Cube& c);
   void drawFace(const MACGrid::Cube& c);
   void drawVelocities();
   void drawVForces();
   vec4 getRenderColor(int i, int j, int k);
   vec4 getRenderColor(const vec3& pt);
   void drawZSheets(bool backToFront);
   void drawXSheets(bool backToFront);

   // Simulation
   vec3 traceBack(const vec3& pt, double dt);
   vec3 getVelocity(const vec3& pt);
   double getVelocityX(const vec3& pt);  
   double getVelocityY(const vec3& pt); 
   double getVelocityZ(const vec3& pt);  
   double getTemperature(const vec3& pt);
   double getDensity(const vec3& pt);

   enum Direction { X, Y, Z };
   vec3 getCenter(int i, int j, int k);
   vec3 getLeftFace(int i, int j, int k);
   vec3 getRightFace(int i, int j, int k);
   vec3 getTopFace(int i, int j, int k);
   vec3 getBottomFace(int i, int j, int k);
   vec3 getFrontFace(int i, int j, int k);
   vec3 getBackFace(int i, int j, int k);
   void getCell(int index, int& i, int& j, int& k);
   int getIndex(int i, int j, int k);
   bool isNeighbor(int i0, int j0, int k0, int i1, int j1, int k1);
   bool isFace(int i, int j, int k, Direction d);
   int isSolidCell(int i, int j, int k); // Returns 1 if true, else otherwise
   int isSolidFace(int i, int j, int k, Direction d); // Returns 1 if true, else otherwise
   bool inSolid(const vec3& pt);
   bool inSolid(const vec3& pt, int& i, int& j, int& k);
   bool intersects(const vec3& pt, const vec3& dir, int i, int j, int k, double& time);
   int numSolidCells();

   // Compute pressure and divergence
   void constructB(ublas::vector<double>& b, unsigned int numCells, double dt);
   void constructA(ublas::compressed_matrix<double>& A);
   void constructPrecon(ublas::compressed_matrix<double>& A, ublas::vector<double>& precon);

   double getPressureCoeffBetweenCells(
      int i0, int j0, int k0, int i1, int j1, int k1);
   double getDivergence(int i, int j, int k);  // At center
   double checkDivergence(int i, int j, int k);  // At center
   bool checkDivergence();

   // Compute forces
   void computeBouyancy(double dt);
   void computeVorticityConfinement(double dt);
   double getBoussinesqForce(const vec3& pt);
   vec3 getVorticityN(int i, int j, int k);
   vec3 getVorticity(int i, int j, int k);
   vec3 getConfinementForce(int i, int j, int k);

   // Testing
   void testInterpolation();

   GridDataX mU;
   GridDataY mV;
   GridDataZ mW;
   CubicGridData mD;
   CubicGridData mT;
 
   GridData mSolid;

public:

   enum RenderMode { CUBES, SHEETS };
   static RenderMode theRenderMode;
   static bool theDisplayVel;
   static double theAirDensity;
   static double theAmbientTemp;
   static double theBoussinesqAlpha;
   static double theBoussinesqBeta;
   static double theVorticityConst;
   static bool theVConfEnabled;
   static bool theSourceEnabled;
};

#endif