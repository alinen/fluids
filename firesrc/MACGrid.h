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

enum RenderMode { NONE, SMOKE, TEMPERATURE, PRESSURE , LEVELSET, BLUECORE, ALL};
enum SolidMode { SNONE, FOUNTAIN, GRILL };

struct index
{
	index(int i, int j, int k)
		: i(i), j(j), k(k){}
	bool operator== (const index &rhs) const;
	int i;
	int j;
	int k;
};

class Camera;
class MACGrid
{
   friend MACGrid;
public:
   MACGrid(bool createA = true);
   ~MACGrid();
   MACGrid(const MACGrid& orig);
   MACGrid& operator=(const MACGrid& orig);

   void reset();
   void reset(SolidMode mode); 

   void draw(const Camera& c);
   void updateSources();
   void evolveFlameFront(double dt);
   void reconditionFlameFront();
   void advectVelocity(double dt);
   void addExternalForces(double dt);
   void project(double dt);
   void advectTemperature(double dt);
   void advectSmokeDensity(double dt);
   void advectColor(double dt);

   vec4 getFireRenderColor(const vec3& pt);
   double getTemperature(const vec3& pt);
   double getSmokeDensity(const vec3& pt);
   bool insideLevelSet(const vec3& pt);

protected:

   // Setup
   void initialize();
   void createSolids(SolidMode m);

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
   vec3 getAirVelocity(const vec3& pt);
   vec3 getFuelVelocity(const vec3& pt);
   
   
   vec3   getColor(const vec3& pt);
   vec3   getGhostFuelVelocity(const vec3& newpos);
   vec3   getGhostAirVelocity(const vec3& newpos);

   enum Direction { X, Y, Z };
   vec3 getCenter(int i, int j, int k);
   vec3 getLeftFace(int i, int j, int k);
   vec3 getRightFace(int i, int j, int k);
   vec3 getTopFace(int i, int j, int k);
   vec3 getBottomFace(int i, int j, int k);
   vec3 getFrontFace(int i, int j, int k);
   vec3 getBackFace(int i, int j, int k);
   void getCell(int index, int& i, int& j, int& k);
   bool isNeighbor(int i0, int j0, int k0, int i1, int j1, int k1);
   bool isFace(int i, int j, int k, Direction d);
   int isSolidCell(int i, int j, int k); // Returns 1 if true, else otherwise
   int isSolidFace(int i, int j, int k, Direction d); // Returns 1 if true, else otherwise
   bool inSolid(const vec3& pt);
   bool inSolid(const vec3& pt, int& i, int& j, int& k);
   bool intersects(const vec3& pt, const vec3& dir, int i, int j, int k, double& time);
   int numSolidCells();

   // Compute pressure and divergence for projection
   void project(double dt, ublas::vector<double>& b, double density, 
      const GridData& U, const GridData& V, const GridData& W,
      GridData& targetU, GridData& targetV, GridData& targetW);
   void constructAirB(ublas::vector<double>& b, unsigned int numCells, double dt);
   void constructFuelB(ublas::vector<double>& b, unsigned int numCells, double dt);
   void constructA();
   void setANeighbor(int ri, int rj, int rk, int ci, int cj, int ck);
   void constructCinv();
   void saveP(ublas::vector<double>& p, unsigned int numCells);

   double getGhostAirPressureConstant(int i, int j, int k);
   double getGhostFuelPressureConstant(int i, int j, int k);
   double getPressureCoeffBetweenCells(
      int i0, int j0, int k0, int i1, int j1, int k1);
   double getDivergence(const GridData& U, const GridData& V, const GridData& W,int i, int j, int k);  // At center
   double checkDivergence(const GridData& U, const GridData& V, const GridData& W,int i, int j, int k);  // At center
   bool checkDivergence();

   // Compute forces
   void computeBouyancy(double dt);
   double getBoussinesqForce(const vec3& pt, double constant);
   void computeVorticityConfinement(double dt, double constant,
     GridData& U, GridData& V, GridData& W,
     GridData& targetU, GridData& targetV, GridData& targetW);
   vec3 getVorticityN(GridData& U, GridData& V, GridData& W,int i,int j,int k);
   vec3 getVorticity(GridData& U, GridData& V, GridData& W,int i,int j,int k);
   vec3 getConfinementForce(double constant, GridData& U, GridData& V, 
      GridData& W,int i, int j, int k);

   // Compute level set
   vec3 getLevelSetNormal(int i, int j, int k);
   bool isFlameFront(int i, int j, int k);
   void InitTestFlameFrontSphere();
   void InitTestFlameFrontWall();
   void InitTestFlameFrontBottom();
   vec3 getUpwindGradient(int i, int j, int k, const vec3& w);
   bool outsideLevelSet(const vec3& pt);
   
   double computeInitDistance(int i, int j, int k);
   void updateDistance(std::multimap<double, index> &acceptedValues, std::pair<double, index> &target);
   bool outsideLevelSet(int i, int j, int k);
   bool insideLevelSet(int i, int j, int k);

   // Testing
   void testInterpolation();

   GridDataX mAirU;
   GridDataY mAirV;
   GridDataZ mAirW;

   GridDataX mFuelU;
   GridDataY mFuelV;
   GridDataZ mFuelW;

   GridData mP;
   CubicGridData mD;
   CubicGridData mT;
   CubicGridData mColorR;
   CubicGridData mColorG;
   CubicGridData mColorB;
   ClampedGridData mLevelSet;
 
   GridData mSolid;
   ublas::compressed_matrix<double> mA;
   ublas::compressed_matrix<double> mCinv;
   bool mCreateA;

public:

   static RenderMode theRenderMode;
   static double theFuelDensity;
   static double theAirDensity;
   static double theAmbientTemp;
   static double theMaxTemp;
   static double theCoolingFactor;
   static double theMaxDensity;
   static double theAirBouyancyConst;
   static double theFuelBouyancyConst;
   static double theAirVorticityConst;
   static double theFuelVorticityConst;
   static bool theVConfEnabled;
   static bool thePreconEnabled;
   static bool theDisplayVel;
   static bool theDisplayVForces;
   static bool theDisplayGrid;
   static bool theSourceEnabled;
   static bool theDisplayFlameFront;
   static double theFlameSpeed;
   static double theLevelSetSpeedRatio;
};

#endif