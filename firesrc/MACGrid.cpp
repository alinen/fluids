#include "MACGrid.h"
#include "GL/glut.h"
#include "camera.h"
#include "ConjGrad.h"
#include <math.h>
#include <map>
#include <stdio.h>

// Globals
MACGrid target(false);
extern int theDim[3];
extern double theCellSize;
double theMaxPressure = 0.0;  // Temp var for rendering

// NOTE: x -> cols, z -> rows, y -> stacks
RenderMode MACGrid::theRenderMode = ALL;
double MACGrid::theAirDensity = 0.1; // Reacted fuel, not gaseous products, from Nguyen paper
double MACGrid::theFuelDensity = 1.0; // Unreacted fuel, from Nguyen paper
double MACGrid::theAmbientTemp = 0.0;
double MACGrid::theMaxTemp = 1.0;
double MACGrid::theCoolingFactor = 0.1;
double MACGrid::theMaxDensity = 0.25;
double MACGrid::theAirBouyancyConst = 800.0;
double MACGrid::theFuelBouyancyConst = 800.0;
double MACGrid::theAirVorticityConst = 60.0;
double MACGrid::theFuelVorticityConst = 16.0;
bool MACGrid::theVConfEnabled = true;
bool MACGrid::thePreconEnabled = false;
bool MACGrid::theDisplayVel = false;
bool MACGrid::theDisplayVForces = false; 
bool MACGrid::theDisplayGrid = true;
bool MACGrid::theSourceEnabled = true;
bool MACGrid::theDisplayFlameFront = true;
double MACGrid::theFlameSpeed = 10.5; 
double MACGrid::theLevelSetSpeedRatio = 0.7;

#define FOR_EACH_CELL \
   for(int k = 0; k < theDim[MACGrid::Z]; k++)  \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_FACE \
   for(int k = 0; k < theDim[MACGrid::Z]+1; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]+1; j++) \
         for(int i = 0; i < theDim[MACGrid::X]+1; i++) 

vec3& clamp(vec3& pos)
{
   pos[0] = max(0.0, min(theDim[0]*theCellSize, pos[0]));
   pos[1] = max(0.0, min(theDim[1]*theCellSize, pos[1]));
   pos[2] = max(0.0, min(theDim[2]*theCellSize, pos[2]));
   return pos;
}

bool index::operator ==(const index &rhs) const
{
	return ((this->i == rhs.i) && (this->j == rhs.j) && (this->k == rhs.k));
}

MACGrid::MACGrid(bool createA): mCreateA(createA)
{
   initialize();
}

MACGrid::MACGrid(const MACGrid& orig)
{
   mAirU = orig.mAirU;
   mAirV = orig.mAirV;
   mAirW = orig.mAirW;
   mFuelU = orig.mFuelU;
   mFuelV = orig.mFuelV;
   mFuelW = orig.mFuelW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;
   mSolid = orig.mSolid;   
   mColorR = orig.mColorR;
   mColorG = orig.mColorG;
   mColorB = orig.mColorB;
}

MACGrid& MACGrid::operator=(const MACGrid& orig)
{
   if (&orig == this)
   {
      return *this;
   }
   mAirU = orig.mAirU;
   mAirV = orig.mAirV;
   mAirW = orig.mAirW;
   mFuelU = orig.mFuelU;
   mFuelV = orig.mFuelV;
   mFuelW = orig.mFuelW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;   
   mSolid = orig.mSolid;
   mColorR = orig.mColorR;
   mColorG = orig.mColorG;
   mColorB = orig.mColorB;

   return *this;
}

MACGrid::~MACGrid()
{
}

void MACGrid::reset()
{
   mAirU.initialize();
   mAirV.initialize();
   mAirW.initialize();

   mFuelU.initialize();
   mFuelV.initialize();
   mFuelW.initialize();

   mP.initialize();
   mD.initialize();
   mT.initialize(theAmbientTemp);
   mColorR.initialize(0.1);
   mColorG.initialize(0.1);
   mColorB.initialize(0.1);
   mLevelSet.initialize();

   //InitTestFlameFrontSphere();
   InitTestFlameFrontWall();
   //InitTestFlameFrontBottom();
}

void MACGrid::reset(SolidMode mode)
{
   reset();
   createSolids(mode);
}

void MACGrid::createSolids(SolidMode mode)
{
   mSolid.initialize();

   if (mode == FOUNTAIN)
   {
      int j = theDim[1]/2;
      int quarter1 = theDim[0]/4;
      int quarter2 = 3*theDim[0]/4;
      for (int i = 4; i < theDim[0]-4; i++)
      {
         if (i > quarter1-4 && i < quarter1+4) continue; // Create holes
         if (i > quarter2-4 && i < quarter2+4) continue; // Create holes
         for (int k = 0; k < theDim[2]; k++)
         {
            mSolid(i, j, k) = 1;
         }
      }
   }

   if (mCreateA)
   {
     constructA();
     //constructCinv();
   }
}

void MACGrid::initialize()
{
   reset();
   createSolids(SNONE);

   //testInterpolation();
   //assert(checkDivergence());
}

void MACGrid::updateSources()
{
   if (!theSourceEnabled) return;


// Initialize bottom fule source
	vec3 center = vec3(theDim[0]*theCellSize/2, 0.0, theDim[2]*theCellSize/2);
	double rad = theDim[0]/8*theCellSize; //assumes x-direction is largest dimension

	FOR_EACH_CELL
	{
      if (j == 0)
      {
		   vec3 pos = getCenter(i, j, k);
		   double d = Distance(pos, center) - rad;
		   if(d < 0)
		   {
			   //mFuelV(i, j, k) =0.5;
			   mT(i, j, k) = theMaxTemp;
			   mD(i, j, k) = 0.5;
			   mLevelSet(i, j, k) = -1.0;
		   }
      }
	}

/* Initialize fire ball
	vec3 center = vec3(theDim[0]*theCellSize/2, theDim[1]*theCellSize/2, theDim[2]*theCellSize/2);\
	double rad = theDim[0]/8*theCellSize; //assumes x-direction is largest dimension

	FOR_EACH_CELL
	{
		vec3 pos = getCenter(i, j, k);
		double d = Distance(pos, center) - rad;
		if(d < 0)
		{
			//mFuelV(i, j, k) =0.5;
			mT(i, j, k) = theMaxTemp;
			mD(i, j, k) = 0.5;
			mLevelSet(i, j, k) = -1.0;
		}
	}
*/
}

void MACGrid::evolveFlameFront(double dt)
{
	//double phi = mLevelSet.interpolate(vec3(0.0, 0.0, 0.0));
	FOR_EACH_CELL
	{
		// Implicit surface normal
		vec3 N = getLevelSetNormal(i, j, k);
		// Fuel velocity @ cell center
		vec3 uf = getFuelVelocity(getCenter(i, j, k));
		// Implicit surface velocity
		vec3 w = theLevelSetSpeedRatio*uf - theFlameSpeed*N;
		// Spatial derivative using upwind difference
		vec3 grad = getUpwindGradient(i, j, k, w);
		double curLevelSet = mLevelSet(i, j, k);
		double deltaLevlSet = dt*(Dot(w, grad));
		double newLevelSet = curLevelSet - deltaLevlSet;
		target.mLevelSet(i, j, k) = newLevelSet;
	}

	mLevelSet = target.mLevelSet;
	reconditionFlameFront();
}

vec3 MACGrid::getLevelSetNormal(int i, int j, int k)
{
	vec3 n;
	n[0] = (mLevelSet(i+1, j, k) - mLevelSet(i-1, j, k))/(2*theCellSize);
	n[1] = (mLevelSet(i, j+1, k) - mLevelSet(i, j-1, k))/(2*theCellSize);
	n[2] = (mLevelSet(i, j, k+1) - mLevelSet(i, j, k-1))/(2*theCellSize);
	n.Normalize();
	return n;
}

vec3 MACGrid::getUpwindGradient(int i, int j, int k, const vec3& w)
{
	vec3 grad;
	if (w[0] >= 0 && i > 0)
		grad[0] = (mLevelSet(i, j, k) - mLevelSet(i-1, j, k))/theCellSize;
	else
		grad[0] = (mLevelSet(i+1, j, k) - mLevelSet(i, j, k))/theCellSize;
	if (w[1] >= 0 && j > 0)
		grad[1] = (mLevelSet(i, j, k) - mLevelSet(i, j-1, k))/theCellSize;
	else
		grad[1] = (mLevelSet(i, j+1, k) - mLevelSet(i, j, k))/theCellSize;
	if (w[2] >= 0 && k > 0)
		grad[2] = (mLevelSet(i, j, k) - mLevelSet(i, j, k-1))/theCellSize;
	else
		grad[2] = (mLevelSet(i, j, k+1) - mLevelSet(i, j, k))/theCellSize;
	return grad;
}

void MACGrid::reconditionFlameFront()
{
	std::multimap<double, index> posTentative;
	std::multimap<double, index> posAccepted;
	std::multimap<double, index> negTentative;
	std::multimap<double, index> negAccepted;

	// Initialize accepted band 
	FOR_EACH_CELL
	{
		if(isFlameFront(i, j, k))
		{
			double dist = computeInitDistance(i, j, k);
			if (mLevelSet(i, j, k) > 0)
				posAccepted.insert(std::pair<double, index>(dist, index(i, j, k)));
			else
				negAccepted.insert(std::pair<double, index>(dist, index(i, j, k)));
		}
		else
		{
			if (mLevelSet(i, j, k) > 0)
				posTentative.insert(std::pair<double, index>(mLevelSet(i, j, k), index(i, j, k)));
			else
				negTentative.insert(std::pair<double, index>(-mLevelSet(i, j, k), index(i, j, k)));
		}
	}

	std::vector<index> neighbors;
	if(!posAccepted.empty())
	{
		while(!posTentative.empty())
		{
			std::pair<double, index> closestTentative = *(posTentative.begin());
			updateDistance(posAccepted, closestTentative);
			posAccepted.insert(closestTentative);
			posTentative.erase(posTentative.begin());
			index updatedCell = closestTentative.second;
			neighbors.clear();
			if((updatedCell.i+1) < theDim[0])
				neighbors.push_back(index(updatedCell.i+1, updatedCell.j, updatedCell.k));
			if((updatedCell.i-1) >= 0)
				neighbors.push_back(index(updatedCell.i-1, updatedCell.j, updatedCell.k));
			if((updatedCell.j+1) < theDim[1])
				neighbors.push_back(index(updatedCell.i, updatedCell.j+1, updatedCell.k));
			if((updatedCell.j-1) >= 0)
				neighbors.push_back(index(updatedCell.i, updatedCell.j-1, updatedCell.k));
			if((updatedCell.k+1) < theDim[2])
				neighbors.push_back(index(updatedCell.i, updatedCell.j, updatedCell.k+1));
			if((updatedCell.k-1) >= 0)
				neighbors.push_back(index(updatedCell.i, updatedCell.j, updatedCell.k-1));

			for(unsigned int n=0; n < neighbors.size(); n++)
			{
				for(std::multimap<double, index>::iterator tentativeIt = posTentative.begin(); tentativeIt != posTentative.end(); ++tentativeIt)
				{
					if(tentativeIt->second == neighbors[n])
					{
						std::pair<double, index> updatedTentative = *(tentativeIt);
						updateDistance(posAccepted, updatedTentative);
						posTentative.erase(tentativeIt);
						posTentative.insert(updatedTentative);
						break;
					}
				}
			}
		}
	}

	for(std::multimap<double, index>::iterator it = posAccepted.begin(); it != posAccepted.end(); ++it)
		mLevelSet(it->second.i, it->second.j, it->second.k) = it->first;

	if(!negAccepted.empty())
	{
		while(!negTentative.empty())
		{
			std::pair<double, index> closestTentative = *(negTentative.begin());
			updateDistance(negAccepted, closestTentative);
			negAccepted.insert(closestTentative);
			negTentative.erase(negTentative.begin());
			index updatedCell = closestTentative.second;
			neighbors.clear();
			if((updatedCell.i+1) < theDim[0])
				neighbors.push_back(index(updatedCell.i+1, updatedCell.j, updatedCell.k));
			if((updatedCell.i-1) >= 0)
				neighbors.push_back(index(updatedCell.i-1, updatedCell.j, updatedCell.k));
			if((updatedCell.j+1) < theDim[1])
				neighbors.push_back(index(updatedCell.i, updatedCell.j+1, updatedCell.k));
			if((updatedCell.j-1) >= 0)
				neighbors.push_back(index(updatedCell.i, updatedCell.j-1, updatedCell.k));
			if((updatedCell.k+1) < theDim[2])
				neighbors.push_back(index(updatedCell.i, updatedCell.j, updatedCell.k+1));
			if((updatedCell.k-1) >= 0)
				neighbors.push_back(index(updatedCell.i, updatedCell.j, updatedCell.k-1));

			for(unsigned int n=0; n < neighbors.size(); n++)
			{
				for(std::multimap<double, index>::iterator tentativeIt = negTentative.begin(); tentativeIt != negTentative.end(); ++tentativeIt)
				{
					if(tentativeIt->second == neighbors[n])
					{
						std::pair<double, index> updatedTentative = *(tentativeIt);
						updateDistance(negAccepted, updatedTentative);
						negTentative.erase(tentativeIt);
						negTentative.insert(updatedTentative);
						break;
					}
				}
			}
		}
	}

	for(std::multimap<double, index>::iterator it = negAccepted.begin(); it != negAccepted.end(); ++it)
		mLevelSet(it->second.i, it->second.j, it->second.k) = -it->first;
}

double MACGrid::computeInitDistance(int i, int j, int k)
{
	std::vector<double> dists;
	double p1 = mLevelSet(i, j, k);
	double p2;
	if(p1 > 0)
	{
		// Outside
		p2 = mLevelSet(i+1, j, k);
		if(p2 <= 0)
			dists.push_back((p1/(p1-p2))*theCellSize);
		p2 = mLevelSet(i-1, j, k);
		if(p2 <=0)
			dists.push_back((p1/(p1-p2))*theCellSize);
		p2 = mLevelSet(i, j+1, k);
		if(p2 <=0)
			dists.push_back((p1/(p1-p2))*theCellSize);
		p2 = mLevelSet(i, j-1, k);
		if(p2 <=0)
			dists.push_back((p1/(p1-p2))*theCellSize);
		p2 = mLevelSet(i, j, k+1);
		if(p2 <=0)
			dists.push_back((p1/(p1-p2))*theCellSize);
		p2 = mLevelSet(i, j, k-1);
		if(p2 <=0)
			dists.push_back((p1/(p1-p2))*theCellSize);
	}
	else
	{
		p2 = mLevelSet(i+1, j, k);
		if(p2 >= 0)
			dists.push_back((p1/(p1-p2))*theCellSize);
		p2 = mLevelSet(i-1, j, k);
		if(p2 >=0)
			dists.push_back((p1/(p1-p2))*theCellSize);
		p2 = mLevelSet(i, j+1, k);
		if(p2 >=0)
			dists.push_back((p1/(p1-p2))*theCellSize);
		p2 = mLevelSet(i, j-1, k);
		if(p2 >=0)
			dists.push_back((p1/(p1-p2))*theCellSize);
		p2 = mLevelSet(i, j, k+1);
		if(p2 >=0)
			dists.push_back((p1/(p1-p2))*theCellSize);
		p2 = mLevelSet(i, j, k-1);
		if(p2 >=0)
			dists.push_back((p1/(p1-p2))*theCellSize);
	}
	return *(std::min_element(dists.begin(), dists.end()));
}

void MACGrid::updateDistance(std::multimap<double, index> &acceptedValues, std::pair<double, index> &target)
{
	index ind = target.second;
	double targetVal = target.first;
	index rightInd(ind.i+1, ind.j, ind.k);
	index leftInd(ind.i-1, ind.j, ind.k);
	index topInd(ind.i, ind.j+1, ind.k);
	index bottomInd(ind.i, ind.j-1, ind.k);
	index frontInd(ind.i, ind.j, ind.k+1);
	index backInd(ind.i, ind.j, ind.k-1);
	double rightVal = DBL_MAX, leftVal = DBL_MAX, topVal = DBL_MAX, bottomVal = DBL_MAX, frontVal = DBL_MAX, backVal = DBL_MAX;
	for(std::multimap<double, index>::iterator it = acceptedValues.begin(); it != acceptedValues.end(); it++)
	{
		rightVal = (it->second == rightInd) ? it->first : rightVal;
		leftVal = (it->second == leftInd) ? it->first : leftVal;
		topVal = (it->second == topInd) ? it->first : topVal;
		bottomVal = (it->second == bottomInd) ? it->first : bottomVal;
		frontVal = (it->second == frontInd) ? it->first : frontVal;
		backVal = (it->second == backInd) ? it->first : backVal;
	}
	double A = 0.0, B = 0.0, C = -1 * theCellSize * theCellSize;
	if((rightVal != DBL_MAX) || (leftVal != DBL_MAX))
	{
		A += 1.0;
		double minVal = min(rightVal, leftVal);
		B -= (2*minVal);
		C += (minVal*minVal);
	}
	if((topVal != DBL_MAX) || (bottomVal != DBL_MAX))
	{
		A += 1.0;
		double minVal = min(topVal, bottomVal);
		B -= (2*minVal);
		C += (minVal*minVal);
	}
	if((frontVal != DBL_MAX) || (backVal != DBL_MAX))
	{
		A += 1.0;
		double minVal = min(frontVal, backVal);
		B -= (2*minVal);
		C += (minVal*minVal);
	}

// TODO: think about if these return statements are acceptable	
	if (A <= 0)
		return;

	// Quadratic Formula
	double discrim = B*B - 4*A*C;

	if (discrim < 0)
		return;
	double dist = (-B + sqrt(discrim))/(2*A);
	target.first = dist;
}

vec3 MACGrid::getGhostAirVelocity(const vec3& newpos)
{
   assert(insideLevelSet(newpos));

   // If newpos is inside the fuel region, this function returns the 
   // corresponding ghost velocity for an air particle
   int i,j,k; mLevelSet.getCell(newpos, i,j,k);
   vec3 normal = getLevelSetNormal(i,j,k);

   vec3 fuelvel = getFuelVelocity(newpos);
   double v_fuel = Dot(fuelvel, normal);
   double v_ghost = v_fuel + (theAirDensity/theFuelDensity - 1)*theFlameSpeed;
   vec3 ghostvel = v_ghost*normal + (fuelvel - v_fuel*normal);
   return ghostvel;
}

vec3 MACGrid::getGhostFuelVelocity(const vec3& newpos)
{
   assert(outsideLevelSet(newpos));

   // If newpos is inside the air region, this function returns the 
   // corresponding ghost velocity for a fuel particle
   int i,j,k; mLevelSet.getCell(newpos, i,j,k);
   vec3 normal = -getLevelSetNormal(i,j,k);

   vec3 airvel = getAirVelocity(newpos);
   double v_air = Dot(airvel, normal);
   double v_ghost = v_air - (theAirDensity/theFuelDensity - 1)*theFlameSpeed;
   vec3 ghostvel = v_ghost*normal + (airvel - v_air*normal);
   return ghostvel;      
}

void MACGrid::advectVelocity(double dt)
{
   // Update each face to update the air velocities
   FOR_EACH_FACE
   {
      vec3 airvel, fuelvel;
      if (isFace(i,j,k,X))
      {
         vec3 pos = getLeftFace(i,j,k);
         vec3 newpos = traceBack(pos, dt);

         if (insideLevelSet(newpos)) airvel = getGhostAirVelocity(newpos);
         else airvel = getAirVelocity(newpos);

         if (outsideLevelSet(newpos)) fuelvel = getGhostFuelVelocity(newpos);
         else fuelvel = getFuelVelocity(newpos);

         target.mAirU(i, j, k) = airvel[X];
         target.mFuelU(i, j, k) = fuelvel[X];
      }
      if (isFace(i,j,k,Y))
      {
         vec3 pos = getBottomFace(i,j,k);
         vec3 newpos = traceBack(pos, dt);

         if (insideLevelSet(newpos)) airvel = getGhostAirVelocity(newpos);
         else airvel = getAirVelocity(newpos);

         if (outsideLevelSet(newpos)) fuelvel = getGhostFuelVelocity(newpos);
         else fuelvel = getFuelVelocity(newpos);

         target.mAirV(i, j, k) = airvel[Y];
         target.mFuelV(i, j, k) = fuelvel[Y];
      }
      if (isFace(i,j,k,Z))
      {
         vec3 pos = getBackFace(i,j,k);
         vec3 newpos = traceBack(pos, dt);

         if (insideLevelSet(newpos)) airvel = getGhostAirVelocity(newpos);
         else airvel = getAirVelocity(newpos);

         if (outsideLevelSet(newpos)) fuelvel = getGhostFuelVelocity(newpos);
         else fuelvel = getFuelVelocity(newpos);

         target.mAirW(i, j, k) = airvel[Z];
         target.mFuelW(i, j, k) = fuelvel[Z];
      }
   }

   mAirU = target.mAirU;
   mAirV = target.mAirV;
   mAirW = target.mAirW;

   mFuelU = target.mFuelU;
   mFuelV = target.mFuelV;
   mFuelW = target.mFuelW;
}

void MACGrid::advectTemperature(double dt)
{
   FOR_EACH_CELL
   {
      vec3 pos = getCenter(i,j,k);
      vec3 newpos = traceBack(pos, dt);

      // 1. Get temperature from old location
      double newt = theMaxTemp; 
      if (outsideLevelSet(newpos))
      {
		   newt = getTemperature(newpos);

         // 2. Cool it down
         double coolingRatio = (newt - theAmbientTemp)/(theMaxTemp - theAmbientTemp);
         newt = newt - theCoolingFactor*pow(coolingRatio, 4);
      } 
      target.mT(i,j,k) = newt;
   }
   mT = target.mT;
}

void MACGrid::advectSmokeDensity(double dt)
{
   FOR_EACH_CELL
   {
      vec3 pos = getCenter(i,j,k);
      vec3 newpos = traceBack(pos, dt);
      double newd = theMaxDensity; 
      if (outsideLevelSet(newpos)) newd = getSmokeDensity(newpos);
      if (getTemperature(newpos) < 0.35) newd = 0.0;
      
      target.mD(i,j,k) = newd;
   }
   mD = target.mD;
}

void MACGrid::advectColor(double dt)
{
   FOR_EACH_CELL
   {
      vec3 pos = getCenter(i,j,k);
      vec3 newpos = traceBack(pos, dt);
      vec3 newc = getColor(newpos);
      target.mColorR(i,j,k) = newc[0];
      target.mColorG(i,j,k) = newc[1];
      target.mColorB(i,j,k) = newc[2];
   }
   mColorR = target.mColorR;
   mColorG = target.mColorG;
   mColorB = target.mColorB;
}

double MACGrid::getBoussinesqForce(const vec3& pos, double constant)
{
   // Use Boussinesq approximation
   // f = [0, -alpha*smokeDensity + beta*(T - T_amb), 0]  // For smoke
   // f = [0, beta*(T - T_amb), 0] // for fire
   double temperature = getTemperature(pos); 
   double yforce = constant*(temperature - theAmbientTemp);
   return yforce;
}

void MACGrid::computeBouyancy(double dt)
{
   FOR_EACH_FACE
   {
      if (isFace(i,j,k,Y))
      {
         vec3 pos = getBottomFace(i,j,k);
         if (insideLevelSet(pos))
         {
            double yforce = getBoussinesqForce(pos, theFuelBouyancyConst);
            double vel = mFuelV(i,j,k);
            vel = vel + dt*yforce;         
            target.mFuelV(i, j, k) = vel;            
         }
         else
         {
            double yforce = getBoussinesqForce(pos, theAirBouyancyConst);
            double vel = mAirV(i,j,k);
            vel = vel + dt*yforce;         
            target.mAirV(i, j, k) = vel;            
         }
      }
   }
   mFuelV = target.mFuelV;
   mAirV = target.mAirV;
}

vec3 MACGrid::getVorticityN(GridData& U, GridData& V, GridData& W, 
   int i, int j, int k)
{
   vec3 right  = getVorticity(U,V,W,i+1,j,k);
   vec3 left   = getVorticity(U,V,W,i-1,j,k);
   vec3 top    = getVorticity(U,V,W,i,j+1,k);
   vec3 bottom = getVorticity(U,V,W,i,j-1,k);
   vec3 front  = getVorticity(U,V,W,i,j,k+1);
   vec3 back   = getVorticity(U,V,W,i,j,k-1);

   double scale = 1.0/(2*theCellSize);
   double x = scale*(right.Length() - left.Length());
   double y = scale*(top.Length() - bottom.Length());
   double z = scale*(front.Length() - back.Length());

   vec3 N(x,y,z);
   return N.Normalize();
}

vec3 getVel(GridData& U, GridData& V, GridData& W, const vec3& pt)
{
   vec3 vel;
   vel[0] = U.interpolate(pt); 
   vel[1] = V.interpolate(pt); 
   vel[2] = W.interpolate(pt); 
   return vel;
}

vec3 MACGrid::getVorticity(GridData& U, GridData& V, GridData& W, 
   int i, int j, int k)
{
   vec3 right  = getVel(U,V,W,getRightFace(i+1,j,k));
   vec3 left   = getVel(U,V,W,getLeftFace(i,j,k));
   vec3 top    = getVel(U,V,W,getTopFace(i,j+1,k));
   vec3 bottom = getVel(U,V,W,getBottomFace(i,j,k));
   vec3 front  = getVel(U,V,W,getFrontFace(i,j,k+1));
   vec3 back   = getVel(U,V,W,getBackFace(i,j,k));

   double scale = 1.0/(theCellSize);
   double x = scale*(top[Z] - bottom[Z] - (front[Y] - back[Y]));
   double y = scale*(front[X] - back[X] - (right[Z] - left[Z]));
   double z = scale*(right[Y] - left[Y] - (top[X] - bottom[X]));

   return vec3(x, y, z);
}

vec3 MACGrid::getConfinementForce(double constant, 
   GridData& U, GridData& V, GridData& W, int i, int j, int k)
{
   vec3 N = getVorticityN(U,V,W, i,j,k);
   vec3 w = getVorticity(U,V,W, i,j,k);
   return constant*theCellSize*N.Cross(w);
}

void MACGrid::computeVorticityConfinement(double dt, double constant,
   GridData& U, GridData& V, GridData& W,
   GridData& targetU, GridData& targetV, GridData& targetW)
{
   GridData forcesX, forcesY, forcesZ;
   forcesX.initialize();
   forcesY.initialize();
   forcesZ.initialize();

   FOR_EACH_CELL // Calculate confinement forces
   {
      // TODO: Should I store this better?  Inverse interp?
      vec3 force = getConfinementForce(constant,U,V,W,i,j,k);
      forcesX(i,j,k) = force[0];
      forcesY(i,j,k) = force[1];
      forcesZ(i,j,k) = force[2];
   }

   //std::cout << "XForce: " << forcesX.data() << std::endl;
   //std::cout << "YForce: " << forcesY.data() << std::endl;
   //std::cout << "ZForce: " << forcesZ.data() << std::endl;

   // Update all velocities
   FOR_EACH_FACE
   {
      if (isFace(i,j,k,X))
      {
         vec3 pos = getLeftFace(i,j,k);
         double vel = U(i,j,k);
         double xforce = 0.5*(forcesX(i,j,k) - forcesX(i-1,j,k));
         vel = vel + dt*xforce;
         targetU(i, j, k) = vel;
      }

      if (isFace(i,j,k,Y))
      {
         vec3 pos = getBottomFace(i,j,k);
         double yforce = 0.5*(forcesY(i,j,k) - forcesY(i,j-1,k));
         double vel = V(i,j,k);
         vel = vel + dt*yforce;
         targetV(i, j, k) = vel;
      }

      if (isFace(i,j,k,Z))
      {
         vec3 pos = getBackFace(i,j,k);
         double zforce = 0.5*(forcesZ(i,j,k) - forcesZ(i,j,k-1));
         double vel = W(i,j,k);
         vel = vel + dt*zforce;
         targetW(i, j, k) = vel;
      }
   }
   //std::cout << "XVel: " << mU.data() << std::endl;
   //std::cout << "YVel: " << mV.data() << std::endl;
   //std::cout << "ZVel: " << mW.data() << std::endl;

   U = targetU;
   V = targetV;
   W = targetW;
}

void MACGrid::addExternalForces(double dt)
{
   computeBouyancy(dt);

   if (theVConfEnabled)
   {
      computeVorticityConfinement(dt, theAirVorticityConst, mAirU, mAirV, mAirW, 
         target.mAirU, target.mAirV, target.mAirW);

      computeVorticityConfinement(dt, theFuelVorticityConst, mFuelU, mFuelV, mFuelW, 
         target.mFuelU, target.mFuelV, target.mFuelW);
   }
}

double MACGrid::checkDivergence(const GridData& U, const GridData& V, const GridData& W, 
   int i, int j, int k)
{
   double x1 = U(i+1,j,k);
   double x0 = U(i,j,k);

   double y1 = V(i,j+1,k);
   double y0 = V(i,j,k);

   double z1 = W(i,j,k+1);
   double z0 = W(i,j,k);

   double xdiv = x1 - x0;
   double ydiv = y1 - y0;
   double zdiv = z1 - z0;
   double div = (xdiv + ydiv + zdiv)/theCellSize;
   return div;
}

double MACGrid::getDivergence(const GridData& U, const GridData& V, const GridData& W, 
   int i, int j, int k)
{
   double x1 = isSolidCell(i+1, j, k)? 0.0 : U(i+1,j,k);
   double x0 = isSolidCell(i-1, j, k)? 0.0 : U(i,j,k);

   double y1 = isSolidCell(i, j+1, k)? 0.0 : V(i,j+1,k);
   double y0 = isSolidCell(i, j-1, k)? 0.0 : V(i,j,k);

   double z1 = isSolidCell(i, j, k+1)? 0.0 : W(i,j,k+1);
   double z0 = isSolidCell(i, j, k-1)? 0.0 : W(i,j,k);

   double xdiv = x1 - x0;
   double ydiv = y1 - y0;
   double zdiv = z1 - z0;
   double div = (xdiv + ydiv + zdiv)/theCellSize;

//   printf("Cell %d %d %d = %.2f\n", i, j, k, div);
   return div;
}

double MACGrid::getGhostAirPressureConstant(int i, int j, int k)
{
   // Add discontinuity constants for each cell + neighbors that are fuel cells
   double constant = 0.0;
   double Ssqr = theFlameSpeed*theFlameSpeed;
   double C = -Ssqr*theFuelDensity*(theFuelDensity/theAirDensity - 1);

   if (insideLevelSet(getCenter(i,j,k))) constant += -6*C;
   if (insideLevelSet(getCenter(i-1,j,k))) constant += C;
   if (insideLevelSet(getCenter(i+1,j,k))) constant += C;
   if (insideLevelSet(getCenter(i,j-1,k))) constant += C;
   if (insideLevelSet(getCenter(i,j+1,k))) constant += C;
   if (insideLevelSet(getCenter(i,j,k-1))) constant += C;
   if (insideLevelSet(getCenter(i,j,k+1))) constant += C;

   return constant;
}

void MACGrid::constructAirB(ublas::vector<double>& b, unsigned int numCells, double dt)
{
   int cell = 0;
   double constant = -(theAirDensity*theCellSize*theCellSize)/dt;
   for (unsigned int index = 0; index < numCells; index++)
   {
      int i,j,k; getCell(index, i, j, k);
      if (isSolidCell(i,j,k)) continue;

      double ghostPressure = getGhostAirPressureConstant(i,j,k);
      b(cell) = constant*getDivergence(mAirU,mAirV,mAirW, i,j,k) + ghostPressure;
      cell++;
   }

   //std::cout << "Air B: " << b << std::endl;
}

double MACGrid::getGhostFuelPressureConstant(int i, int j, int k)
{
   // Add discontinuity constants for each cell + neighbors that are fuel cells
   double constant = 0.0;
   double Ssqr = theFlameSpeed*theFlameSpeed;
   double C = Ssqr*theFuelDensity*(theFuelDensity/theAirDensity - 1);

   if (outsideLevelSet(getCenter(i,j,k))) constant += -6*C;
   if (outsideLevelSet(getCenter(i-1,j,k))) constant += C;
   if (outsideLevelSet(getCenter(i+1,j,k))) constant += C;
   if (outsideLevelSet(getCenter(i,j-1,k))) constant += C;
   if (outsideLevelSet(getCenter(i,j+1,k))) constant += C;
   if (outsideLevelSet(getCenter(i,j,k-1))) constant += C;
   if (outsideLevelSet(getCenter(i,j,k+1))) constant += C;

   return constant;
}

void MACGrid::constructFuelB(ublas::vector<double>& b, unsigned int numCells, double dt)
{
   int cell = 0;
   double constant = -(theFuelDensity*theCellSize*theCellSize)/dt;
   for (unsigned int index = 0; index < numCells; index++)
   {
      int i,j,k; getCell(index, i, j, k);
      if (isSolidCell(i,j,k)) continue;

      double ghostPressure = getGhostFuelPressureConstant(i,j,k); 
      b(cell) = constant*getDivergence(mFuelU,mFuelV,mFuelW, i,j,k) + ghostPressure;
      cell++;
   }

   //std::cout << "Fuel B: " << b << std::endl;
}

void MACGrid::constructCinv()
{
   int numFluid = mA.size1();
   mCinv.resize(numFluid, numFluid, false);
   for (unsigned int i = 0; i < mA.size1(); i++)
   {
      for (unsigned int j = 0; j < mA.size2(); j++)
      {
         if (i == j) mCinv(i,j) = 1.0/mA(i,j);
         //else mCinv(i,j) = 0.0;
      }
   }
/*
  printf("Cinv---------------------------\n");
   for (unsigned int i = 0; i < mCinv.size1(); i++)
   {
      for (unsigned int j = 0; j < mCinv.size2(); j++)
      {
         int test = mCinv(i, j);
         std::cout << mCinv(i, j) << " " ; 
      }
      std::cout  << std::endl; 
   }
*/
}

int getIndex(int i, int j, int k)
{
   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];
   return col+row+stack;
}

void MACGrid::setANeighbor(int ri, int rj, int rk, int ci, int cj, int ck)
{
   if (isSolidCell(ci, cj, ck)) return;
   int r = getIndex(ri, rj, rk);
   int c = getIndex(ci, cj, ck);
   double coeff = getPressureCoeffBetweenCells(ri, rj, rk, ci, cj, ck);
   if (fabs(coeff) > 0.0001) 
   {
      mA(r,c) = coeff;         
   }
}

void MACGrid::constructA()
{
   unsigned int numCells = mSolid.data().size();
   unsigned int numFluid = numCells - numSolidCells();
   mA.resize(numFluid, numFluid, false);

   int r = 0;
   for (unsigned int row = 0; row < numCells; row++)
   {
	   if(row %100 == 0)
		std::cout << "Initializing row " << row << "/" << numCells << std::endl;
      int ri, rj, rk; getCell(row, ri, rj, rk); // Each row corresponds to a cell
      if (isSolidCell(ri, rj, rk)) continue;

      // WARNING: Test with solid cells
      setANeighbor(ri, rj, rk, ri, rj, rk);
      setANeighbor(ri, rj, rk, ri+1, rj, rk);
      setANeighbor(ri, rj, rk, ri-1, rj, rk);
      setANeighbor(ri, rj, rk, ri, rj+1, rk);
      setANeighbor(ri, rj, rk, ri, rj-1, rk);
      setANeighbor(ri, rj, rk, ri, rj, rk+1);
      setANeighbor(ri, rj, rk, ri, rj, rk-1);      
      r++;
   }
/*
   printf("A---------------------------\n");
   for (unsigned int i = 0; i < mA.size1(); i++)
   {
      for (unsigned int j = 0; j < mA.size2(); j++)
      {
         std::cout << mA(i, j) << " " ; 
      }
      std::cout  << std::endl; 
   }
*/
}

void MACGrid::saveP(ublas::vector<double>& p, unsigned int numCells)
{
   int cell = 0;
   for (unsigned int index = 0; index < numCells; index++)
   {
      int i,j,k; getCell(index, i, j, k);

	  if(!_isnan(p(cell)))  // LOOK - TravG added this to fix issue caused by cgsolve returning NaN values
	  {
		  if (!isSolidCell(i,j,k)) mP(i,j,k) = p(cell++);
		  else mP(i,j,k) = 0.0;
	  }
   }

   //std::cout << "P: " << mP.data() << std::endl;
}

void MACGrid::project(double dt)
{
   unsigned int numCells = theDim[0]*theDim[1]*theDim[2];
   unsigned int numFluid = numCells - numSolidCells();
   ublas::vector<double> b(numFluid);

   constructAirB(b, numCells, dt);
   project(dt, b, theAirDensity, mAirU, mAirV, mAirW, target.mAirU, target.mAirV, target.mAirW);

   constructFuelB(b, numCells, dt);
   project(dt, b, theFuelDensity, mFuelU, mFuelV, mFuelW, target.mFuelU, target.mFuelV, target.mFuelW);

   //std::cout << "u: " << target.mU.data() << std::endl;
   //std::cout << "v: " << target.mV.data() << std::endl;
   //std::cout << "w: " << target.mW.data() << std::endl;

   mAirU = target.mAirU;
   mAirV = target.mAirV;
   mAirW = target.mAirW;

   mFuelU = target.mFuelU;
   mFuelV = target.mFuelV;
   mFuelW = target.mFuelW;

   //checkDivergence();
/*
   printf("AFTER PROJECT: ");
   FOR_EACH_FACE
   {
      if (i == 0) printf("\n");
      if (isFace(i,j,k,Y))printf("(%.2f %.2f %.2f) ", getVelocity(getBottomFace(i,j,k))); 
   }
   printf("\n");
*/
}

void MACGrid::project(double dt, ublas::vector<double>& b, double density, 
   const GridData& U, const GridData& V, const GridData& W,
   GridData& targetU, GridData& targetV, GridData& targetW)
{
   // Solve Ax = b for pressure
   unsigned int numCells = theDim[0]*theDim[1]*theDim[2];
   unsigned int numFluid = numCells - numSolidCells();

   // 1. Contruct b - Done before function call
   // 2. Construct A - Already done since our grid is static 
   // 3. Solve for p
   ublas::vector<double> p(numFluid);
   if (thePreconEnabled) cg_psolve(mA, mCinv, b, p, 500, 0.005);
   else cg_solve(mA, b, p, 500, 0.005);
   saveP(p, numCells);
   //std::cout << "P: " << p << std::endl;

   // Subtract pressure from our velocity and save in target
   // u_new = u - dt*(1/theAirPressure)*((p_i+1-p_i)/theCellSize)
   double pressureChange; 
   double scaleConstant = dt/(theCellSize*density); 
   double C = theFlameSpeed*theFlameSpeed*theFuelDensity*(theFuelDensity/theAirDensity - 1);
   theMaxPressure = 0.0;

   FOR_EACH_FACE
   {
      if (isFace(i,j,k,X))
      {
         if (isSolidFace(i, j, k, X)) 
         {
            targetU(i, j, k) = 0.0;
         }
         else
         {
            if (mP(i,j,k) > theMaxPressure) theMaxPressure = mP(i,j,k);
            pressureChange = (mP(i,j,k) - mP(i-1,j,k));
            double vel = U(i,j,k); 
            vel = vel - scaleConstant*pressureChange;
            targetU(i, j, k) = vel;
         }           
         
      }
      if (isFace(i,j,k,Y))   
      {
         // Hard-code boundary condition for now
         if (isSolidFace(i,j,k,Y))
         {
            targetV(i, j, k) = 0.0;
         }
         else
         {
            pressureChange = (mP(i,j,k) - mP(i,j-1,k));
            double vel = V(i,j,k);
            vel = vel - scaleConstant*pressureChange;
            targetV(i, j, k) = vel;
         }
      }

      if (isFace(i,j,k,Z))
      {
         // Hard-code boundary condition for now
         if (isSolidFace(i,j,k,Z))
         {
            targetW(i, j, k) = 0.0;
         }
         else
         {
            pressureChange = (mP(i,j,k) - mP(i,j,k-1));
            double vel = W(i,j,k);
            vel = vel - scaleConstant*pressureChange;
            targetW(i, j, k) = vel;
         }
      }
   }
}

bool MACGrid::checkDivergence()
{
   FOR_EACH_CELL
   {
      double div = checkDivergence(mAirU, mAirV, mAirW, i, j, k);
      if (fabs(div) > 0.01) 
      {
         printf("Air Divergence(%d,%d,%d) = %.2f\n", i, j, k, div);
         //return false;
      }

      div = checkDivergence(mFuelU, mFuelV, mFuelW, i, j, k);
      if (fabs(div) > 0.01) 
      {
         printf("Fuel Divergence(%d,%d,%d) = %.2f\n", i, j, k, div);
         //return false;
      }
   }
   return true;
}

vec3 MACGrid::traceBack(const vec3& pt, double dt)
{
   vec3 vel = getVelocity(pt);
   vec3 pos = pt - vel*dt;

   // 1. Clamp pos to insides of our container
   pos[0] = max(0.0, min((theDim[0]-1)*theCellSize, pos[0]));
   pos[1] = max(0.0, min((theDim[1]-1)*theCellSize, pos[1]));
   pos[2] = max(0.0, min((theDim[2]-1)*theCellSize, pos[2]));

   // 2. Push point outside of our interior solids
   int i, j, k;
   if (inSolid(pt, i, j, k))
   {
      double t = 0;
      if (intersects(pt, vel, i, j, k, t))
      {
         pos = pt - vel*t;
      }
      else
      {
         cout << "WARNING: Shouldn't get here." << std::endl;
      }
   }
   return pos;
}

bool MACGrid::intersects(const vec3& wPos, const vec3& wDir, int i, int j, int k, double& time)
{
   // Transform pos/dir to local coordinates
   // Cell needs to be translated to origin and scaled by 1/theCellSize
   vec3 pos = getCenter(i, j, k);

   vec3 rayStart = wPos - pos;  // * => transform vector
   vec3 rayDir = wDir;  // ^ => transform vector; the symbol was chosen arbitrarily

   // Calculate ray/box intersection test using slabs method
   double tmin = -9999999999.0;
   double tmax = 9999999999.0;

   // Min/Max is the same in every direction for our cube
   double min = -0.5*theCellSize;
   double max =  0.5*theCellSize;

   // For X,Y,Z planes, find intersection with the minimum/maximum boundary values
   // The clever part: the maximum closest intersection will be our cube intersection
   for (int i = 0; i < 3; i++)
   {
      double e = rayStart[i];
      double f = rayDir[i];
      if (fabs(f) > 0.000000001) // Not in ith plane
      {
         double t1 = (min - e)/f;
         double t2 = (max - e)/f; 
         if (t1 > t2) swap(t1, t2);  // Always make t1 the closest for simplicity
         if (t1 > tmin) tmin = t1; 
         if (t2 < tmax) tmax = t2;
         if (tmin > tmax) return false;
         if (tmax < 0) return false;
      }
      // the ray is parallel: check if we're inside the slabs or outside
      else if (e < min || e > max) return false;
   }

   if (tmin >= 0) 
   {
      time = tmin;
      return true;
   }
   else
   {
      time = tmax;
      return true;
   }
   return false;
}

void MACGrid::getCell(int index, int& i, int& j, int& k)
{
   j = (int) index/(theDim[0]*theDim[2]);           // stack
   k = (int) (index - j*theDim[0]*theDim[2])/theDim[0];       // row
   i = index - j*theDim[0]*theDim[2] - k*theDim[0]; // col
}

vec3 MACGrid::getCenter(int i, int j, int k)
{
   double xstart = theCellSize/2.0;
   double ystart = theCellSize/2.0;
   double zstart = theCellSize/2.0;

   double x = xstart + i*theCellSize;
   double y = ystart + j*theCellSize;
   double z = zstart + k*theCellSize;
   return vec3(x, y, z);
}

vec3 MACGrid::getLeftFace(int i, int j, int k)
{
   return getCenter(i,j,k) - vec3(theCellSize*0.5, 0.0, 0.0);
}

vec3 MACGrid::getRightFace(int i, int j, int k)
{
   return getCenter(i,j,k) + vec3(theCellSize*0.5, 0.0, 0.0);
}

vec3 MACGrid::getTopFace(int i, int j, int k)
{
   return getCenter(i,j,k) + vec3(0.0, theCellSize*0.5, 0.0);
}

vec3 MACGrid::getBottomFace(int i, int j, int k)
{
   return getCenter(i,j,k) - vec3(0.0, theCellSize*0.5, 0.0);
}

vec3 MACGrid::getFrontFace(int i, int j, int k)
{
   return getCenter(i,j,k) + vec3(0.0, 0.0, theCellSize*0.5);
}

vec3 MACGrid::getBackFace(int i, int j, int k)
{
   return getCenter(i,j,k) - vec3(0.0, 0.0, theCellSize*0.5);
}

vec3 MACGrid::getVelocity(const vec3& pt)
{
   if (inSolid(pt))
   {
      //pt.Print("WARNING: Velocity given point in solid");
      return vec3Zero;
   }
   
   vec3 vel;
   if (insideLevelSet(pt))
   {
      vel = getFuelVelocity(pt);
   }
   else
   {
      vel = getAirVelocity(pt);
   }
   return vel;
}

vec3 MACGrid::getAirVelocity(const vec3& pt)
{
   if (inSolid(pt))
   {
      //pt.Print("WARNING: Velocity given point in solid");
      return vec3Zero;
   }
   
   vec3 vel;
   vel[0] = mAirU.interpolate(pt); 
   vel[1] = mAirV.interpolate(pt); 
   vel[2] = mAirW.interpolate(pt); 
   return vel;
}

vec3 MACGrid::getFuelVelocity(const vec3& pt)
{
   if (inSolid(pt))
   {
      //pt.Print("WARNING: Velocity given point in solid");
      return vec3Zero;
   }
   
   vec3 vel;
   vel[0] = mFuelU.interpolate(pt); 
   vel[1] = mFuelV.interpolate(pt); 
   vel[2] = mFuelW.interpolate(pt); 
   return vel;
}

double MACGrid::getTemperature(const vec3& pt)
{
   return mT.interpolate(pt);
}

double MACGrid::getSmokeDensity(const vec3& pt)
{
   return mD.interpolate(pt);
}

vec3 MACGrid::getColor(const vec3& pt)
{
   vec3 color;
   color[0] = mColorR.interpolate(pt);
   color[1] = mColorG.interpolate(pt);
   color[2] = mColorB.interpolate(pt);
   return color;
}

int MACGrid::numSolidCells()
{
   int numSolid = 0;
   FOR_EACH_CELL { numSolid += mSolid(i,j,k); }
   return numSolid;
}

bool MACGrid::inSolid(const vec3& pt)
{
   int i,j,k;
   mSolid.getCell(pt, i,j,k);
   return isSolidCell(i,j,k) == 1;
}

bool MACGrid::inSolid(const vec3& pt, int& i, int& j, int& k)
{
   mSolid.getCell(pt, i,j,k);
   return isSolidCell(i,j,k) == 1;
}

int MACGrid::isSolidCell(int i, int j, int k)
{
   bool containerBoundary = (i < 0 || i > theDim[0]-1) ||
                            (j < 0 || j > theDim[1]-1) ||
                            (k < 0 || k > theDim[2]-1);


   // Check interior boundaries too
   bool objectBoundary = (mSolid(i,j,k) == 1);

   return containerBoundary || objectBoundary? 1 : 0;
}

int MACGrid::isSolidFace(int i, int j, int k, MACGrid::Direction d) 
{
   if (d == X && (i==0 || i == theDim[0])) return 1;
   else if (d == Y && (j==0 || j == theDim[1])) return 1;
   else if (d == Z && (k==0 || k == theDim[2])) return 1;

   if (d == X && (mSolid(i,j,k) || mSolid(i-1,j,k))) return 1;
   if (d == Y && (mSolid(i,j,k) || mSolid(i,j-1,k))) return 1;
   if (d == Z && (mSolid(i,j,k) || mSolid(i,j,k-1))) return 1;

   return 0;
}

bool MACGrid::isNeighbor(int i0, int j0, int k0, int i1, int j1, int k1)
{
   if (abs(i0-i1) == 1 && j0 == j1 && k0 == k1) return true;
   if (abs(j0-j1) == 1 && i0 == i1 && k0 == k1) return true;
   if (abs(k0-k1) == 1 && j0 == j1 && i0 == i1) return true;
   return false;
} 

double MACGrid::getPressureCoeffBetweenCells(
   int i, int j, int k, int pi, int pj, int pk)
{
   if (i == pi && j == pj && k == pk) // self
   {
      int numSolidNeighbors = (isSolidCell(i+1,j,k) +
                               isSolidCell(i-1,j,k) +
                               isSolidCell(i,j+1,k) +
                               isSolidCell(i,j-1,k) +
                               isSolidCell(i,j,k+1) +
                               isSolidCell(i,j,k-1));                               
      // Return number of non-solid boundaries around cel ijk
      return 6.0 - numSolidNeighbors;
   }
   if (isNeighbor(i, j, k, pi, pj, pk) && !isSolidCell(pi, pj, pk)) return -1.0;
   return 0.0;
}

// *********************
//	Rendering Functions
// *********************

void MACGrid::draw(const Camera& c)
{   
   if (theDisplayGrid) drawWireGrid();
   if (theDisplayVel) drawVelocities();
   if (theDisplayVForces) drawVForces();
   drawSolids(c);
   drawSmoke(c);
   //drawSmokeCubes(c);
}

void MACGrid::drawVelocities()
{
   // draw line at each center
   glBegin(GL_LINES);
      FOR_EACH_CELL
      {
         vec3 pos = getCenter(i,j,k);
         vec3 vel = getVelocity(pos);
         if (vel.Length() > 0.0001)
         {
           vel.Normalize();
           vel *= theCellSize/2.0;
           vel += pos;

           glColor4f(1.0, 1.0, 1.0, 1.0);
           glVertex3dv(pos.n);
           glColor4f(0.0, 1.0, 0.0, 1.0);
           glVertex3dv(vel.n);
         }
      }
   glEnd();
}

void MACGrid::drawVForces()
{
   // draw line at each center
   glColor4f(1.0, 1.0, 0.0, 1.0);
   glBegin(GL_LINES);
      FOR_EACH_CELL
      {
         vec3 pos = getCenter(i,j,k);
         vec3 force = vec3(0.0,0.0,0.0); // ASN TODO: getConfinementForce(i,j,k);
         if (force.Length() > 0.0001)
         {
           force.Normalize();
           force *= theCellSize/2.0;
           force += pos;
           glVertex3dv(pos.n);
           glVertex3dv(force.n);
         }
      }
   glEnd();
}

static std::map<double, vec4> theColors;
void initBlackBodyColors()
{
   double scale = 1/255.0;
   theColors[1000] = scale*vec4(248,38,0,128);
   theColors[2000] = scale*vec4(255,160,0,128);
   theColors[3000] = scale*vec4(255,171,31,128);
   theColors[4000] = scale*vec4(255,196,98,128);
   theColors[5000] = scale*vec4(255,221,165,128);
   //theColors[6000] = scale*vec4(255,246,232,128);
   //theColors[7000] = scale*vec4(245,247,255,128);
}

vec4 MACGrid::getFireRenderColor(const vec3& pt)
{   
   double temp = getTemperature(pt);
   double value = getSmokeDensity(pt);
   if (temp < 0.5)
   {
      return vec4(0.5, 0.25, 0.1, value); // OR 0.0, 0.0, 0.0 for black
   }

   temp = temp*2.0 - 1;
   // Scale temperature between 1000K and 5000K
   double tempK = 1000*(1-temp) + 5000*(temp);

   if (theColors.size() == 0) initBlackBodyColors();
   std::map<double, vec4>::const_iterator it;
   for (it = theColors.begin(); it != theColors.end(); ++it)
   {
      if (tempK < it->first)
      {
         vec4 color;
         double etemp = it->first;
         vec4 ecolor = it->second;

         --it;
         vec4 bcolor = it->second;
         double btemp = it->first;
         double fract = (tempK - btemp)/(etemp - btemp);
         color = bcolor*(1-fract) + ecolor*fract;  

         color[3] = value+0.3;
         return color;
      }
   }

   return theColors[5000];
}

vec4 MACGrid::getRenderColor(int i, int j, int k)
{
   vec3 pt = getCenter(i,j,k);
   return getRenderColor(pt);
}

vec4 MACGrid::getRenderColor(const vec3& pt)
{
   switch (theRenderMode)
   {
   case ALL:
        {
        return getFireRenderColor(pt);
        }
   case SMOKE:
        {  
        vec3 color = vec3(0.5, 0.5, 0.5);
        double value = getSmokeDensity(pt);
        return vec4(color[0], color[1], color[2], value);
        }
   case PRESSURE: 
        {
        double value = mP.interpolate(pt)/theMaxPressure; 
        if (value < 0) return vec4(0.0, -value, -value, 0.5);
        else return vec4(value, 0.0, value, 0.5);
        }
   case TEMPERATURE: 
        {
        double value = getTemperature(pt); 
        return vec4(value, 1.0-value, 0.0, 0.5);
		
        }
   case LEVELSET:
	   {
		double value = mLevelSet.interpolate(pt);
		if(value < 0)
		{
			value = 1.0 + max(value, -1.0);
			return vec4(0.0, 1.0-value, value, 0.5);
		}
		else
		{
			value = 1.0 - min(value, 1.0);
			return vec4(1.0-value, 0.0f, value, 0.5);
		}
	   }
   case BLUECORE:
	   {
		   double value = mLevelSet.interpolate(pt);
		   vec4 color;
		   color = (fabs(value) <= theCellSize) ? vec4(0.0, 0.0, 1.0, 0.5) : vec4(0.0, 0.0, 0.0, 0.0);
		   return color;
	   }
   }
   return vec4(0,0,0,0);
}

void MACGrid::drawZSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double startk = back - stepsize;
   double endk = 0;
   double stepk = -theCellSize;

   if (!backToFront)
   {
      startk = 0;
      endk = back;   
      stepk = theCellSize;
   }

   for (double k = startk; backToFront? k > endk : k < endk; k += stepk)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double i = 0.0; i <= right; i += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double i = right; i >= 0.0; i -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawXSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double starti = right - stepsize;
   double endi = 0;
   double stepi = -theCellSize;

   if (!backToFront)
   {
      starti = 0;
      endi = right;   
      stepi = theCellSize;
   }

   for (double i = starti; backToFront? i > endi : i < endi; i += stepi)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double k = 0.0; k <= back; k += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double k = back; k >= 0.0; k -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}


void MACGrid::drawSmoke(const Camera& c)
{
   if (theRenderMode == NONE) return;

   vec3 eyeDir = c.getBackward();
   double zresult = fabs(Dot(eyeDir, vec3(1,0,0)));
   double xresult = fabs(Dot(eyeDir, vec3(0,0,1)));
   double yresult = fabs(Dot(eyeDir, vec3(0,1,0)));

   if (zresult < xresult)
   {      
      drawZSheets(c.getPosition()[2] < 0);
   }
   else 
   {
      drawXSheets(c.getPosition()[0] < 0);
   }
}

void MACGrid::drawSmokeCubes(const Camera& c)
{
   if (theRenderMode == NONE) return;

   std::multimap<double, MACGrid::Cube, std::greater<double>> cubes;
   FOR_EACH_CELL
   {
      MACGrid::Cube cube;
      cube.color = getRenderColor(i,j,k);
      cube.pos = getCenter(i,j,k);
      cube.dist = DistanceSqr(cube.pos, c.getPosition());
      cubes.insert(make_pair(cube.dist, cube));
   } 

   // Draw cubes from back to front
   std::multimap<double, MACGrid::Cube, std::greater<double>>::const_iterator it;
   for (it = cubes.begin(); it != cubes.end(); ++it)
   {
      drawFace(it->second);
   }
}

void MACGrid::drawSolids(const Camera& c)
{
   FOR_EACH_CELL
   {
      MACGrid::Cube cube;
      cube.color.set(0.5, 0.5, 0.5, 1.0);
      if (mSolid(i,j,k))
      {
         cube.pos = getCenter(i,j,k);
         drawCube(cube);
      }
   }
}

void MACGrid::drawWireGrid()
{
   // Display grid in light grey, draw top & bottom

   double xstart = 0.0;
   double ystart = 0.0;
   double zstart = 0.0;
   double xend = theDim[0]*theCellSize;
   double yend = theDim[1]*theCellSize;
   double zend = theDim[2]*theCellSize;

   glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
      glDisable(GL_LIGHTING);
      glColor3f(0.25, 0.25, 0.25);

      glBegin(GL_LINES);
      for (int i = 0; i <= theDim[0]; i++)
      {
         double x = xstart + i*theCellSize;
         glVertex3d(x, ystart, zstart);
         glVertex3d(x, ystart, zend);

         glVertex3d(x, yend, zstart);
         glVertex3d(x, yend, zend);
      }

      for (int i = 0; i <= theDim[2]; i++)
      {
         double z = zstart + i*theCellSize;
         glVertex3d(xstart, ystart, z);
         glVertex3d(xend, ystart, z);

         glVertex3d(xstart, yend, z);
         glVertex3d(xend, yend, z);
      }

      glVertex3d(xstart, ystart, zstart);
      glVertex3d(xstart, yend, zstart);

      glVertex3d(xend, ystart, zstart);
      glVertex3d(xend, yend, zstart);

      glVertex3d(xstart, ystart, zend);
      glVertex3d(xstart, yend, zend);

      glVertex3d(xend, ystart, zend);
      glVertex3d(xend, yend, zend);
      glEnd();
   glPopAttrib();

   glEnd();
}

#define LEN 0.5
void MACGrid::drawFace(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);
      glEnd();
   glPopMatrix();
}

bool MACGrid::isFace(int i, int j, int k, MACGrid::Direction d)
{
   switch (d)
   {
   case X: return (i >= 0 && i < theDim[X]+1 &&
                   j >= 0 && j < theDim[Y] &&
                   k >= 0 && k < theDim[Z]);
   case Y: return (i >= 0 && i < theDim[X] &&
                   j >= 0 && j < theDim[Y]+1 &&
                   k >= 0 && k < theDim[Z]);
   case Z: return (i >= 0 && i < theDim[X] &&
                   j >= 0 && j < theDim[Y] &&
                   k >= 0 && k < theDim[Z]+1);
   }
   printf("Error: bad direction passed to isFace\n");
   return false;
}

void MACGrid::testInterpolation()
{
   printf("WARNING: testing interpolations\n");
   FOR_EACH_FACE
   {
      if (isFace(i,j,k,X)) mFuelU(i, j, k) = i*theCellSize;
      if (isFace(i,j,k,Y)) mFuelV(i, j, k) = j*theCellSize;
      if (isFace(i,j,k,Z)) mFuelW(i, j, k) = k*theCellSize; 
   }

   FOR_EACH_CELL
   {
      mD(i, j, k) = i + j + k;
   }
   std::cout << "D: " << mD.data() << std::endl;

   vec3 test = getFuelVelocity(vec3(0.25,0.25,0.2));  
   vec3 test1 = getFuelVelocity(vec3(1.5,2.2,0.1));  
   vec3 test2 = getFuelVelocity(vec3(0.0,0.1,0.1));  

   double test3 = getSmokeDensity(vec3(0.0, 0.0, 0.0)); // should be 0
   double test4 = getSmokeDensity(vec3(0.5, 0.5, 0.0));  // should be 0
   double test5 = getSmokeDensity(vec3(1.0, 1.0, 0.0));  // should be 1
}


void MACGrid::drawCube(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0, -1.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN, -LEN);         

         glNormal3d( 0.0,  0.0, -0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN,  LEN, -LEN);
         glVertex3d( LEN,  LEN, -LEN);
         glVertex3d( LEN, -LEN, -LEN);

         glNormal3d(-1.0,  0.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d(-LEN,  LEN,  LEN);
         glVertex3d(-LEN,  LEN, -LEN);

         glNormal3d( 0.0, 1.0,  0.0);
         glVertex3d(-LEN, LEN, -LEN);
         glVertex3d(-LEN, LEN,  LEN);
         glVertex3d( LEN, LEN,  LEN);
         glVertex3d( LEN, LEN, -LEN);

         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);

         glNormal3d(1.0,  0.0,  0.0);
         glVertex3d(LEN, -LEN, -LEN);
         glVertex3d(LEN, -LEN,  LEN);
         glVertex3d(LEN,  LEN,  LEN);
         glVertex3d(LEN,  LEN, -LEN);
      glEnd();
   glPopMatrix();
}

bool MACGrid::outsideLevelSet(const vec3& pt)
{
   return mLevelSet.interpolate(pt) > 0;
}

bool MACGrid::insideLevelSet(const vec3& pt)
{
   return mLevelSet.interpolate(pt) <= 0;
}

bool MACGrid::outsideLevelSet(int i, int j, int k)
{
   return mLevelSet(i, j, k) > 0;
}

bool MACGrid::insideLevelSet(int i, int j, int k)
{
   return mLevelSet(i, j, k) <= 0;
}

bool MACGrid::isFlameFront(int i, int j, int k)
{
	if(mLevelSet(i, j, k) > 0)
	{
		// Outside
		bool rightInside = mLevelSet(i+1, j, k) <= 0;
		bool leftInside = mLevelSet(i-1, j, k) <= 0;
		bool topInside = mLevelSet(i, j+1, k) <= 0;
		bool bottomInside = mLevelSet(i, j-1, k) <= 0;
		bool frontInside = mLevelSet(i, j, k+1) <= 0;
		bool backInside = mLevelSet(i, j, k-1) <= 0;
		return rightInside || leftInside || topInside || bottomInside || frontInside || backInside;
	}
	else if(mLevelSet(i, j, k) < 0)
	{
		// Inside
		bool rightOutside = mLevelSet(i+1, j, k) >= 0;
		bool leftOutside = mLevelSet(i-1, j, k) >= 0;
		bool topOutside = mLevelSet(i, j+1, k) >= 0;
		bool bottomOutside = mLevelSet(i, j-1, k) >= 0;
		bool frontOutside = mLevelSet(i, j, k+1) >= 0;
		bool backOutside = mLevelSet(i, j, k-1) >= 0;
		return rightOutside || leftOutside || topOutside || bottomOutside || frontOutside || backOutside;
	}
	else
		return true;
}

void MACGrid::InitTestFlameFrontBottom()
{
	vec3 center = vec3(theDim[0]*theCellSize/2, 0.0, theDim[2]*theCellSize/2);
	double rad = theDim[0]/6*theCellSize; //assumes x-direction is largest dimension

	FOR_EACH_CELL
	{
		vec3 pos = getCenter(i, j, k);
		double d = Distance(pos, center) - rad;
		mLevelSet(i, j, k) = d;
	}
}

void MACGrid::InitTestFlameFrontSphere()
{
	vec3 center = vec3(theDim[0]*theCellSize/2, theDim[1]*theCellSize/2, theDim[2]*theCellSize/2);
	double rad = theDim[0]/6*theCellSize; //assumes x-direction is largest dimension

	FOR_EACH_CELL
	{
		vec3 pos = getCenter(i, j, k);
		double d = Distance(pos, center) - rad;
		mLevelSet(i, j, k) = d;
	}
}

void MACGrid::InitTestFlameFrontWall()
{
	double height = theDim[1]*theCellSize/2;
	
	FOR_EACH_CELL
	{
		//if((i != 0) && (i != theDim[0]-1))
		//{
			vec3 pos = getCenter(i, j, k);
			double d = pos[1]-height;
			mLevelSet(i, j, k) = d;
		//}
	}
   std::cout << "Levelset: " << mLevelSet.data() << std::endl;
}