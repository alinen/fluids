#include "MACGrid.h"
#include "GL/glut.h"
#include "camera.h"
#include "ConjGrad.h"
#include <math.h>
#include <map>
#include <stdio.h>

ublas::compressed_matrix<double> A;
ublas::vector<double> precon;

// Globals
MACGrid target;
extern int theDim[3];
extern double theCellSize;

// NOTE: x -> cols, z -> rows, y -> stacks
MACGrid::RenderMode MACGrid::theRenderMode = SHEETS;
bool MACGrid::theDisplayVel = true;
double MACGrid::theAirDensity = 1.0; //(1.3 kg/m^2) in reality,
double MACGrid::theAmbientTemp = 0.0;
double MACGrid::theBoussinesqAlpha = 500.0;
double MACGrid::theBoussinesqBeta = 2500.0;
double MACGrid::theVorticityConst = 100.0;
bool MACGrid::theVConfEnabled = true;
bool MACGrid::theSourceEnabled = true;

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

MACGrid::MACGrid()
{
   initialize();
}

MACGrid::MACGrid(const MACGrid& orig)
{
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mD = orig.mD;
   mT = orig.mT;
   mSolid = orig.mSolid;   
}

MACGrid& MACGrid::operator=(const MACGrid& orig)
{
   if (&orig == this)
   {
      return *this;
   }
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mD = orig.mD;
   mT = orig.mT;   
   mSolid = orig.mSolid;

   return *this;
}

MACGrid::~MACGrid()
{
}

void MACGrid::reset()
{
   mU.initialize();
   mV.initialize();
   mW.initialize();
   mD.initialize();
   mT.initialize(theAmbientTemp);

   // Set default vel to make things more interesting and avoid degenerate fluid case
//#ifndef _DEBUG
   int sourcei = (int) theDim[0]/2;
   int sourcek = (int) theDim[2]/2;
   mU(sourcei, 0, sourcek) = cos(0.3);
   mV(sourcei, 0, sourcek) = sin(0.3);
//#endif
}

void MACGrid::createSolids()
{
    mSolid.initialize();

    int j = theDim[1]/2;
    int quarter1 = theDim[0]/4;
    int quarter2 = 3*theDim[0]/4;
    int size = 4; //theDim[0]/10;
    for (int i = size; i < theDim[0]-size; i++)
    {
        if (i > quarter1-size && i < quarter1+size) continue; // Create holes
        if (i > quarter2-size && i < quarter2+size) continue; // Create holes
        for (int k = 0; k < theDim[2]; k++)
        {
            mSolid(i, j, k) = 1;
        }
    }
}

void MACGrid::initialize()
{
   reset();
   createSolids();
   constructA(A);
   constructPrecon(A, precon);

   //testInterpolation();
   //assert(checkDivergence());
}

void MACGrid::updateSources()
{
   /*
#ifdef _DEBUG
    mV(0,1,0) = 1.0;
    return;
#endif
    */
   // Set heat and smoke sources from bottom of container
   int size = 2; // neighborhood size
   int sourcei = (int) theDim[0]/5;
   int sourcek = (int) theDim[2]/2;

   for (int k = sourcek-size; k < sourcek+size; k++)
   {        
      mT(sourcei,0,k) = 2.0;
      mD(sourcei,0,k) = 2.0;
      mV(sourcei,0,k) = 1.0;   
   }
}

void MACGrid::advectVelocity(double dt)
{
   // Update each face
   FOR_EACH_FACE
   {
      if (isFace(i,j,k,X))
      {
         vec3 pos = getLeftFace(i,j,k);
         vec3 newpos = traceBack(pos, dt);
         vec3 newvel = getVelocity(newpos);
         target.mU(i, j, k) = newvel[X];
      }
      if (isFace(i,j,k,Y))
      {
         vec3 pos = getBottomFace(i,j,k);
         vec3 newpos = traceBack(pos, dt);
         vec3 newvel = getVelocity(newpos);
         target.mV(i, j, k) = newvel[Y];
      }
      if (isFace(i,j,k,Z))
      {
         vec3 pos = getBackFace(i,j,k);
         vec3 newpos = traceBack(pos, dt);
         vec3 newvel = getVelocity(newpos);
         target.mW(i, j, k) = newvel[Z];
      }
   }

   mU = target.mU;
   mV = target.mV;
   mW = target.mW;

}

void MACGrid::advectTemperature(double dt)
{
   FOR_EACH_CELL
   {
      vec3 pos = getCenter(i,j,k);
      vec3 newpos = traceBack(pos, dt);
      double newt = getTemperature(newpos);
      target.mT(i,j,k) = newt;
   }
   mT = target.mT;
}

void MACGrid::advectDensity(double dt)
{
   FOR_EACH_CELL
   {
      vec3 pos = getCenter(i,j,k);
      vec3 newpos = traceBack(pos, dt);
      double newd = getDensity(newpos);
      target.mD(i,j,k) = newd;
   }
   mD = target.mD;
}

double MACGrid::getBoussinesqForce(const vec3& pos)
{
   // Use Boussinesq approximation
   // f = [0, -alpha*smokeDensity + beta*(T - T_amb), 0]
   double temperature = getTemperature(pos); 
   double smokeDensity = getDensity(pos);   

   double yforce = -theBoussinesqAlpha*smokeDensity + 
      theBoussinesqBeta*(temperature - theAmbientTemp);

   return yforce;
}

void MACGrid::computeBouyancy(double dt)
{
   FOR_EACH_FACE
   {
      if (isFace(i,j,k,Y))
      {
         vec3 pos = getBottomFace(i,j,k);
         double yforce = getBoussinesqForce(pos);
         double vel = mV(i,j,k);
         vel = vel + dt*yforce;
         target.mV(i, j, k) = vel;
      }
   }
   mV = target.mV;
}

vec3 MACGrid::getVorticityN(int i, int j, int k)
{
   vec3 right  = getVorticity(i+1,j,k);
   vec3 left   = getVorticity(i-1,j,k);
   vec3 top    = getVorticity(i,j+1,k);
   vec3 bottom = getVorticity(i,j-1,k);
   vec3 front  = getVorticity(i,j,k+1);
   vec3 back   = getVorticity(i,j,k-1);

   double scale = 1.0/(2*theCellSize);
   double x = scale*(right.Length() - left.Length());
   double y = scale*(top.Length() - bottom.Length());
   double z = scale*(front.Length() - back.Length());

   vec3 N(x,y,z);
   return N.Normalize();
}

vec3 MACGrid::getVorticity(int i, int j, int k)
{
   vec3 right  = getVelocity(getRightFace(i+1,j,k));
   vec3 left   = getVelocity(getLeftFace(i,j,k));
   vec3 top    = getVelocity(getTopFace(i,j+1,k));
   vec3 bottom = getVelocity(getBottomFace(i,j,k));
   vec3 front  = getVelocity(getFrontFace(i,j,k+1));
   vec3 back   = getVelocity(getBackFace(i,j,k));

   double scale = 1.0/(theCellSize);
   double x = scale*(top[Z] - bottom[Z] - (front[Y] - back[Y]));
   double y = scale*(front[X] - back[X] - (right[Z] - left[Z]));
   double z = scale*(right[Y] - left[Y] - (top[X] - bottom[X]));

   return vec3(x, y, z);
}

vec3 MACGrid::getConfinementForce(int i, int j, int k)
{
   vec3 N = getVorticityN(i,j,k);
   vec3 w = getVorticity(i,j,k);
   return theVorticityConst*theCellSize*N.Cross(w);
}

void MACGrid::computeVorticityConfinement(double dt)
{
   GridData forcesX, forcesY, forcesZ;
   forcesX.initialize();
   forcesY.initialize();
   forcesZ.initialize();

   FOR_EACH_CELL // Calculate confinement forces
   {
      // TODO: Should I store this better?  Inverse interp?
      vec3 force = getConfinementForce(i,j,k);
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
         double vel = mU(i,j,k);
         double xforce = 0.5*(forcesX(i,j,k) - forcesX(i-1,j,k));
         vel = vel + dt*xforce;
         target.mU(i, j, k) = vel;
      }

      if (isFace(i,j,k,Y))
      {
         vec3 pos = getBottomFace(i,j,k);
         double yforce = 0.5*(forcesY(i,j,k) - forcesY(i,j-1,k));
         double vel = mV(i,j,k);
         vel = vel + dt*yforce;
         target.mV(i, j, k) = vel;
      }

      if (isFace(i,j,k,Z))
      {
         vec3 pos = getBackFace(i,j,k);
         double zforce = 0.5*(forcesZ(i,j,k) - forcesZ(i,j,k-1));
         double vel = mW(i,j,k);
         vel = vel + dt*zforce;
         target.mW(i, j, k) = vel;
      }
   }
   //std::cout << "XVel: " << mU.data() << std::endl;
   //std::cout << "YVel: " << mV.data() << std::endl;
   //std::cout << "ZVel: " << mW.data() << std::endl;

   mU = target.mU;
   mV = target.mV;
   mW = target.mW;
}


void MACGrid::addExternalForces(double dt)
{
   computeBouyancy(dt);
   if (theVConfEnabled) computeVorticityConfinement(dt);
}

double MACGrid::checkDivergence(int i, int j, int k)
{
   double x1 = mU(i+1,j,k);
   double x0 = mU(i,j,k);

   double y1 = mV(i,j+1,k);
   double y0 = mV(i,j,k);

   double z1 = mW(i,j,k+1);
   double z0 = mW(i,j,k);

   double xdiv = x1 - x0;
   double ydiv = y1 - y0;
   double zdiv = z1 - z0;
   double div = (xdiv + ydiv + zdiv)/theCellSize;
   return div;
}

double MACGrid::getDivergence(int i, int j, int k)
{
   // TODO: Can we fix this so that we calculate the true divergence
   // here and subtract arbitrary coundary conditions later?

   double x1 = isSolidCell(i+1, j, k)? 0.0 : mU(i+1,j,k);
   double x0 = isSolidCell(i-1, j, k)? 0.0 : mU(i,j,k);

   double y1 = isSolidCell(i, j+1, k)? 0.0 : mV(i,j+1,k);
   double y0 = isSolidCell(i, j-1, k)? 0.0 : mV(i,j,k);

   double z1 = isSolidCell(i, j, k+1)? 0.0 : mW(i,j,k+1);
   double z0 = isSolidCell(i, j, k-1)? 0.0 : mW(i,j,k);

   double xdiv = x1 - x0;
   double ydiv = y1 - y0;
   double zdiv = z1 - z0;
   double div = (xdiv + ydiv + zdiv)/theCellSize;

//   printf("Cell %d %d %d = %.2f\n", i, j, k, div);
   return div;
}

void MACGrid::constructB(ublas::vector<double>& b, unsigned int numCells, double dt)
{
   double constant = -(theAirDensity*theCellSize*theCellSize)/dt;
   for (unsigned int index = 0; index < numCells; index++)
   {
      int i,j,k; getCell(index, i, j, k);
      if (!isSolidCell(i,j,k))
         b(index) = constant*getDivergence(i, j, k);
      else b(index) = 0;
   }
#ifdef _DEBUG
   //std::cout << "B: " << b << std::endl;
#endif
}

#define VALA(r,c) (r != -1 && c != -1)? A(r,c) : 0
void MACGrid::constructPrecon(ublas::compressed_matrix<double>& A,
                              ublas::vector<double>& precon)
{
    precon.resize(A.size1());
    std::fill(precon.begin(), precon.end(), 0);

    double tau = 0.0; // Disable MIC(0) 0.97;
    for (unsigned int index = 0; index < A.size1(); index++)
    {
        int i,j,k; getCell(index, i, j, k);
        if (isSolidCell(i,j,k)) continue;

        int neighbori = getIndex(i-1,j,k);
        int neighborj = getIndex(i,j-1,k);
        int neighbork = getIndex(i,j,k-1);
        double termi = neighbori != -1? A(index,neighbori) * precon(neighbori) : 0;
        double termj = neighborj != -1? A(index,neighborj) * precon(neighborj) : 0;
        double termk = neighbork != -1? A(index,neighbork) * precon(neighbork) : 0;

        double termii = 0;
        if (neighbori != -1)
        {
            int neighborij = getIndex(i-1,j+1,k);
            int neighborik = getIndex(i-1,j,k+1);
            double termii0 = (VALA(neighbori,neighborij) + VALA(neighbori,neighborik));
            double termii1 = precon(neighbori)*precon(neighbori);
            termii = VALA(index,neighbori) * termii0 / termii1;
        }

        double termjj = 0;
        if (neighborj != -1)
        {
            int neighborji = getIndex(i+1,j-1,k);
            int neighborjk = getIndex(i,j-1,k+1);
            double termjj0 = (VALA(neighborj,neighborji) + VALA(neighborj,neighborjk));
            double termjj1 = precon(neighborj)*precon(neighborj);
            termjj = VALA(index,neighborj) * termjj0 / termjj1;
        }

        double termkk = 0;
        if (neighbork != -1)
        {
            int neighborki = getIndex(i+1,j,k-1);
            int neighborkj = getIndex(i,j+1,k-1);
            double termkk0 = (VALA(neighbork,neighborki) + VALA(neighbork,neighborkj));
            double termkk1 = precon(neighbork)*precon(neighbork);
            termkk = VALA(index,neighbork) * termkk0 / termkk1;
        }

        double e = A(index,index) - termi*termi - termj*termj - termk*termk
            - tau * (termii + termjj + termkk);

        precon(index) = 1/sqrt(e);
    }

#ifdef _DEBUG
    /*
   printf("E---------------------------\n");
   for (unsigned int i = 0; i < precon.size(); i++)
   {
      for (unsigned int j = 0; j < precon.size(); j++)
      {
         if (i == j) std::cout << precon(i) << " " ; 
         else std::cout << "0 " ;
      }
      std::cout  << std::endl; 
   }*/
#endif
}

void MACGrid::constructA(ublas::compressed_matrix<double>& A)
{
   unsigned int numCells = mSolid.data().size();
   A.resize(numCells, numCells, false);

   for (unsigned int row = 0; row < numCells; row++)
   {
      //std::cout << "Initializing row " << row << "/" << numCells << std::endl;
      int ri, rj, rk; getCell(row, ri, rj, rk); // Each row corresponds to a cell
      if (isSolidCell(ri, rj, rk)) continue;

      for (unsigned int col = 0; col < numCells; col++)
      {         
         int ci, cj, ck; getCell(col, ci, cj, ck); // Each col corresponds to a possible neighbor
         if (isSolidCell(ci, cj, ck)) continue;
         double coeff = getPressureCoeffBetweenCells(ri, rj, rk, ci, cj, ck);
         if (fabs(coeff) > 0.0001) 
         {
            A(row,col) = coeff;         
         }
      }
   } 

#ifdef _DEBUG
   /*
   printf("A---------------------------\n");
   for (unsigned int i = 0; i < A.size1(); i++)
   {
      for (unsigned int j = 0; j < A.size2(); j++)
      {
         int test = A(i, j);
         std::cout << A(i, j) << " " ; 
      }
      std::cout  << std::endl; 
   }*/
#endif
}

void MACGrid::project(double dt)
{
   // Solve Ax = b for pressure
   unsigned int numCells = theDim[0]*theDim[1]*theDim[2];

   // 1. Contruct b
   ublas::vector<double> b(numCells);
   constructB(b, numCells, dt);

   // 2. Construct A - currently cached

   // 3. Solve for p
   ublas::vector<double> p(numCells);
   //cg_solve(A, b, p, 500, 0.005);
   cg_psolve(A, precon, b, p, 500, 0.005);
   //std::cout << "P: " << p << std::endl;

   // Subtract pressure from our velocity and save in target
   // u_new = u - dt*(1/theAirPressure)*((p_i+1-p_i)/theCellSize)
   double scaleConstant = dt/theAirDensity;
   double pressureChange; 

   FOR_EACH_FACE
   {
      if (isFace(i,j,k,X))
      {
         if (isSolidFace(i, j, k, X)) 
         {
            target.mU(i, j, k) = 0.0;
         }
         else
         {
            int index1 = getIndex(i,j,k);
            int index2 = getIndex(i-1,j,k);
            pressureChange = (p(index1) - p(index2))/theCellSize;
            double vel = mU(i,j,k);            
            vel = vel - scaleConstant*pressureChange;
            target.mU(i, j, k) = vel;
         }           
         
      }
      if (isFace(i,j,k,Y))   
      {
         // Hard-code boundary condition for now
         if (isSolidFace(i,j,k,Y))
         {
            target.mV(i, j, k) = 0.0;
         }
         else
         {
            int index1 = getIndex(i,j,k);
            int index2 = getIndex(i,j-1,k);
            pressureChange = (p(index1) - p(index2))/theCellSize;
            double vel = mV(i,j,k);
            vel = vel - scaleConstant*pressureChange;
            target.mV(i, j, k) = vel;
         }
      }

      if (isFace(i,j,k,Z))
      {
         // Hard-code boundary condition for now
         if (isSolidFace(i,j,k,Z))
         {
            target.mW(i, j, k) = 0.0;
         }
         else
         {
            int index1 = getIndex(i,j,k);
            int index2 = getIndex(i,j,k-1);
            pressureChange = (p(index1) - p(index2))/theCellSize;
            double vel = mW(i,j,k);
            vel = vel - scaleConstant*pressureChange;
            target.mW(i, j, k) = vel;
         }
      }
   }

   //std::cout << "u: " << target.mU.data() << std::endl;
   //std::cout << "v: " << target.mV.data() << std::endl;
   //std::cout << "w: " << target.mW.data() << std::endl;

   mU = target.mU;
   mV = target.mV;
   mW = target.mW;
   assert (checkDivergence());
}

bool MACGrid::checkDivergence()
{
   FOR_EACH_CELL
   {
      double div = checkDivergence(i, j, k);
      if (fabs(div) > 0.01) 
      {
         printf("Divergence(%d,%d,%d) = %.2f\n", i, j, k, div);
         return false;
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

int MACGrid::getIndex(int i, int j, int k)
{
   if (i < 0 || i > theDim[0]-1) return -1;
   if (j < 0 || j > theDim[1]-1) return -1;
   if (k < 0 || k > theDim[2]-1) return -1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];
   return col+row+stack;
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
   vel[0] = getVelocityX(pt); 
   vel[1] = getVelocityY(pt); 
   vel[2] = getVelocityZ(pt); 
   return vel;
}

double MACGrid::getVelocityX(const vec3& pt)
{
   return mU.interpolate(pt);
}

double MACGrid::getVelocityY(const vec3& pt)
{
   return mV.interpolate(pt);
}

double MACGrid::getVelocityZ(const vec3& pt)
{
   return mW.interpolate(pt);
}

double MACGrid::getTemperature(const vec3& pt)
{
   return mT.interpolate(pt);
}

double MACGrid::getDensity(const vec3& pt)
{
   return mD.interpolate(pt);
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

void MACGrid::draw(const Camera& c)
{   
   drawWireGrid();
   if (theDisplayVel) drawVelocities();   
   drawSolids(c);
   if (theRenderMode == CUBES) drawSmokeCubes(c);
   else drawSmoke(c);
}

void MACGrid::drawVelocities()
{
   // draw line at each center
   glColor4f(0.0, 1.0, 0.0, 1.0);
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
           glVertex3dv(pos.n);
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
         vec3 force = getConfinementForce(i,j,k);
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

vec4 MACGrid::getRenderColor(int i, int j, int k)
{
    double value = mD(i, j, k); 
    return vec4(1.0, 0.9, 1.0, value);
}

vec4 MACGrid::getRenderColor(const vec3& pt)
{
    double value = getDensity(pt); 
    return vec4(1.0, 1.0, 1.0, value);
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
      //cube.color.set(1.0, 0.4, 0.8, 1.0);
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
      if (isFace(i,j,k,X)) mU(i, j, k) = i*theCellSize;
      if (isFace(i,j,k,Y)) mV(i, j, k) = j*theCellSize;
      if (isFace(i,j,k,Z)) mW(i, j, k) = k*theCellSize; 
   }

   FOR_EACH_CELL
   {
      mD(i, j, k) = i + j + k;
   }
   //std::cout << "D: " << mD.data() << std::endl;

   vec3 test = getVelocity(vec3(0.25,0.25,0.2));  
   vec3 test1 = getVelocity(vec3(1.5,2.2,0.1));  
   vec3 test2 = getVelocity(vec3(0.0,0.1,0.1));  

   double test3 = getDensity(vec3(0.0, 0.0, 0.0)); // should be 0
   double test4 = getDensity(vec3(0.5, 0.5, 0.0));  // should be 0
   double test5 = getDensity(vec3(1.0, 1.0, 0.0));  // should be 1
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