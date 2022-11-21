#include "smokeSim.h"
#include <GL/glut.h>

SmokeSim::SmokeSim() : mFrameNum(0), mRecordEnabled(false)
{
   reset();
}

SmokeSim::~SmokeSim()
{
}

void SmokeSim::reset()
{
   mGrid.reset();
}

void SmokeSim::setGridDimensions(int x, int y, int z)
{
   extern int theDim[3]; // Naughty globals...
   theDim[0] = x;
   theDim[1] = y;
   theDim[2] = z;
   reset();
}

void SmokeSim::step()
{
   double dt = 0.01;

   // Step0: Gather user forces
   mGrid.updateSources();

   // Step1: Calculate new velocities
   mGrid.advectVelocity(dt);
   mGrid.addExternalForces(dt);
   mGrid.project(dt);

   // Step2: Calculate new temperature
   mGrid.advectTemperature(dt);

   // Step3: Calculate new density 
   mGrid.advectDensity(dt);
}

void SmokeSim::setRecording(bool on)
{
   if (on && ! mRecordEnabled)  // reset counter
   {
      mFrameNum = 0;
   }
   mRecordEnabled = on;
}

bool SmokeSim::isRecording()
{
   return mRecordEnabled;
}

void SmokeSim::draw(const Camera& c)
{
   drawAxes(); 
   mGrid.draw(c);
}

void SmokeSim::drawAxes()
{
  glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
      glDisable(GL_LIGHTING);

      glLineWidth(2.0); 
      glBegin(GL_LINES);
         glColor3f(1.0, 0.0, 0.0);
         glVertex3f(0.0, 0.0, 0.0);
         glVertex3f(1.0, 0.0, 0.0);

         glColor3f(0.0, 1.0, 0.0);
         glVertex3f(0.0, 0.0, 0.0);
         glVertex3f(0.0, 1.0, 0.0);

         glColor3f(0.0, 0.0, 1.0);
         glVertex3f(0.0, 0.0, 0.0);
         glVertex3f(0.0, 0.0, 1.0);
      glEnd();
  glPopAttrib();
}
