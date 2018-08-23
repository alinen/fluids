#include "fireSim.h"
#include "renderman.h"
#include <GL/glut.h>
#include <IL/il.h>
#include <IL/ilu.h>
#include <IL/ilut.h>

FireSim::FireSim() : mFrameNum(0), mStepNum(0), mRecordEnabled(false)
{
   ilInit();
   iluInit();
   ilEnable(IL_FILE_OVERWRITE);
   ilutRenderer(ILUT_OPENGL);

   reset();
}

FireSim::~FireSim()
{
}

void FireSim::reset()
{
   mGrid.reset();
   extern MACGrid target;  // HACK - laziness!
   target.reset();   
}

void FireSim::reset(SolidMode m)
{
   mGrid.reset(m);
   extern MACGrid target;  // HACK - laziness!
   target.reset(m);   
}

void FireSim::setGridDimensions(int x, int y, int z)
{
   extern int theDim[3]; // Naughty globals...
   theDim[0] = x;
   theDim[1] = y;
   theDim[2] = z;
   reset(SNONE);
}

void FireSim::step()
{
   double dt = 0.01;
   mStepNum++;

   // Step0: Gather user forces
   mGrid.updateSources();

   // Step1: Evolve flame front
   mGrid.evolveFlameFront(dt);

   // Step2: Calculate new velocities
   mGrid.advectVelocity(dt);
   mGrid.addExternalForces(dt);
   mGrid.project(dt);

   // Step3: Calculate new temperature
   mGrid.advectTemperature(dt);

   // Step4: Calculate new density and color
   mGrid.advectSmokeDensity(dt);   
}

void FireSim::setRecording(bool on)
{
   if (on && ! mRecordEnabled)  // reset counter
   {
      mFrameNum = 0;
   }
   mRecordEnabled = on;
}

bool FireSim::isRecording()
{
   return mRecordEnabled;
}

void FireSim::draw(const Camera& c)
{
   drawAxes(); 
   mGrid.draw(c);
   if (mRecordEnabled) 
   {
	   grabScreen();
	   generateRib();
      mFrameNum++;
   }

}

void FireSim::drawAxes()
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

void FireSim::grabScreen()  // Code adapted from asst#1
{
	unsigned int image;
    ilGenImages(1, &image);
	ilBindImage(image);

	ILenum error = ilGetError();
	assert(error == IL_NO_ERROR);

	ilTexImage(SCREEN_WIDTH, SCREEN_HEIGHT, 1, 3, IL_RGB, IL_UNSIGNED_BYTE, NULL);   

	error = ilGetError();
	assert(error == IL_NO_ERROR);

	unsigned char* data = ilGetData();

	error = ilGetError();
	assert(error == IL_NO_ERROR);

	for (int i=SCREEN_HEIGHT-1; i>=0; i--) 
	{
		glReadPixels(0,i,SCREEN_WIDTH,1,GL_RGB, GL_UNSIGNED_BYTE, 
			data + (SCREEN_WIDTH * 3 * i));
	}

	char anim_filename[2048];
	sprintf_s(anim_filename, 2048, "output/hwfire_%04d.png", mFrameNum); 

	ilSave(IL_PNG, anim_filename);

	error = ilGetError();
	assert(error == IL_NO_ERROR);

	ilDeleteImages(1, &image);

	error = ilGetError();
	assert(error == IL_NO_ERROR);
}

void FireSim::generateRib()
{
#if HAVE_PRMAN
#ifndef _DEBUG
	char rib_filename[2048];
	char image_filename[2048];

	sprintf_s(rib_filename, 2048, "output/fireblob_%04d.rib", mFrameNum);
	sprintf_s(image_filename, 2048, "fireblob_%04d.png", mFrameNum);
	render(rib_filename, image_filename, mGrid, true);

	sprintf_s(rib_filename, 2048, "output/fireparticle_%04d.rib", mFrameNum);
	sprintf_s(image_filename, 2048, "fireparticle_%04d.png", mFrameNum);
	render(rib_filename, image_filename, mGrid, false);
#endif
#endif
}