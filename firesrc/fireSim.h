#ifndef FireSim_H_
#define FireSim_H_

#include "MACGrid.h"
#define SCREEN_WIDTH 640
#define SCREEN_HEIGHT 480

class Camera;
class FireSim
{
public:
   FireSim();
   virtual ~FireSim();

   virtual void reset();
   virtual void reset(SolidMode m);
   virtual void step();
   virtual void draw(const Camera& c);
   virtual void setGridDimensions(int x, int y, int z);
   virtual void setRecording(bool on);
   virtual bool isRecording();
   int getStepNum(){return mStepNum;}

protected:
   virtual void drawAxes();
   virtual void grabScreen();
   virtual void generateRib();

protected:
   MACGrid mGrid;
   bool mRecordEnabled;
   int mFrameNum;
   int mStepNum;
};

#endif