#include <stdio.h>
#include <math.h>
#include "fireSim.h"
#include "camera.h"
#include "fps.h"
#include <GL/glut.h>
#include <math.h>

// Geometry and whatnot
FireSim theSmokeSim;
Camera theCamera;
mmc::FpsTracker theFpsTracker;

// UI Helpers
int lastX = 0, lastY = 0;
int theMenu = 0;
int theButtonState = 0;
int theModifierState = 0;
bool isRunning = true;

void initCamera()
{
   extern int theDim[3]; // Naughty globals...
   extern double theCellSize;
   double w = theDim[0]*theCellSize;   
   double h = theDim[1]*theCellSize;   
   double d = theDim[2]*theCellSize;   
   double angle = 0.5*theCamera.dfltVfov*M_PI/180.0;
   double dist;
   if (w > h) dist = w*0.5/tan(angle);  // aspect is 1, so i can do this
   else dist = h*0.5/tan(angle);
   theCamera.dfltEye.set(w*0.5, h*0.5, -(dist+d*0.5));
   theCamera.dfltLook.set(w*0.5, h*0.5, 0.0);
   theCamera.reset();
}

void onMouseMotionCb(int x, int y)
{
   int deltaX = lastX - x;
   int deltaY = lastY - y;
   bool moveLeftRight = abs(deltaX) > abs(deltaY);
   bool moveUpDown = !moveLeftRight;

   if (theButtonState == GLUT_LEFT_BUTTON)  // Rotate
   {
      if (moveLeftRight && deltaX > 0) theCamera.orbitLeft(deltaX);
      else if (moveLeftRight && deltaX < 0) theCamera.orbitRight(-deltaX);
      else if (moveUpDown && deltaY > 0) theCamera.orbitUp(deltaY);
      else if (moveUpDown && deltaY < 0) theCamera.orbitDown(-deltaY);
   }
   else if (theButtonState == GLUT_MIDDLE_BUTTON) // Zoom
   {
      if (moveUpDown && deltaY > 0) theCamera.moveForward(deltaY);
      else if (moveUpDown && deltaY < 0) theCamera.moveBack(-deltaY);
   }    

   if (theModifierState & GLUT_ACTIVE_ALT) // camera move
   {
      if (theButtonState == GLUT_RIGHT_BUTTON) // Pan
      {
         if (moveLeftRight && deltaX > 0) theCamera.moveLeft(deltaX);
         else if (moveLeftRight && deltaX < 0) theCamera.moveRight(-deltaX);
         else if (moveUpDown && deltaY > 0) theCamera.moveUp(deltaY);
         else if (moveUpDown && deltaY < 0) theCamera.moveDown(-deltaY);
      }   
   }
 
   lastX = x;
   lastY = y;
   glutPostRedisplay();
}

void onMouseCb(int button, int state, int x, int y)
{
   theButtonState = button;
   theModifierState = glutGetModifiers();
   lastX = x;
   lastY = y;

   glutSetMenu(theMenu);
   if (theModifierState & GLUT_ACTIVE_ALT)
   {
      glutDetachMenu(GLUT_RIGHT_BUTTON);
   }
   else
   {
      glutAttachMenu(GLUT_RIGHT_BUTTON);
   }

   onMouseMotionCb(x, y);
}


void onKeyboardCb(unsigned char key, int x, int y)
{
   if (key == ' ') theCamera.reset();
   else if (key == '0') MACGrid::theRenderMode = NONE;
   else if (key == '1') MACGrid::theRenderMode = SMOKE;
   else if (key == '2') MACGrid::theRenderMode = TEMPERATURE;
   else if (key == '3') MACGrid::theRenderMode = PRESSURE;
   else if (key == '4') MACGrid::theRenderMode = LEVELSET;
   else if (key == '5') MACGrid::theRenderMode = BLUECORE;
   else if (key == '6') MACGrid::theRenderMode = ALL;
   else if (key == 'v') MACGrid::theDisplayVel = !MACGrid::theDisplayVel;
   else if (key == 'f') MACGrid::theDisplayVForces = !MACGrid::theDisplayVForces;
   else if (key == 'g') MACGrid::theDisplayGrid = !MACGrid::theDisplayGrid;
   else if (key == 's') MACGrid::theSourceEnabled = !MACGrid::theSourceEnabled;
   else if (key == 'r') theSmokeSim.setRecording(!theSmokeSim.isRecording());
   else if (key == '>') isRunning = true;
   else if (key == '=') isRunning = false;
   else if (key == '<') theSmokeSim.reset();
   else if (key == '.') theSmokeSim.step();
   else if (key == 27) exit(0); // ESC Key
   glutPostRedisplay();
}

void onMenuCb(int value)
{
   switch (value)
   {
   case -1: exit(0);
   case -2: theSmokeSim.setGridDimensions(40,30,2); initCamera(); break;
   case -3: theSmokeSim.setGridDimensions(10,20,10); initCamera(); break;
   case -4: MACGrid::theVConfEnabled = !MACGrid::theVConfEnabled; break; 
   case -5: MACGrid::thePreconEnabled = !MACGrid::thePreconEnabled; break; 
   case -6: theSmokeSim.reset(SNONE); break;
   default: onKeyboardCb(value, 0, 0); break;
   }
}

void onKeyboardSpecialCb(int key, int x, int y)
{
}

void onTimerCb(int value)
{
   if (isRunning) theSmokeSim.step();
   glutTimerFunc(100, onTimerCb, 0);
   glutPostRedisplay();
}

void onResizeCb(int width, int height)
{
   // Update viewport
   glViewport(0, 0, width, height);

   // Update camera projection's aspect ratio
   float vfov, aspect, zNear, zFar;
   theCamera.getProjection(&vfov, &aspect, &zNear, &zFar);
   theCamera.setProjection(vfov, ((GLfloat) width)/height, zNear, zFar);
}

void drawOverlay()
{
  // Draw Overlay
  glColor4f(1.0, 1.0, 1.0, 1.0);
  glPushAttrib(GL_LIGHTING_BIT);
     glDisable(GL_LIGHTING);

     glMatrixMode(GL_PROJECTION);
     glLoadIdentity();
     gluOrtho2D(0.0, 1.0, 0.0, 1.0);

     glMatrixMode(GL_MODELVIEW);
     glLoadIdentity();
     glRasterPos2f(0.01, 0.01);
     
     char* modeString;
     switch (MACGrid::theRenderMode)
     {
     case NONE: modeString = "None"; break; 
     case SMOKE: modeString = "Smoke"; break; 
     case TEMPERATURE: modeString = "Temperature"; break; 
     case PRESSURE: modeString = "Pressure"; break; 
     case LEVELSET: modeString = "Levelset"; break;
     case BLUECORE: modeString = "BlueCore"; break;
     case ALL: modeString = "Fire"; break;
     }

     char info[1024];
	 sprintf(info, "Framerate: %3.1f   Step: %i   Display: %s   Sources: %s  %s", 
         theFpsTracker.fpsAverage(), 
         theSmokeSim.getStepNum(),
         modeString,
         MACGrid::theSourceEnabled? "On" : "Off",
         theSmokeSim.isRecording()? "Recording..." : "");
 
     for (unsigned int i = 0; i < strlen(info); i++)
     {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, info[i]);
     }
  glPopAttrib();
}

void onDrawCb()
{
  // Keep track of time
  theFpsTracker.timestamp();

  // Draw Scene and overlay
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  theCamera.draw();
  theSmokeSim.draw(theCamera);
  drawOverlay();
  glutSwapBuffers();
}

void init(void)
{
   initCamera(); 
	glClearColor(0.0, 0.0, 0.0, 1.0); // Or 0.5, 0.5, 0.5 to show off black smoke

   glEnable(GL_BLEND);
   glEnable(GL_ALPHA_TEST);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glEnable(GL_DEPTH_TEST);
   glDepthFunc(GL_LEQUAL);
	glShadeModel(GL_SMOOTH);

   glEnable(GL_NORMALIZE);
   glDisable(GL_LIGHTING);
   glCullFace(GL_BACK);

   theSmokeSim.setRecording(true);
}

int main(int argc, char **argv)
{
   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(SCREEN_WIDTH, SCREEN_HEIGHT);
	glutInitWindowPosition(100, 100);
   glutCreateWindow("FireSim by TravisG & AlineN");
   glutDisplayFunc(onDrawCb);
   glutKeyboardFunc(onKeyboardCb);
   glutSpecialFunc(onKeyboardSpecialCb);
   glutMouseFunc(onMouseCb);
   glutMotionFunc(onMouseMotionCb); 
   glutTimerFunc(100, onTimerCb, 0);
   glutReshapeFunc(onResizeCb);

   int viewMenu = glutCreateMenu(onMenuCb);
   glutAddMenuEntry("Display smoke\t'1'", '1');
   glutAddMenuEntry("Display temperature\t'2'", '2');
   glutAddMenuEntry("Display pressure\t'3'", '3');
   glutAddMenuEntry("Toggle velocities\t'v'", 'v');
   glutAddMenuEntry("Toggle forces\t'f'", 'f');
   glutAddMenuEntry("Toggle grid\t'g'", 'g');

   int gridMenu = glutCreateMenu(onMenuCb);
   glutAddMenuEntry("Make 40x30x2", -2); 
   glutAddMenuEntry("Make 10x20x10", -3);

   int simMenu = glutCreateMenu(onMenuCb);
   glutAddMenuEntry("Toggle vorticity confinement", -4);
   glutAddMenuEntry("Toggle preconditioner", -5);
   glutAddMenuEntry("Hide solid", -6);
   glutAddMenuEntry("Toggle sources\t's'", 's');

   theMenu = glutCreateMenu(onMenuCb);
   glutAddMenuEntry("Start\t'>'", '>');
   glutAddMenuEntry("Pause\t'='", '=');
   glutAddMenuEntry("Reset\t'<'", '<');
   glutAddMenuEntry("Reset camera\t' '", ' ');
   glutAddMenuEntry("Record\t'r'", 'r');
   glutAddSubMenu("Display", viewMenu);
   glutAddSubMenu("Grids", gridMenu);
   glutAddSubMenu("Simulation", simMenu);
   glutAddMenuEntry("_________________", -1);
   glutAddMenuEntry("Exit", 27);
   glutAttachMenu(GLUT_RIGHT_BUTTON);

   init();

   glutMainLoop();
   return 0;             
}

