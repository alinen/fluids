#include <stdio.h>
#include <math.h>
#include "camera.h"
#include "fps.h"
#include "LevelSet.h"
#include "container.h"
#include <GL/glut.h>
#include <math.h>

// Geometry and whatnot
int theDim[3] = {10, 10, 10}; 
double theCellSize = 0.5;
Camera theCamera;
mmc::FpsTracker theFpsTracker;
Container contain(NX, NY, NZ, theCellSize);

// UI Helpers
int lastX = 0, lastY = 0;
int theMenu = 0;
int theButtonState = 0;
int theModifierState = 0;
bool isRunning = false;

void initCamera()
{
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
   else if (key == '>') isRunning = true;
   else if (key == '=') isRunning = false;
   else if (key == 27) exit(0); // ESC Key
   glutPostRedisplay();
}

void onMenuCb(int value)
{
   switch (value)
   {
   case -1: exit(0);
   default: onKeyboardCb(value, 0, 0); break;
   }
}

void onKeyboardSpecialCb(int key, int x, int y)
{
}

void onTimerCb(int value)
{
   if (isRunning)
   {
       contain.Update();
   }
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
     
     char info[1024];
     sprintf(info, "Framerate: %3.1f", 
         theFpsTracker.fpsAverage());
 
     for (unsigned int i = 0; i < strlen(info); i++)
     {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, info[i]);
     }
  glPopAttrib();
}

void drawWireGrid()
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

void drawLevelset()
{
    glBegin(GL_LINES);
    FOR_LS
        // Draw normals outside the surface for one sheet
        if (k != NZ/2) continue;

        double dist = contain.lset.LinearSample(Vector(i,j,k));
        if (dist < 0.0) continue;

		Vector pos(i * theCellSize, j * theCellSize, k * theCellSize);  // Grid position in world pos
        double distx = contain.lset.LinearSample(Vector(i+1,j,k)) - contain.lset.LinearSample(Vector(i-1,j,k));;
        double disty = contain.lset.LinearSample(Vector(i,j+1,k)) - contain.lset.LinearSample(Vector(i,j-1,k));
        double distz = contain.lset.LinearSample(Vector(i,j,k+1)) - contain.lset.LinearSample(Vector(i,j,k-1));
        Vector normal(distx, disty, distz);
        normal = (0.5*theCellSize) * normal;
        normal.Normalize();
        glColor3f(0.0, 0.0, 1.0);
        glVertex3f(pos[0], pos[1], pos[2]);
        glColor3f(1.0, 0.0, 0.0);
        glVertex3f(pos[0]+0.25*normal[0], pos[1]+0.25*normal[1], pos[2]+0.25*normal[2]);
    END_FOR_THREE
    glEnd();

    glColor3f(1.0, 1.0, 1.0);
    glPointSize(2.0);
    glBegin(GL_POINTS);
    FOR_LS
        double dist = contain.lset.LinearSample(Vector(i,j,k));
        // Draw points in the interior of the level set
		Vector pos(i * theCellSize, j * theCellSize, k * theCellSize);  // Grid position in world pos
        if (dist < 0.1) glVertex3f(pos[0], pos[1], pos[2]);
    END_FOR_THREE
    glEnd();
}

void drawScene()
{
    drawWireGrid();
    drawLevelset();

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
}

void onDrawCb()
{
  // Keep track of time
  theFpsTracker.timestamp();

  // Draw Scene and overlay
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  theCamera.draw();
  drawScene();
  drawOverlay();
  glutSwapBuffers();
}

void init(void)
{
    initCamera();
    glClearColor(0.1, 0.1, 0.1, 1.0);

    glEnable(GL_BLEND);
    glEnable(GL_ALPHA_TEST);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glShadeModel(GL_SMOOTH);

    glEnable(GL_NORMALIZE);
    glDisable(GL_LIGHTING);
    glCullFace(GL_BACK);
}


int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    glutInitWindowSize(640, 480);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("LevelSetTest - CIS563");
    glutDisplayFunc(onDrawCb);
    glutKeyboardFunc(onKeyboardCb);
    glutSpecialFunc(onKeyboardSpecialCb);
    glutMouseFunc(onMouseCb);
    glutMotionFunc(onMouseMotionCb); 
    glutTimerFunc(100, onTimerCb, 0);
    glutReshapeFunc(onResizeCb);

    theMenu = glutCreateMenu(onMenuCb);
    glutAddMenuEntry("Start\t'>'", '>');
    glutAddMenuEntry("Pause\t'='", '=');
    glutAddMenuEntry("Reset camera\t' '", ' ');
    glutAddMenuEntry("_________________", -1);
    glutAddMenuEntry("Exit", 27);
    glutAttachMenu(GLUT_RIGHT_BUTTON);

    init();

    glutMainLoop();
    return 0;             
}

