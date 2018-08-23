#include "GL/glut.h"
#include "main.h"
#include "container2d.h"
#include "timer.h"
#include "particle2d.h"
#include "particleset2d.h"
#include "vec2d.h"

/* global variables */

Container2D container(NX, NY, HH, DT);

//static int win_x = 6*NX, win_y = 6*NY;
static int win_x = 600, win_y = 600;
static int dvel = 3, pause = 1, area = 0, rate = 0;

static int win_id;

static int mouse_down[3];
static int omx, omy, mx, my;
Double dt = 0;

static Double rot_time = 0.0;
static Double fps_time = 0.0;
static int frames = 0;
Timer timer;

/*
  ----------------------------------------------------------------------
   OpenGL specific drawing routines
  ----------------------------------------------------------------------
*/

static void pre_display ( void )
{
	glViewport ( 0, 0, win_x, win_y );
	glMatrixMode ( GL_PROJECTION );
	glLoadIdentity ();
	//gluOrtho2D ( 0.0, 1.0, 0.0, 1.0 );
	gluOrtho2D ( HH, (HH*NX), HH, (HH*NY) );
	glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
	glClear ( GL_COLOR_BUFFER_BIT );
}

static void post_display ( void )
{
	glutSwapBuffers ();
}

static void draw_velocity ( void )
{
	static Float x, y, u, v;

	glColor3f ( 1.0f, 1.0f, 1.0f );
	glLineWidth ( 1.0f );

	glBegin ( GL_LINES );

		for (int i=1 ; i<=NX ; i++ ) {
			x = (i-0.5f)*HH;
			for (int j=1 ; j<=NY ; j++ ) {
				y = (j-0.5f)*HH;
				container.GetVelocity(i,j,u,v);
				glVertex2f ( x, y );
				glVertex2f ( x+u, y+v );
			}
		}

	glEnd ();
}

static void draw_density ( void )
{
	static Float x, y, d00, d01, d10, d11;
	glBegin ( GL_QUADS );

		for (int i=1 ; i<=NX ; i++ ) {
			x = (i-0.5)*HH;
			for (int j=1 ; j<=NY ; j++ ) {
				y = (j-0.5)*HH;

                d00 = container(i,j)     > 0. ? 0. : 200;
				d01 = container(i,j+1)   > 0. ? 0. : 200;
				d10 = container(i+1,j)   > 0. ? 0. : 200;
				d11 = container(i+1,j+1) > 0. ? 0. : 200;

				glColor3f ( d00, d00, d00 ); glVertex2f ( x, y );
				glColor3f ( d10, d10, d10 ); glVertex2f ( x+HH, y );
				glColor3f ( d11, d11, d11 ); glVertex2f ( x+HH, y+HH );
				glColor3f ( d01, d01, d01 ); glVertex2f ( x, y+HH );
			}
		}

	glEnd ();
}

static void draw_points ( void )
{
	static Float x, y, d00;
	glBegin ( GL_POINTS );

		for (int i=1 ; i<=NX ; i++ ) {
			x = (i-0.5)*HH;
			for (int j=1 ; j<=NY ; j++ ) {
				y = (j-0.5)*HH;
                d00 = container(i,j)     > 0. ? 0. : 1.0;
				glColor3f ( d00, d00, d00 ); glVertex2f ( x, y );
			}
		}

	glEnd ();
}

static void FindArea()
{
    int a = 0;
	for (int i=1 ; i<=NX ; i++ ) 
		for (int j=1 ; j<=NY ; j++ ) 
            a += container(i,j) > 0. ? 0 : 1;
    cout << endl << endl << "Area: " << a << endl << endl;
    area = 0;
}

static void draw_particles (void)
{
    static Float x,y;
    glBegin ( GL_POINTS );

        ParticleSet2D::cIterator pit, end = container.pset.end();
        for(pit = container.pset.begin(); pit != end; ++pit) {
            Vec2D pos;
            Particle2D& particle = *(*pit);
		    particle.GetPosition(pos);
		    Double phi = container.lset.LinearSample(pos);
		    Double sign = particle.Sign();

            x = (pos[0] - .5) * HH;
            y = (pos[1] - .5) * HH;

            if(sign * phi < 0) {
                if(sign > 0) glColor3f ( 0.0, 1.0, 0.0 );
                else glColor3f ( 1.0, 1.0, 0.0 );
            }
            else {
                if(sign > 0) glColor3f ( 1.0, 0.0, 0.0 );
                else glColor3f ( 0.0, 0.0, 1.0 );
            }
            glVertex2f ( x, y );
        }

    glEnd();
}

/*
  ----------------------------------------------------------------------
   GLUT callback routines
  ----------------------------------------------------------------------
*/

static void key_func ( unsigned char key, int x, int y )
{
	switch ( key )
	{
		case 'c':
		case 'C':
			container.Clear();
			break;

		case 'q':
		case 'Q':
			exit ( 0 );
			break;

		case 'x':
		case 'X':
			dvel = dvel++ % 4;
			break;
        case 'p':
        case 'P':
            pause = pause ? 0 : 1;
            break;
        case 'a':
        case 'A':
            area = 1;
            break;
        case 'f':
        case 'F':
            rate = rate ? 0 : 1;
            break;
	}
}

static void reshape_func ( int width, int height )
{
	glutSetWindow ( win_id );
	glutReshapeWindow ( width, height );

	win_x = width;
	win_y = height;
}

static void idle_func ( void )
{
    if(!pause) { container.Update(); pause = 0; dt += DT; }
    if(dt > (628 * HH)) { pause = 0; dt = 0; cout << rot_time << endl; rot_time = 0; }
    if(area) FindArea();

	glutSetWindow ( win_id );
	glutPostRedisplay ();
}

static void display_func ( void )
{	
	fps_time += timer.GetElapsedTime();
    rot_time += timer.GetElapsedTime();
	timer.Reset();
	frames++;

	if(fps_time >= 10.0 && rate)
	{
		cout<<"Framerate:"<<(frames/fps_time)<<endl;
		frames = 0;
		fps_time = 0;
	}
	pre_display ();

		if      ( dvel==0 ) draw_velocity ();
		else if	( dvel==1 )	draw_density ();
        else if ( dvel==2 ) draw_points ();
        else if ( dvel==3 ) draw_particles ();

	post_display ();
}


/*
  ----------------------------------------------------------------------
   open_glut_window --- open a glut compatible window and set callbacks
  ----------------------------------------------------------------------
*/

static void open_glut_window ( void )
{
	glutInitDisplayMode ( GLUT_RGBA | GLUT_DOUBLE );

	glutInitWindowPosition ( 0, 0 );
	glutInitWindowSize ( win_x, win_y );
	win_id = glutCreateWindow ( "2D Level Set Library" );

	glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
	glClear ( GL_COLOR_BUFFER_BIT );
	glutSwapBuffers ();
	glClear ( GL_COLOR_BUFFER_BIT );
	glutSwapBuffers ();

	pre_display ();

	glutKeyboardFunc ( key_func );
	glutReshapeFunc ( reshape_func );
	glutIdleFunc ( idle_func );
	glutDisplayFunc ( display_func );
}


/*
  ----------------------------------------------------------------------
   main --- main routine
  ----------------------------------------------------------------------
*/

int main ( int argc, char ** argv )
{
	glutInit ( &argc, argv );

	printf ( "\n\nHow to use this demo:\n\n" );
	printf ( "\t Start and Stop the demo with the 'p' key\n" );
	printf ( "\t Toggle density/velocity display with the 'x' key\n" );
	printf ( "\t Clear the simulation by pressing the 'c' key\n" );
	printf ( "\t Quit by pressing the 'q' key\n" );

	open_glut_window ();

	timer.Reset();
	glutMainLoop ();

	exit ( 0 );
}