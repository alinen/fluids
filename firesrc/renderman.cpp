/* Copyrighted Pixar 1989 */
/* From the RenderMan Companion p. 168 */
/* Listing 8.5 An improved boilerplate viewing program   */
#if HAVE_PRMAN

#include <ri.h>
#include <stdio.h>
#include <math.h>
#include "renderman.h"

//link against  libprmansdk.lib prman.lib 

// External simulation globals
extern int theDim[3];
extern double theCellSize;

// Picture and camera settings
#define DATATYPE RI_RGBA         /* Pixels have RGB and coverage */
#define PICXRES 640//1280.0            /* Horizontal output resolution */
#define PICYRES 480//720.0            /* Vertical output resolution   */

#define CROPMINX 0.0             /* RiCropWindow() parameters    */
#define CROPMAXX  1.0 
#define CROPMINY  0.0
#define CROPMAXY  1.0 
#define CAMXFROM 0.0             /* Camera position  */
#define CAMYFROM 2.5
#define CAMZFROM -18.0

#define CAMXTO 0.0               /* Camera direction */
#define CAMYTO 0.5
#define CAMZTO 4.0
#define CAMROLL 0.0              /* Camera roll      */

#define CAMZOOM 1.0              /* Camera zoom rate */

RtPoint CameraFrom  = { CAMXFROM, CAMYFROM, CAMZFROM }, 
        CameraTo    = { CAMXTO, CAMYTO, CAMZTO };

//#define min(a,b) ((a)<(b)?(a):(b))

void FrameCamera(float focallength, float framewidth, float frameheight)
{
    RtFloat fov;
    /* A nonzero focal length is taken to be a "normal" lens */
    if(focallength != 0.0) {
		fov = 2 * atan((min(framewidth,frameheight)*.5f)/focallength) *
	      180.0f/3.14159f;
        RiProjection("perspective", RI_FOV, (RtPointer)&fov, RI_NULL);
    } 
	else 
        RiProjection("orthographic", RI_NULL);    
    RiFrameAspectRatio((RtFloat)(framewidth/frameheight));
}

/*  AimZ(): rotate the world so the directionvector points in 
 *    positive z by rotating about the y axis, then x. The cosine 
 *    of each rotation is given by components of the normalized 
 *    direction vector. Before the y rotation the direction vector 
 *    might be in negative z, but not afterward.
*/

#define PI 3.14159265359
void AimZ(RtPoint direction)
{
    RtFloat xzlen, yzlen, yrot, xrot;

    if (direction[0]==0 && direction[1]==0 && direction[2]==0)
        return;
    /*
     * The initial rotation about the y axis is given by the projection of
     * the direction vector onto the x,z plane: the x and z components
     * of the direction. 
     */
    xzlen = sqrt(direction[0]*direction[0]+direction[2]*direction[2]);
    if (xzlen == 0)
        yrot = (direction[1] < 0) ? 180 : 0;
    else
        yrot = 180*acos(direction[2]/xzlen)/PI;

    /*
     * The second rotation, about the x axis, is given by the projection on
     * the y,z plane of the y-rotated direction vector: the original y
     * component, and the rotated x,z vector from above. 
    */
    yzlen = sqrt(direction[1]*direction[1]+xzlen*xzlen);
    xrot = 180*acos(xzlen/yzlen)/PI;       /* yzlen should never be 0 */

    if (direction[1] > 0)
        RiRotate(xrot, 1.0f, 0.0f, 0.0f);
    else
        RiRotate(-xrot, 1.0f, 0.0f, 0.0f);
    /* The last rotation declared gets performed first */
    if (direction[0] > 0)
        RiRotate(-yrot, 0.0f, 1.0f, 0.0f);
    else
        RiRotate(yrot, 0.0f, 1.0f, 0.0f);
}

/* PlaceCamera(): establish a viewpoint, viewing direction and orientation
 *    for a scene. This routine must be called before RiWorldBegin(). 
 *    position: a point giving the camera position
 *    direction: a point giving the camera direction relative to position
 *    roll: an optional rotation of the camera about its direction axis 
 */

void PlaceCamera(RtPoint position, RtPoint direction, float roll)
{
    RiIdentity();                 /* Initialize the camera transformation */
    RiRotate(-roll, 0.0, 0.0, 1.0);
    AimZ(direction);
    RiTranslate(-position[0], -position[1], -position[2]);
}

void render(char* ribName, char* imageName, MACGrid& grid, bool blobs)
{
	RiBegin(ribName);

		/* Output image characteristics */
		RiDisplay(imageName, RI_FILE, RI_RGBA, RI_NULL);
		RiShadingRate(1.0);
		RiFormat((RtInt)PICXRES, (RtInt)PICYRES, -1.0);  /* Image resolution */
		/* Region of image rendered */
		RiCropWindow(CROPMINX, CROPMAXX, CROPMINY, CROPMAXY);
		FrameCamera(PICXRES*CAMZOOM, PICXRES, PICYRES);
		/* Camera position and orientation */
		RtPoint CamTo;
		CamTo[0] = CameraTo[0] - CameraFrom[0];
		CamTo[1] = CameraTo[1] - CameraFrom[1];
		CamTo[2] = CameraTo[2] - CameraFrom[2];
		PlaceCamera(CameraFrom, CamTo, CAMROLL);
		RiClipping((float)RI_EPSILON, (float)RI_INFINITY);            /* Clipping planes*/

		/* Background Color */
		RtColor bg = {0,0,0};
		RiImager("background", "background", (RtPointer)bg, RI_NULL);

		/* Lighting */
		RtPoint from;
		//from[0] = 2.0;     
		//from[1] = 4.0;   
		//from[2] = -1.0;
		from[0] = 0.0;     
		from[1] = 0.0;   
		from[2] = 0.0;
		float intensity = 22.0;
		RiLightSource("pointlight", (RtToken)"intensity", (RtPointer)&intensity,
					  (RtToken)"from", (RtPointer)from, RI_NULL);

      /* Point Falloff */
      RtFloat falloff[1] = {0.9};
      int falloffOn = 1;
      RiAttribute((RtToken) "stochastic", (RtToken) "int pointfalloff", (RtPointer) &falloffOn, RI_NULL);
      RiHider((RtToken) "stochastic", (RtToken) "float pointfalloffpower", (RtPointer) &falloff, RI_NULL);

		/* Now describe the world */
		RiWorldBegin();
		/* Scene */	

			if (blobs) renderBlobbyBands(grid);
         else renderPoints(grid);
			
			//renderBlueCore(grid);

			// Add floor
			RiSurface("plastic", RI_NULL);
			RtFloat floorWidth = 100.0;
			RtColor floorColor = {1.0, 1.0, 1.0};
			RiColor(floorColor);
				
			RtPoint floor[4] = { 
				{floorWidth,0,floorWidth},
				{floorWidth,0,-floorWidth},
				{-floorWidth,0,-floorWidth},
				{-floorWidth,0,floorWidth} 
				};
			RiPolygon((RtInt) 4, RI_P, (RtPointer) floor, RI_NULL);

		RiWorldEnd();
        
    RiEnd();

}

void renderBlobbies(MACGrid& grid)
{
	double spacing = theCellSize;
	double gridMaxX = theDim[0]*theCellSize;
	double gridMaxY = theDim[1]*theCellSize;
	double gridMaxZ = theDim[2]*theCellSize;

	int blobbyCount = 0;
	float blobbyScale = 2*theCellSize;
	std::vector<RtInt> codeArray;
	std::vector<RtFloat> fltArray;
	std::vector<RtFloat> colorArray;
	std::vector<RtFloat> opacArray;
	for(double z = 0.0; z <= gridMaxZ; z+=spacing)
		for(double y = 0.0; y <= gridMaxY; y+=spacing)
			for(double x = 0.0; x <= gridMaxX; x+=spacing)
			{
				codeArray.push_back(1001); // OpCode for ellipsoid blobby
				codeArray.push_back(blobbyCount*16); // This blobby's index into the float array (each blobby uses 16 floats)
				fltArray.push_back(blobbyScale); fltArray.push_back(0.0); fltArray.push_back(0.0); fltArray.push_back(0.0);
				fltArray.push_back(0.0); fltArray.push_back(blobbyScale); fltArray.push_back(0.0); fltArray.push_back(0.0);
				fltArray.push_back(0.0); fltArray.push_back(0.0); fltArray.push_back(blobbyScale); fltArray.push_back(0.0);
				fltArray.push_back(x); fltArray.push_back(y); fltArray.push_back(z); fltArray.push_back(1.0);
				vec4 cv = grid.getFireRenderColor(vec3(x, y, z));
				colorArray.push_back(cv[0]);
				colorArray.push_back(cv[1]);
				colorArray.push_back(cv[2]);
				opacArray.push_back(cv[3]);
				opacArray.push_back(cv[3]);
				opacArray.push_back(cv[3]);
				blobbyCount++;
			}

	// Append 'add' OpCode (code=0), causes implicit field functions to be summed
	codeArray.push_back(0); 
	codeArray.push_back(blobbyCount);
	for(int i=0; i < blobbyCount; i++)
		codeArray.push_back(i);
	
	RiSurface("matte", RI_NULL);
	RiTranslate(-gridMaxX/2, -gridMaxY/2, (gridMaxX+gridMaxY)*0.4);
	RiBlobby(blobbyCount, codeArray.size(), &codeArray[0], fltArray.size(), &fltArray[0], 0, (RtToken*)RI_NULL, 
			(RtToken)"vertex color Cs", (RtPointer)&colorArray[0],
			(RtToken)"vertex color Os", (RtPointer)&opacArray[0], 
			RI_NULL);
}

void renderPoints(MACGrid& grid)
{
	double spacing = theCellSize*0.5;
	double gridMaxX = theDim[0]*theCellSize;
	double gridMaxY = theDim[1]*theCellSize;
	double gridMaxZ = theDim[2]*theCellSize;

	const RtFloat scale = 2.5*spacing;

	std::vector<RtFloat> firePts;
	std::vector<RtFloat> fireColorArray;
	std::vector<RtFloat> fireOpacArray;
	std::vector<RtFloat> smokePts;
	std::vector<RtFloat> smokeColorArray;
	std::vector<RtFloat> smokeOpacArray;

	for(double z = 0.0; z <= gridMaxZ; z+=spacing)
		for(double y = 0.0; y <= gridMaxY; y+=spacing)
			for(double x = 0.0; x <= gridMaxX; x+=spacing)
			{
				vec4 cv = grid.getFireRenderColor(vec3(x, y, z));
            if (cv[3] <= 0.01) continue;

				double t = grid.getTemperature(vec3(x, y, z));
            if (t < 0.5)
            {
		        smokePts.push_back(x);
				  smokePts.push_back(y);
				  smokePts.push_back(z);

				  smokeColorArray.push_back(cv[0]);
				  smokeColorArray.push_back(cv[1]);
				  smokeColorArray.push_back(cv[2]);

				  smokeOpacArray.push_back(cv[3]);
				  smokeOpacArray.push_back(cv[3]);
				  smokeOpacArray.push_back(cv[3]);
            }
            else
            {
		        firePts.push_back(x);
				  firePts.push_back(y);
				  firePts.push_back(z);

				  fireColorArray.push_back(cv[0]);
				  fireColorArray.push_back(cv[1]);
				  fireColorArray.push_back(cv[2]);

				  fireOpacArray.push_back(cv[3]);
				  fireOpacArray.push_back(cv[3]);
				  fireOpacArray.push_back(cv[3]);
            }
			}

	RiSurface("fireParticle", RI_NULL);
	RiTranslate(-gridMaxX/2, -gridMaxY/2, (gridMaxX+gridMaxY)*0.4);
	RiPoints(firePts.size()/3, (RtToken)"P", (RtPointer)&firePts[0],
				(RtToken)"constantwidth", (RtPointer)&scale, 
				(RtToken)"vertex color Cs", (RtPointer)&fireColorArray[0],
				(RtToken)"vertex color Os", (RtPointer)&fireOpacArray[0],
				RI_NULL);

	RiSurface("smokeParticle", RI_NULL);
	RiTranslate(-gridMaxX/2, -gridMaxY/2, (gridMaxX+gridMaxY)*0.4);
	RiPoints(smokePts.size()/3, (RtToken)"P", (RtPointer)&smokePts[0],
				(RtToken)"constantwidth", (RtPointer)&scale, 
				(RtToken)"vertex color Cs", (RtPointer)&smokeColorArray[0],
				(RtToken)"vertex color Os", (RtPointer)&smokeOpacArray[0],
				RI_NULL);
}

void renderBlobbyBands(MACGrid& grid)
{
	double spacing = theCellSize/4;
	double gridMaxX = theDim[0]*theCellSize;
	double gridMaxY = theDim[1]*theCellSize;
	double gridMaxZ = theDim[2]*theCellSize;

	int blobbyCount1 = 0;
	int blobbyCount2 = 0;
	int blobbyCount3 = 0;
	float blobbyScale = spacing*2.2;
	std::vector<RtInt> codeArray1;
	std::vector<RtInt> codeArray2;
	std::vector<RtInt> codeArray3;
	std::vector<RtFloat> fltArray1;
	std::vector<RtFloat> fltArray2;
	std::vector<RtFloat> fltArray3;
	std::vector<RtFloat> colorArray1;
	std::vector<RtFloat> colorArray2;
	std::vector<RtFloat> colorArray3;
	std::vector<RtFloat> opacArray1;
	std::vector<RtFloat> opacArray2;
	std::vector<RtFloat> opacArray3;
	for(double z = 0.0; z <= gridMaxZ; z+=spacing)
		for(double y = 0.0; y <= gridMaxY; y+=spacing)
			for(double x = 0.0; x <= gridMaxX; x+=spacing)
			{
				double temp = grid.getTemperature(vec3(x, y, z));
				double dens = grid.getSmokeDensity(vec3(x, y, z));
				if(temp < 0.4 && dens > 0.05)
				{
					codeArray1.push_back(1001); // OpCode for ellipsoid blobby
					codeArray1.push_back(blobbyCount1*16); // This blobby's index into the float array (each blobby uses 16 floats)
					fltArray1.push_back(blobbyScale); fltArray1.push_back(0.0); fltArray1.push_back(0.0); fltArray1.push_back(0.0);
					fltArray1.push_back(0.0); fltArray1.push_back(blobbyScale); fltArray1.push_back(0.0); fltArray1.push_back(0.0);
					fltArray1.push_back(0.0); fltArray1.push_back(0.0); fltArray1.push_back(blobbyScale); fltArray1.push_back(0.0);
					fltArray1.push_back(x); fltArray1.push_back(y); fltArray1.push_back(z); fltArray1.push_back(1.0);
					vec4 cv = grid.getFireRenderColor(vec3(x, y, z));
					colorArray1.push_back(cv[0]);
					colorArray1.push_back(cv[1]);
					colorArray1.push_back(cv[2]);
					opacArray1.push_back(cv[3]);
					opacArray1.push_back(cv[3]);
					opacArray1.push_back(cv[3]);
					blobbyCount1++;
				}
				else if((temp >= 0.4) && (temp < 0.9))
				{
					codeArray2.push_back(1001); // OpCode for ellipsoid blobby
					codeArray2.push_back(blobbyCount2*16); // This blobby's index into the float array (each blobby uses 16 floats)
					fltArray2.push_back(blobbyScale); fltArray2.push_back(0.0); fltArray2.push_back(0.0); fltArray2.push_back(0.0);
					fltArray2.push_back(0.0); fltArray2.push_back(blobbyScale); fltArray2.push_back(0.0); fltArray2.push_back(0.0);
					fltArray2.push_back(0.0); fltArray2.push_back(0.0); fltArray2.push_back(blobbyScale); fltArray2.push_back(0.0);
					fltArray2.push_back(x); fltArray2.push_back(y); fltArray2.push_back(z); fltArray2.push_back(1.0);
					vec4 cv = grid.getFireRenderColor(vec3(x, y, z));
					colorArray2.push_back(cv[0]);
					colorArray2.push_back(cv[1]);
					colorArray2.push_back(cv[2]);
					opacArray2.push_back(cv[3]);
					opacArray2.push_back(cv[3]);
					opacArray2.push_back(cv[3]);
					blobbyCount2++;
				}
				else if (temp >= 0.9)
				{
					codeArray3.push_back(1001); // OpCode for ellipsoid blobby
					codeArray3.push_back(blobbyCount3*16); // This blobby's index into the float array (each blobby uses 16 floats)
					fltArray3.push_back(blobbyScale); fltArray3.push_back(0.0); fltArray3.push_back(0.0); fltArray3.push_back(0.0);
					fltArray3.push_back(0.0); fltArray3.push_back(blobbyScale); fltArray3.push_back(0.0); fltArray3.push_back(0.0);
					fltArray3.push_back(0.0); fltArray3.push_back(0.0); fltArray3.push_back(blobbyScale); fltArray3.push_back(0.0);
					fltArray3.push_back(x); fltArray3.push_back(y); fltArray3.push_back(z); fltArray3.push_back(1.0);
					vec4 cv = grid.getFireRenderColor(vec3(x, y, z));
					colorArray3.push_back(cv[0]);
					colorArray3.push_back(cv[1]);
					colorArray3.push_back(cv[2]);
					opacArray3.push_back(cv[3]);
					opacArray3.push_back(cv[3]);
					opacArray3.push_back(cv[3]);
					blobbyCount3++;
				}
			}

	// Append 'add' OpCode (code=0), causes implicit field functions to be summed
	codeArray1.push_back(0); 
	codeArray1.push_back(blobbyCount1);
	for(int i=0; i < blobbyCount1; i++)
		codeArray1.push_back(i);

	codeArray2.push_back(0); 
	codeArray2.push_back(blobbyCount2);
	for(int i=0; i < blobbyCount2; i++)
		codeArray2.push_back(i);

	codeArray3.push_back(0); 
	codeArray3.push_back(blobbyCount3);
	for(int i=0; i < blobbyCount3; i++)
		codeArray3.push_back(i);
	
	RtColor op = {0.3, 0.3, 0.3};
	RiOpacity(op);
	RiSurface("matte", RI_NULL);
	RiTranslate(-gridMaxX/2, -gridMaxY/2, 0);
	RiBlobby(blobbyCount1, codeArray1.size(), &codeArray1[0], fltArray1.size(), &fltArray1[0], 0, (RtToken*)RI_NULL, 
			(RtToken)"vertex color Cs", (RtPointer)&colorArray1[0],
			//(RtToken)"vertex color Os", (RtPointer)&opacArray1[0], 
			RI_NULL);
	RtColor op2 = {0.6, 0.6, 0.6};
	RiOpacity(op2);
	RiSurface("particle", RI_NULL);
	RiBlobby(blobbyCount2, codeArray2.size(), &codeArray2[0], fltArray2.size(), &fltArray2[0], 0, (RtToken*)RI_NULL, 
			(RtToken)"vertex color Cs", (RtPointer)&colorArray2[0],
			//(RtToken)"vertex color Os", (RtPointer)&opacArray2[0], 
			RI_NULL);
	RiBlobby(blobbyCount3, codeArray3.size(), &codeArray3[0], fltArray3.size(), &fltArray3[0], 0, (RtToken*)RI_NULL, 
			(RtToken)"vertex color Cs", (RtPointer)&colorArray3[0],
			//(RtToken)"vertex color Os", (RtPointer)&opacArray3[0], 
			RI_NULL);
}

void renderBlueCore(MACGrid& grid)
{
	double spacing = theCellSize/6;
	double gridMaxX = theDim[0]*theCellSize;
	double gridMaxY = theDim[1]*theCellSize;
	double gridMaxZ = theDim[2]*theCellSize;

	int blobbyCount = 0;
	float blobbyScale = 2.2*spacing;
	std::vector<RtInt> codeArray;
	std::vector<RtFloat> fltArray;
	for(double z = 0.0; z <= gridMaxZ; z+=spacing)
		for(double y = 0.0; y <= gridMaxY; y+=spacing)
			for(double x = 0.0; x <= gridMaxX; x+=spacing)
			{
				if(grid.insideLevelSet(vec3(x, y, z)))
				{
					codeArray.push_back(1001); // OpCode for ellipsoid blobby
					codeArray.push_back(blobbyCount*16); // This blobby's index into the float array (each blobby uses 16 floats)
					fltArray.push_back(blobbyScale); fltArray.push_back(0.0); fltArray.push_back(0.0); fltArray.push_back(0.0);
					fltArray.push_back(0.0); fltArray.push_back(blobbyScale); fltArray.push_back(0.0); fltArray.push_back(0.0);
					fltArray.push_back(0.0); fltArray.push_back(0.0); fltArray.push_back(blobbyScale); fltArray.push_back(0.0);
					fltArray.push_back(x); fltArray.push_back(y); fltArray.push_back(z); fltArray.push_back(1.0);
					blobbyCount++;
				}
			}

	// Append 'add' OpCode (code=0), causes implicit field functions to be summed
	codeArray.push_back(0); 
	codeArray.push_back(blobbyCount);
	for(int i=0; i < blobbyCount; i++)
		codeArray.push_back(i);
	
	RiSurface("plastic", RI_NULL);
	RtColor c = {0.0, 0.0, 1.0};
	RiColor(c);
	RtColor o = {0.7, 0.7, 0.7};
	RiOpacity(o);
	RiTranslate(-gridMaxX/2, -gridMaxY/2, 0);
	RiBlobby(blobbyCount, codeArray.size(), &codeArray[0], fltArray.size(), &fltArray[0], 0, (RtToken*)RI_NULL, 
			RI_NULL);
}

#endif