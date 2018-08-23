#ifndef renderman_h
#ifdef HAVE_PRMAN

#include "MACGrid.h"
void render(char* ribName, char* imageName, MACGrid& grid, bool blob = true);

void renderBlobbies(MACGrid& grid);
void renderPoints(MACGrid& grid);
void renderBlobbyBands(MACGrid& grid);
void renderBlueCore(MACGrid& grid);

#endif
#endif