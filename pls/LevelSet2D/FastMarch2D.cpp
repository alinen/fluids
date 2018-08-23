#include "FastMarch2D.h"

FastMarch2D::FastMarch2D(int nx,int ny, Double hi) 
    : Nx(nx), Ny(ny), h(hi), hInv(1./hi), size((nx+2)*(ny+2)), NxInv(1./(nx+2))
{
	FMHeap.resize(size);
	ClosePoints.resize(size);
    ClosePoints[1] = 10;
	grid = new FMContainer[size];
}
void FastMarch2D::PrintFlags()
{
    cout << endl;
    for(int j=0; j < Ny+2; j++) {
        for(int i=0; i < Nx+2; i++) {
            if(grid[GI(i,j)].DoneFlag >= 0) cout << " ";
            cout << grid[GI(i,j)].DoneFlag;
        }
        cout << endl;
    }
}

void FastMarch2D::PrintValues()
{
    cout << endl;
    for(int j=0; j < Ny+2; j++){
        for(int i=0; i < Nx+2; i++) 
            cout << grid[GI(i,j)].value << " ";
        cout << endl;
    }
}

void FastMarch2D::Reinitialize() {
    //Negative Phi first
	ReinitHalf();
    //Then Positive Phi
    FOR_GRID2D Set(i, grid[i].value);
    ReinitHalf();
}

inline void FastMarch2D::ReinitHalf() {
    heapSize = 0;
	closeSize = 0;
    SetBoundary();
	Initialize();
    //PrintFlags();
	InitHeap();
	FastMarch();
}

inline void FastMarch2D::Set(int index, const Double &value) {
	grid[index].value = -value; // - is to swap values for initializing negative phi
	grid[index].HeapPosition = -1;
	if(grid[index].value < 0.) grid[index].DoneFlag = -1;
	else                       grid[index].DoneFlag = 0;
}

void FastMarch2D::SetBoundary()
{
	for(int i=1 ; i<=Ny; i++ ) grid[GI(0   ,i)].DoneFlag = grid[GI(Nx+1,i)].DoneFlag = -1;
	for(int i=1 ; i<=Nx; i++ ) grid[GI(i,0   )].DoneFlag = grid[GI(i,Ny+1)].DoneFlag = -1;
	grid[GI(0   ,0   )].DoneFlag = grid[GI(0   ,Ny+1)].DoneFlag = -1;
	grid[GI(Nx+1,0   )].DoneFlag = grid[GI(Nx+1,Ny+1)].DoneFlag = -1;
}

void FastMarch2D::Initialize() {
    static int ci, pi, flag;
	for(int i = 1; i <= Nx; i++)
	{
        //if(Done == 0)  f= 0: change --> -1^f = -1;           NO -->  0^f = 0; 1^f = 1;
        //if(Done == -1) f=-1: change -->  0^f = -1; 1^f = -2; NO --> -1^f = 0;
        //if(Done == 1)  f= 1: change --> -1^f = -2;           NO -->  0^f = 1; 1^f = 0;
        flag = grid[GI(i,1)].DoneFlag;
		for(int j = 2; j <= Ny; j++)
		{
            ci = GI(i,j);
            if((flag ^ grid[ci].DoneFlag) < 0)
            {
                pi = GI(i,j-1);
                flag = grid[ci].DoneFlag;
                if(flag >= 0) {
                    grid[ci].DoneFlag = 1;
                    AddClose(GI(i,j+1));
                    grid[ci].value = min(grid[ci].value, abs(h + grid[pi].value));
                }
                else {
                    grid[pi].DoneFlag = 1;
                    AddClose(GI(i,j-2));
                    grid[pi].value = min(grid[pi].value, abs(h + grid[ci].value));
                }
            }
        }
    }
    
	for(int j = 1; j <= Ny; j++)
	{
        flag = grid[GI(1,j)].DoneFlag;
		for(int i = 2; i <= Nx; i++)
		{
            ci = GI(i,j);
            if((flag ^ grid[ci].DoneFlag) < 0)
            {
                pi = GI(i-1,j);
                flag = grid[ci].DoneFlag;
                if(flag >= 0) {
                    grid[ci].DoneFlag = 1;
                    AddClose(GI(i+1,j));
                    grid[ci].value = min(grid[ci].value, abs(h + grid[pi].value));
                }
                else {
                    grid[pi].DoneFlag = 1;
                    AddClose(GI(i-2,j));
                    grid[pi].value = min(grid[pi].value, abs(h + grid[ci].value));
                }
            }
		}
	}
}

inline void FastMarch2D::AddClose(int index)
{
	if(grid[index].DoneFlag == 0) {
		//Add index to list for initialization of close band
		ClosePoints[closeSize] = index;
		closeSize++;
	}
}

void FastMarch2D::InitHeap() {
    static int x,y;
	for(int i = 0; i < closeSize; i++) {
		if(grid[ClosePoints[i]].HeapPosition == -1 && grid[ClosePoints[i]].DoneFlag == 0) {
			GIJ(ClosePoints[i], x, y);
			FindPhi(ClosePoints[i], x, y);
		}
	}
}

void FastMarch2D::FindPhi(int index, int x, int y) {
	static Double phiX, phiY, phiZ, b, quotient, phi;
    static int a;
    static bool flagX, flagY;

    phiX = phiY = phiZ = 0.;
    a = 0;
    flagX = flagY = 0;

	//Find The phiS
	CheckFront (phiX, a, flagX, GI(x+1,  y));
	CheckBehind(phiX, a, flagX, GI(x-1,  y));
	CheckFront (phiY, a, flagY, GI(x  ,y+1));
	CheckBehind(phiY, a, flagY, GI(x  ,y-1));

	//Max Tests
	if(a == 2) {
		if(phiX >= phiY) CheckMax2(a, phiX, phiY);
		else			 CheckMax2(a, phiY, phiX);
	}

	b = phiX + phiY + phiZ;
	quotient = square(b) - 
		Double(a) * (square(phiX) + square(phiY) + square(phiZ) - square(h));
	if(quotient < 0.) cout << "0 ";
	else {
		phi = b + sqrt(quotient);
		phi /= Double(a);
		grid[index].value = phi;
		if(grid[index].HeapPosition == -1) AddToHeap(index);
	    else                               UpdateHeap(index); 
	}
}

inline void FastMarch2D::CheckFront(Double& phi, int& a, bool& flag, int index) 
{
	if(grid[index].DoneFlag == 1) {
		phi = grid[index].value;
		flag = 1;
		a++;
	}
}

inline void FastMarch2D::CheckBehind(Double& phi, int& a, bool& flag, int index)
{
	if(grid[index].DoneFlag == 1) {
		if(!flag) { phi = grid[index].value; a++; }
		else phi = min(grid[index].value, phi);
        flag = 1;
	}
}

inline void FastMarch2D::CheckMax2(int& a, Double& phi1, const Double &phi2) {
	if(square((phi1 - phi2)*hInv) > 1.) { phi1 = 0; a = 1; }
}

void FastMarch2D::AddToHeap(int index) {
	FMHeap[heapSize] = index;
	grid[index].HeapPosition = heapSize;
	int j, i = heapSize;
	for(i; i > 0; i = j) {
		j = (i-1)/2;
		if (grid[ (FMHeap[i]) ].value < grid[ (FMHeap[j]) ].value ){
			FMHeap[i] = FMHeap[j];
			grid[ (FMHeap[j]) ].HeapPosition = i;
			FMHeap[j] = index;
			grid[index].HeapPosition = j;
		}
		else
			break;
	}
	heapSize++;
}

void FastMarch2D::UpdateHeap(int index) {
	int j, i = grid[index].HeapPosition;
	for(i; i > 0; i = j) {
		j = (i-1)/2;
		if (grid[ (FMHeap[i]) ].value < grid[ (FMHeap[j]) ].value ) {
			FMHeap[i] = FMHeap[j];
			grid[ (FMHeap[j]) ].HeapPosition = i;
			FMHeap[j] = index;
			grid[index].HeapPosition = j;
		}
		else
			break;
	}
}

void FastMarch2D::FastMarch() {
	static int x, y; 
	for(int index = PopHeap(); index != -1; index = PopHeap()) {
		if(grid[index].value > FASTMARCH_LIMIT) return;
		GIJ(index, x, y);
        if(grid[GI(x-1,y)].DoneFlag == 0) FindPhi(GI(x-1,y),x-1,y);
        if(grid[GI(x+1,y)].DoneFlag == 0) FindPhi(GI(x+1,y),x+1,y);
        if(grid[GI(x,y-1)].DoneFlag == 0) FindPhi(GI(x,y-1),x,y-1);
        if(grid[GI(x,y+1)].DoneFlag == 0) FindPhi(GI(x,y+1),x,y+1);
	}
}

int FastMarch2D::PopHeap() {
	if(heapSize == 0)
		return -1;
	int j, index = FMHeap[0];
	grid[index].DoneFlag = 1;
	heapSize--;
	FMHeap[0] = FMHeap[heapSize];
	grid[FMHeap[heapSize]].HeapPosition = 0;
	for(int i = 0; i < (heapSize-1); i = j) {
		int lc = 2*i+1;
		int rc = 2*i+2;
		Double current = grid[ (FMHeap[i]) ].value;
		Double lv, rv;
		if(lc < heapSize) {
			lv = grid[ (FMHeap[lc]) ].value;
			if(rc < heapSize) {
				rv = grid[ (FMHeap[rc]) ].value;
				if(lv > rv) {
					lc = rc;
					lv = rv;
				}
			}
			if(current > lv) {
				FMHeap[i] = FMHeap[lc];
				grid[ FMHeap[i] ].HeapPosition = i;
				FMHeap[lc] = FMHeap[heapSize];
				grid[ FMHeap[heapSize] ].HeapPosition = lc;
				j = lc;
			}
			else
				break;
		}
		else
			break;
	}
	return index;
}