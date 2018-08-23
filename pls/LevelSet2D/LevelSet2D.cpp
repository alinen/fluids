#include "LevelSet2D.h"
#include "ParticleSet2D.h"
#include "Particle2D.h"

void LevelSet2D::Update(const Velocity2D& grid, const Double &dt)
{
	//First Order time integration
	FOR_LS2D
		SemiLagrangianStep(i,j,grid,dt);
	END_FOR_TWO

	swap(gridPhi, gridTmp);
	gridPhi.SetBoundarySignedDist();
}

void LevelSet2D::SemiLagrangianStep(int x, int y, const Velocity2D &grid, const Double &dt)
{
    static int r,s;
    static Double a,b;
	if(gridPhi(x,y) > SEMILAGRA_LIMIT) {
		gridTmp(x,y) = gridPhi(x,y);
		return;
	}
	static Vec2D u; //obtain from velocity grid
	grid.GetVelocity(Vec2D(x,y), u);

    r = x - int(ceil(u[0] * dt * hInv));
	s = y - int(ceil(u[1] * dt * hInv));

    //doesn't get from boundary
    r = Clamp(r,0,Nx);
    s = Clamp(s,0,Ny);

	a = (Double(x - r) * h - u[0] * dt) * hInv;
	b = (Double(y - s) * h - u[1] * dt) * hInv;

	gridTmp(x,y) =    a  *    b  * gridPhi(r+1, s+1) +
				   (1-a) *    b  * gridPhi(r  , s+1) +
				      a  * (1-b) * gridPhi(r+1, s  ) +
				   (1-a) * (1-b) * gridPhi(r  , s  );
}

void LevelSet2D::ReInitialize() {
    FOR_GRID2D gridFM.Set(i, gridPhi[i]);
    gridFM.Reinitialize();
    FOR_GRID2D gridPhi[i] = gridFM[i];
    gridPhi.SetBoundarySignedDist();
}
    
void LevelSet2D::Fix(const ParticleSet2D& particleSet)
{	
	static ParticleSet2D::cIterator pit, end;
	static Vec2D pos;
    static Double phi;
    static int sign;

    end = particleSet.end();
	gridPos = gridPhi;
	gridNeg = gridPhi;
	for(pit = particleSet.begin(); pit != end; ++pit)
	{
		Particle2D& particle = *(*pit);
		particle.GetPosition(pos);
		phi = LinearSample(pos);
		sign = particle.Sign();

		if(phi * sign < 0.)
		{
			//particle has crossed the boundary
			if(sign < 0.) FixNeg(particle, int(pos[0]), int(pos[1]));
			else		  FixPos(particle, int(pos[0]), int(pos[1]));
		}
	}
	//Merge gridPos & gridNeg
	static Double phiPos, phiNeg;
    //cout << endl;
	FOR_ALL_LS2D
		phiPos = gridPos(i,j); phiNeg = gridNeg(i,j);
        //if(phiPos * phiNeg < 0) cout << i << " " << j << "  +: " << phiPos << "  -: " << phiNeg << endl; 
		gridPhi(i,j) = abs(phiPos) < abs(phiNeg) ? phiPos : phiNeg;
        //gridPhi(i,j) = phiNeg;
	END_FOR_TWO
}

inline void LevelSet2D::FixNeg(const Particle2D &particle, int i, int j)
{
	static Double particlePhi;
	for(int dx = 0; dx < 2; dx++) {
		for(int dy = 0; dy < 2; dy++) {
			particlePhi = particle.phi(Vec2D(i+dx,j+dy), h);
            gridNeg(i+dx,j+dy) = min(particlePhi, gridNeg(i+dx,j+dy));
		}
	}
}

inline void LevelSet2D::FixPos(const Particle2D &particle, int i, int j)
{
	static Double particlePhi;
	for(int dx = 0; dx < 2; dx++) {
		for(int dy = 0; dy < 2; dy++) {
			particlePhi = particle.phi(Vec2D(i+dx,j+dy), h);
            gridPos(i+dx,j+dy) = max(particlePhi, gridPos(i+dx,j+dy));
		}
    }
}

inline Double LevelSet2D::LinearSample(const Vec2D &p) const
{
	static Double xlerp,ylerp;
    xlerp = p[0]-int(p[0]);
	ylerp = p[1]-int(p[1]);

    return Lerp(ylerp,
                Lerp( xlerp, gridPhi( int(p[0])  , int(p[1])   ), 
                             gridPhi( int(p[0])+1, int(p[1])   ) ),
			    Lerp( xlerp, gridPhi( int(p[0])  , int(p[1])+1 ), 
                             gridPhi( int(p[0])+1, int(p[1])+1 ) ) );	
}

Double LevelSet2D::CubicSample(const Vec2D &pos) const
{
    static Double s, t;
    static int i0,i1,i2,i3,j0,j1,j2,j3;
    
    i0 = int(pos[0]) - 1; i1 = i0+1; i2 = i1+1; i3 = i2+1;
    j0 = int(pos[1]) - 1; j1 = j0+1; j2 = j1+1; j3 = j2+1;
    s = pos[0] - i1;
    t = pos[1] - j1;
    i0 = i0 < 0 ? 0 : i0; i3 = i3 > Nx+1 ? Nx+1 : i3;
    j0 = j0 < 0 ? 0 : j0; j3 = j3 > Ny+1 ? Ny+1 : j3;

	return  MCerp(s, MCerp(t, gridPhi(i0,j0),gridPhi(i0,j1),gridPhi(i0,j2),gridPhi(i0,j3)),
                     MCerp(t, gridPhi(i1,j0),gridPhi(i1,j1),gridPhi(i1,j2),gridPhi(i1,j3)),
                     MCerp(t, gridPhi(i2,j0),gridPhi(i2,j1),gridPhi(i2,j2),gridPhi(i2,j3)),
                     MCerp(t, gridPhi(i3,j0),gridPhi(i3,j1),gridPhi(i3,j2),gridPhi(i3,j3)) );
}

void LevelSet2D::normal(const Vec2D &pos, Vec2D &n) const
{
	gradient(pos, n);
	n.Normalize();
}

void LevelSet2D::gradient(const Vec2D &pos, Vec2D &g) const
{
	static Vec2D dx = (1. , 0.);
	static Vec2D dy = (0. , 1.);
	g[0] = LinearSample(pos+dx) - LinearSample(pos-dx);
	g[1] = LinearSample(pos+dy) - LinearSample(pos-dy);
	g *= 0.5 * hInv;
}

void LevelSet2D::gradient(const Vec2D &pos, const Vec2D &u, Vec2D &g)
{
	static Vec2D dx = (1., 0.);
	static Vec2D dy = (0., 1.);
	Double center = LinearSample(pos);

	if(u[0] < 0)	g[0] = (LinearSample(pos + dx) - center) * hInv;
	else			g[0] = (center - LinearSample(pos - dx)) * hInv;

	if(u[1] < 0)	g[1] = (LinearSample(pos + dy) - center) * hInv;
	else			g[1] = (center - LinearSample(pos - dy)) * hInv;
}