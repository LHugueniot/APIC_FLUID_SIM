#include "mathUtils.h"


namespace pic
{

#define GRAV_Y -9.8

struct ParticleAttributes
{
	std::vector<float> positions_x;
	std::vector<float> positions_y;
	std::vector<float> positions_z;

	std::vector<float> velocities_x;
	std::vector<float> velocities_y;
	std::vector<float> velocities_z;

	std::vector<float> masses;
};

void randomisedParticleBB(ParticleAttributes & particleSim,  int nParticles,
	float bottomCorner_x, float bottomCorner_y, float bottomCorner_z,
	float topCorner_x, float topCorner_y, float topCorner_z);

struct GridAttributes 
{
	struct Constants
	{
		Constants(int dimension_x, int dimension_y, int dimension_z, float cellSize):
			dimension_x(dimension_x + 1),
			dimension_y(dimension_y + 1),
			dimension_z(dimension_z + 1),
			cellSize(cellSize),
			attribSize(this->dimension_x 
					 * this->dimension_y 
					 * this->dimension_z)
			{}

		int dimension_x;
		int dimension_y;
		int dimension_z;

		float cellSize;

		int attribSize;
	};

	GridAttributes(int dimension_x, int dimension_y, int dimension_z, float cellSize):
		constants(dimension_x, dimension_y, dimension_z, cellSize),
		velocities_x(constants.attribSize, 0),
		velocities_y(constants.attribSize, 0),
		velocities_z(constants.attribSize, 0),
		masses(constants.attribSize, 0)
	{}

	Constants const constants;

	std::vector<float> velocities_x;
	std::vector<float> velocities_y;
	std::vector<float> velocities_z;

	std::vector<float> masses;

	int flatten3DIndex(int i, int  j, int k) const
	{
		return (i * (constants.dimension_x * constants.dimension_y)
			+ j * constants.dimension_x + k);
	}
};


std::array<int, 3> findGridIndex(float x, float y, float z, 
	float cellLen_x, float cellLen_y, float cellLen_z);

std::array<float, 3> cellSpacePos(float x, float y, float z, float cellLen_x, 
	float cellLen_y, float cellLen_z, int cell_i,int cell_j, int cell_k);

void transferAttributes(ParticleAttributes const & particleSim, GridAttributes & grid);
void transferAttributes(GridAttributes const & grid, ParticleAttributes & particleSim);

template<typename CollisionProc>
void timeStep(ParticleAttributes & particleSim,  CollisionProc cp, float dt)
{
	int const nParticles = particleSim.velocities_x.size();

	for (int i = 0 ; i < nParticles ; ++i)
	{

		std::array<double,3> particlePos = cp(particleSim.positions_x[i],
				particleSim.positions_y[i], particleSim.positions_z[i],
				particleSim.velocities_x[i] * dt, 
				(particleSim.velocities_y[i] + GRAV_Y) * dt,
				particleSim.velocities_z[i] * dt);

		particleSim.positions_x[i] = particlePos[0];
		particleSim.positions_y[i] = particlePos[1];
		particleSim.positions_z[i] = particlePos[2];
	}
}

}