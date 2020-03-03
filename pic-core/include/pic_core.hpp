#include <array>
#include <cmath>
#include <random>
#include <chrono>
#include <vector>
#include <iostream>

#define DEBUG() std::cout<<"BREAK POINT: LINE "<<__LINE__<<" IN "<<__FILE__<<std::endl


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
				float topCorner_x, float topCorner_y, float topCorner_z)
{
	particleSim.positions_x.resize(nParticles);
	particleSim.positions_y.resize(nParticles);
	particleSim.positions_z.resize(nParticles);
	particleSim.velocities_x.resize(nParticles, 0.f);
	particleSim.velocities_y.resize(nParticles, 0.f);
	particleSim.velocities_z.resize(nParticles, 0.f);
	particleSim.masses.resize(nParticles, 1);

    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_real_distribution<> x_dist(bottomCorner_x, topCorner_x);
    std::uniform_real_distribution<> y_dist(bottomCorner_y, topCorner_y);
    std::uniform_real_distribution<> z_dist(bottomCorner_z, topCorner_z);

    for (int i = 0 ; i < nParticles ; i ++)
    {
        particleSim.positions_x[i] = x_dist(eng);
        particleSim.positions_y[i] = y_dist(eng);
        particleSim.positions_z[i] = z_dist(eng);   
    }
}

struct GridAttributes 
{
	struct Constants
	{
		Constants(int dimension_x, int dimension_y, int dimension_z, float cellSize) :
		dimension_x(dimension_x + 1),
		dimension_y(dimension_y + 1),
		dimension_z(dimension_z + 1),
		cellSize(cellSize),
		attribSize(this->dimension_x * this->dimension_y * this->dimension_z)
		{
		}

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
	{
	}

	Constants const constants;

	std::vector<float> velocities_x;
	std::vector<float> velocities_y;
	std::vector<float> velocities_z;

	std::vector<float> masses;

	int flatten3DIndex(int i, int  j, int k) const
	{
		return (i * (constants.dimension_x * constants.dimension_y) + j * constants.dimension_x + k);
	}
};

void printVec3(std::array<float, 3> v)
{
	std::cout<<"x: "<<v[0]<<", y: "<<v[1]<<", z: "<<v[2]<<std::endl;
}

void printVec2(std::array<float, 2> v)
{
	std::cout<<"y: "<<v[0]<<", z: "<<v[1]<<std::endl;
}

std::array<int, 3> findGridIndex(float x, float y, float z, 
								 float cellLen_x, 
								 float cellLen_y, 
								 float cellLen_z)
{
	return {floor(x/cellLen_x), floor(y/cellLen_y), floor(z/cellLen_z)};
}

float getDiff(float const x0, float const x1, float const x)
{
	if ((x1 - x0) == 0)
		DEBUG();
	return (x - x0)/(x1 - x0);
}

float linearInterpolation(float const c0, float const c1, float const z)
{
	return c0 * (1.f - z) + c1 * z;
}

float bilinearIterpolation(std::array<float, 4> attr,
							float y, float z)
{
	auto [c00, c01, c10, c11] = attr;

	float const c0 = linearInterpolation(c00, c10, y);
	float const c1 = linearInterpolation(c01, c11, y);

	return linearInterpolation(c0, c1, z);
}

float trilinearInterpolation(std::array<float, 8> attr,
							float x, float y, float z)
{
	auto [c000, c100, c001, c101, c010, c110, c011, c111] = attr;

	float const c00 = linearInterpolation(c000, c100, x);
	float const c01 = linearInterpolation(c001, c101, x);
	float const c10 = linearInterpolation(c010, c110, x);
	float const c11 = linearInterpolation(c011, c111, x);

	return bilinearIterpolation({c00, c01, c10, c11}, y, z);
}

std::array<float, 8> getWeights(float x, float y, float z)
{
	float w000 = (1 - z) * (1 - y) *  (1 - x);
	float w100 = z * (1 - y) *  (1 - x);
	float w001 = (1 - z) *  (1 - y) * x;
	float w101 = z *  (1 - y) * x;
	float w010 = (1 - z) * y *  (1 - x);
	float w110 = z * y *  (1 - x);
	float w011 = (1 - z) *  y * x;
	float w111 = z *  y * x;
	
	return {w000, w100, w001, w101, w010, w110, w011, w111};
}

std::array<float, 3> cellSpacePos(float x, float y, float z, 
								  float cellLen_x, float cellLen_y, 
								  float cellLen_z, int cell_i,
								  int cell_j, int cell_k)
{

	float cellPos_x = cell_i * cellLen_x;
	float cellPos_y = cell_j * cellLen_y;
	float cellPos_z = cell_k * cellLen_z;

	auto cellSpace_x = getDiff(cellPos_x, cellPos_x + cellLen_x, x);
	auto cellSpace_y = getDiff(cellPos_y, cellPos_y + cellLen_y, y);
	auto cellSpace_z = getDiff(cellPos_z, cellPos_z + cellLen_z, z);

	return {cellSpace_x, cellSpace_y, cellSpace_z};
}

void transferAttributes(ParticleAttributes const & particleSim, GridAttributes & grid)
{
	//Store particle pos size to reduce verbosity and overhead
	int const nParticles = particleSim.velocities_x.size();

	for (int i = 0 ; i < nParticles ; ++i)
	{
		//Store current particle's position
		float const particle_x = particleSim.positions_x[i];
		float const particle_y = particleSim.positions_y[i];
		float const particle_z = particleSim.positions_z[i];

		float const particle_mass = particleSim.masses[i];

		float const particle_velocity_x = particleSim.velocities_x[i];
		float const particle_velocity_y = particleSim.velocities_y[i];
		float const particle_velocity_z = particleSim.velocities_z[i];

		float const cellSize = grid.constants.cellSize;

		//Get grid index from particle position
		auto [cell_i, cell_j, cell_k] = findGridIndex(particle_x, particle_y, particle_z,
													  cellSize, cellSize, cellSize);
		
		//std::cout<<"cell_i:"<<cell_i<<"cell_j:"<<cell_j<<"cell_k:"<<cell_k<<std::endl;
		
		//get corner indices
		int c000_index = grid.flatten3DIndex(cell_i, cell_j, cell_k);
		int c100_index = grid.flatten3DIndex(cell_i + 1, cell_j, cell_k);
		int c001_index = grid.flatten3DIndex(cell_i, cell_j, cell_k + 1);
		int c101_index = grid.flatten3DIndex(cell_i + 1, cell_j, cell_k + 1);
		int c010_index = grid.flatten3DIndex(cell_i, cell_j + 1, cell_k);
		int c110_index = grid.flatten3DIndex(cell_i + 1, cell_j + 1, cell_k);
		int c011_index = grid.flatten3DIndex(cell_i, cell_j + 1, cell_k + 1);
		int c111_index = grid.flatten3DIndex(cell_i + 1, cell_j + 1, cell_k + 1);

		//Get cell position from grid indices
		auto [cellSpace_x, cellSpace_y, cellSpace_z] = cellSpacePos(particle_x, particle_y, particle_z, 
					 												cellSize, cellSize, cellSize,
					 												cell_i, cell_j, cell_k);

		auto [w000, w100, w001, w101, w010, w110, w011, w111] = getWeights(cellSpace_x, cellSpace_y, cellSpace_z);

		grid.masses[c000_index] += w000 * particle_mass;
		grid.masses[c100_index] += w100 * particle_mass;
		grid.masses[c001_index] += w001 * particle_mass;
		grid.masses[c101_index] += w101 * particle_mass;
		grid.masses[c010_index] += w010 * particle_mass;
		grid.masses[c110_index] += w110 * particle_mass;
		grid.masses[c011_index] += w011 * particle_mass;
		grid.masses[c111_index] += w111 * particle_mass;

		grid.velocities_x[c000_index] += w000 * particle_velocity_x;
		grid.velocities_x[c100_index] += w100 * particle_velocity_x;
		grid.velocities_x[c001_index] += w001 * particle_velocity_x;
		grid.velocities_x[c101_index] += w101 * particle_velocity_x;
		grid.velocities_x[c010_index] += w010 * particle_velocity_x;
		grid.velocities_x[c110_index] += w110 * particle_velocity_x;
		grid.velocities_x[c011_index] += w011 * particle_velocity_x;
		grid.velocities_x[c111_index] += w111 * particle_velocity_x;

		grid.velocities_y[c000_index] += w000 * particle_velocity_y;
		grid.velocities_y[c100_index] += w100 * particle_velocity_y;
		grid.velocities_y[c001_index] += w001 * particle_velocity_y;
		grid.velocities_y[c101_index] += w101 * particle_velocity_y;
		grid.velocities_y[c010_index] += w010 * particle_velocity_y;
		grid.velocities_y[c110_index] += w110 * particle_velocity_y;
		grid.velocities_y[c011_index] += w011 * particle_velocity_y;
		grid.velocities_y[c111_index] += w111 * particle_velocity_y;

		grid.velocities_z[c000_index] += w000 * particle_velocity_z;
		grid.velocities_z[c100_index] += w100 * particle_velocity_z;
		grid.velocities_z[c001_index] += w001 * particle_velocity_z;
		grid.velocities_z[c101_index] += w101 * particle_velocity_z;
		grid.velocities_z[c010_index] += w010 * particle_velocity_z;
		grid.velocities_z[c110_index] += w110 * particle_velocity_z;
		grid.velocities_z[c011_index] += w011 * particle_velocity_z;
		grid.velocities_z[c111_index] += w111 * particle_velocity_z;
	}
}

void transferAttributes(GridAttributes const & grid, ParticleAttributes & particleSim)
{
	//Store particle pos size to reduce verbosity and overhead
	int const nParticles = particleSim.velocities_x.size();

	for (int i = 0 ; i < nParticles ; ++i)
	{
		//Store current particle's position
		float const particle_x = particleSim.positions_x[i];
		float const particle_y = particleSim.positions_y[i];
		float const particle_z = particleSim.positions_z[i];

		//float const particle_mass = particleSim.masses[i];

		float const cellSize = grid.constants.cellSize;

		//Get grid index from particle position
		auto [cell_i, cell_j, cell_k] = findGridIndex(particle_x, particle_y, particle_z,
													  cellSize, cellSize, cellSize);
		
		//std::cout<<"cell_i:"<<cell_i<<"cell_j:"<<cell_j<<"cell_k:"<<cell_k<<std::endl;
		
		//get corner indices
		int c000_index = grid.flatten3DIndex(cell_i, cell_j, cell_k);
		int c100_index = grid.flatten3DIndex(cell_i + 1, cell_j, cell_k);
		int c001_index = grid.flatten3DIndex(cell_i, cell_j, cell_k + 1);
		int c101_index = grid.flatten3DIndex(cell_i + 1, cell_j, cell_k + 1);
		int c010_index = grid.flatten3DIndex(cell_i, cell_j + 1, cell_k);
		int c110_index = grid.flatten3DIndex(cell_i + 1, cell_j + 1, cell_k);
		int c011_index = grid.flatten3DIndex(cell_i, cell_j + 1, cell_k + 1);
		int c111_index = grid.flatten3DIndex(cell_i + 1, cell_j + 1, cell_k + 1);


		//Get cell position from grid indices
		auto [cellSpace_x, cellSpace_y, cellSpace_z] = cellSpacePos(particle_x, particle_y, particle_z, 
					 												cellSize, cellSize, cellSize,
					 												cell_i, cell_j, cell_k);

		particleSim.masses[i] =
		trilinearInterpolation({grid.masses[c000_index], grid.masses[c100_index],
							    grid.masses[c001_index], grid.masses[c101_index],
							    grid.masses[c010_index], grid.masses[c110_index],
							    grid.masses[c011_index], grid.masses[c111_index]},
							    cellSpace_x, cellSpace_y, cellSpace_z);

		particleSim.velocities_x[i] =
		trilinearInterpolation({grid.velocities_x[c000_index], grid.velocities_x[c100_index],
							    grid.velocities_x[c001_index], grid.velocities_x[c101_index],
							    grid.velocities_x[c010_index], grid.velocities_x[c110_index],
							    grid.velocities_x[c011_index], grid.velocities_x[c111_index]},
							    cellSpace_x, cellSpace_y, cellSpace_z);

		particleSim.velocities_y[i] =
		trilinearInterpolation({grid.velocities_y[c000_index], grid.velocities_y[c100_index],
							    grid.velocities_y[c001_index], grid.velocities_y[c101_index],
							    grid.velocities_y[c010_index], grid.velocities_y[c110_index],
							    grid.velocities_y[c011_index], grid.velocities_y[c111_index]},
							    cellSpace_x, cellSpace_y, cellSpace_z);

		particleSim.velocities_z[i] =
		trilinearInterpolation({grid.velocities_z[c000_index], grid.velocities_z[c100_index],
							    grid.velocities_z[c001_index], grid.velocities_z[c101_index],
							    grid.velocities_z[c010_index], grid.velocities_z[c110_index],
							    grid.velocities_z[c011_index], grid.velocities_z[c111_index]},
							    cellSpace_x, cellSpace_y, cellSpace_z);

	}
}

void timeStep(ParticleAttributes & particleSim, float dt)
{
	int const nParticles = particleSim.velocities_x.size();

	for (int i = 0 ; i < nParticles ; ++i)
	{
		particleSim.velocities_y[i] += GRAV_Y;

		particleSim.positions_x[i] += particleSim.velocities_x[i] * dt;
		particleSim.positions_y[i] += particleSim.velocities_y[i] * dt;
		particleSim.positions_z[i] += particleSim.velocities_z[i] * dt;
	}
}

}

/* FUNCTION GRAVEYARD TRASHCAN

std::array<float, 3> cellPos(int cell_i, int cell_j, int cell_k,
							 float cellLen_x, 
							 float cellLen_y, 
							 float cellLen_z)
{
	return {cell_i * cellLen_x, cell_j * cellLen_y, cell_k * cellLen_z};
}

void transferMass(ParticleAttributes const & particleSim, GridAttributes* grid)
{
	//Store particle pos size to reduce verbosity and overhead
	int const nParticles = particleSim.velocities_x.size();

	for (int i = 0 ; i < nParticles ; ++i)
	{
		//Store current particle's position
		float const particle_x = particleSim.positions_x[i];
		float const particle_y = particleSim.positions_y[i];
		float const particle_z = particleSim.positions_z[i];

		//Get grid index from particle position
		auto [cell_i, cell_j, cell_k] = findGridIndex(particle_x, particle_y, particle_z, grid->constants.cellSize);


		//Get cell position from grid indices
		auto [cell_x, cell_y, cell_z] = cellPos(cell_i, cell_j, cell_k, grid->constants.cellSize);

		//Calculate weight using the interpolation kernel 
		//float const weight = trilinearCellInterpolation(particle_x - cell_x,
		//												particle_y - cell_y,
		//												particle_z - cell_z,
		//												cell_x,
		//												cell_y,
		//												cell_z,
		//												grid->constants.cellSize);

		float const weight = trilinearInterpolation({}, 
												   particle_x - cell_x,
												   particle_y - cell_y,
												   particle_z - cell_z);

		//Store the interpolated mass.
		grid->masses[grid->flatten3DIndex(cell_i, cell_j, cell_k)] += weight * particleSim.masses[i];  
	}
}

std::array<float, 8> getWeights(float x, float y, float z,
								float cellSize)
{
	float w0 = (1 - x / cellSize);
	float w1 = (x / cellSize);

	float w00 = (1 - y / cellSize) * w0;
	float w10 = (y / cellSize) * w0;
	float w01 = (1 - y / cellSize) * w1;
	float w11 = (y / cellSize) * w1;

	float w000 = (1 - z / cellSize) * w00;
	float w100 = (z / cellSize) * w00;
	float w001 = (1 - z / cellSize) * w01;
	float w101 = (z / cellSize) * w01;
	float w010 = (1 - z / cellSize) * w10;
	float w110 = (z / cellSize) * w10;
	float w011 = (1 - z / cellSize) * w11;
	float w111 = (z / cellSize) * w11;

	return {w000, w100, w001, w101, w010, w110, w011, w111};
}

float trilinearCellInterpolation(float x, float y, float z,
							     float cellPos_x, float cellPos_y,
							     float cellPos_z, float cellSize)
{
	// Calculate the 8 cube corner positions
	float const halfCellSize = cellSize * 0.5f;

	std::array<float, 3> const c000{cellPos_x - halfCellSize, cellPos_y - halfCellSize, cellPos_z - halfCellSize};
	std::array<float, 3> const c001{cellPos_x - halfCellSize, cellPos_y - halfCellSize, cellPos_z + halfCellSize};
	std::array<float, 3> const c010{cellPos_x - halfCellSize, cellPos_y + halfCellSize, cellPos_z - halfCellSize};
	std::array<float, 3> const c011{cellPos_x - halfCellSize, cellPos_y + halfCellSize, cellPos_z + halfCellSize};
	std::array<float, 3> const c100{cellPos_x + halfCellSize, cellPos_y - halfCellSize, cellPos_z - halfCellSize};
	std::array<float, 3> const c101{cellPos_x + halfCellSize, cellPos_y - halfCellSize, cellPos_z + halfCellSize};
	std::array<float, 3> const c110{cellPos_x + halfCellSize, cellPos_y + halfCellSize, cellPos_z - halfCellSize};
	std::array<float, 3> const c111{cellPos_x + halfCellSize, cellPos_y + halfCellSize, cellPos_z + halfCellSize};

	// Util to scale a vec3
	auto const mul3 = [](auto const& v, float m) -> std::array<float, 3> { return {v[0]*m, v[1]*m, v[2]*m}; };

	// Util to add two vec3
	auto const add3 = [](auto const& a, auto const& b) -> std::array<float, 3> { return {a[0]+b[0], a[1]+b[1], a[2]+b[2]}; };

	// Drop the x component from a vec3
	auto const drop_x = [](auto const& v) -> std::array<float, 2> { return {v[1], v[2]}; };

	// Interpolate the x dimension to produce the four corners of our 2D square
	std::array<float, 2> const c00 = drop_x(add3(mul3(c000, 1.f - x), mul3(c100, x)));
	std::array<float, 2> const c01 = drop_x(add3(mul3(c001, 1.f - x), mul3(c101, x)));
	std::array<float, 2> const c10 = drop_x(add3(mul3(c010, 1.f - x), mul3(c110, x)));
	std::array<float, 2> const c11 = drop_x(add3(mul3(c011, 1.f - x), mul3(c111, x)));

	printVec2(c00);
	printVec2(c01);
	printVec2(c10);
	printVec2(c11);

	// Util to scale a vec3
	auto const mul2 = [](auto const& v, float m) -> std::array<float, 2> { return {v[0]*m, v[1]*m}; };

	// Util to add two vec3
	auto const add2 = [](auto const& a, auto const& b) -> std::array<float, 2> { return {a[0]+b[0], a[1]+b[1]}; };

	// Drop the x component from a vec3
	auto const drop_y = [](auto const& v) -> float { return v[2]; };

	// Pass the square to bilinear interpolation
	float const c0 = drop_y(add2(mul2(c00, 1.f - y), mul2(c10, y)));
	float const c1 = drop_y(add2(mul2(c01, 1.f - y), mul2(c11, y)));

	// Finally pass these on to linear interpolation
	return linearInterpolation(c0, c1, z);
}

*/