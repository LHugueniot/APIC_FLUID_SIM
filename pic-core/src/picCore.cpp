#include "picCore.h"

namespace pic
{

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

std::array<int, 3> findGridIndex(float x, float y, float z, float cellLen_x,
	float cellLen_y, float cellLen_z)
{
	return {floor(x/cellLen_x), floor(y/cellLen_y), floor(z/cellLen_z)};
}


float3 cellSpacePos(float x, float y, float z, float cellLen_x, 
	float cellLen_y, float cellLen_z, int cell_i, int cell_j, int cell_k)
{

	float cellPos_x = cell_i * cellLen_x;
	float cellPos_y = cell_j * cellLen_y;
	float cellPos_z = cell_k * cellLen_z;

	auto cellSpace_x = getDiff(cellPos_x, cellPos_x + cellLen_x, x);
		DEBUG();
	auto cellSpace_y = getDiff(cellPos_y, cellPos_y + cellLen_y, y);
		DEBUG();
	auto cellSpace_z = getDiff(cellPos_z, cellPos_z + cellLen_z, z);
		DEBUG();

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
		
		std::cout<<"cell_i: "<<cell_i<<" cell_j: "<<cell_j<<" cell_k: "<<cell_k<<std::endl;
		
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

		std::cout<<"cell_i: "<<cell_i<<" cell_j: "<<cell_j<<" cell_k: "<<cell_k<<std::endl;

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

//template<typename CollisionProc>
//void timeStep(ParticleAttributes & particleSim,  CollisionProc cp, float dt)
//{
//	int const nParticles = particleSim.velocities_x.size();
//
//	for (int i = 0 ; i < nParticles ; ++i)
//	{
//
//		std::array<float,3> particlePos = cp(particleSim.positions_x[i],
//				particleSim.positions_y[i], particleSim.positions_z[i],
//				particleSim.velocities_x[i] * dt, 
//				(particleSim.velocities_y[i] + GRAV_Y) * dt,
//				particleSim.velocities_z[i] * dt);
//
//		particleSim.positions_x[i] = particlePos[0];
//		particleSim.positions_y[i] = particlePos[1];
//		particleSim.positions_z[i] = particlePos[2];
//	}
//}

}

/* FUNCTION GRAVEYARD TRASHCAN

float3 cellPos(int cell_i, int cell_j, int cell_k,
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

	float3 const c000{cellPos_x - halfCellSize, cellPos_y - halfCellSize, cellPos_z - halfCellSize};
	float3 const c001{cellPos_x - halfCellSize, cellPos_y - halfCellSize, cellPos_z + halfCellSize};
	float3 const c010{cellPos_x - halfCellSize, cellPos_y + halfCellSize, cellPos_z - halfCellSize};
	float3 const c011{cellPos_x - halfCellSize, cellPos_y + halfCellSize, cellPos_z + halfCellSize};
	float3 const c100{cellPos_x + halfCellSize, cellPos_y - halfCellSize, cellPos_z - halfCellSize};
	float3 const c101{cellPos_x + halfCellSize, cellPos_y - halfCellSize, cellPos_z + halfCellSize};
	float3 const c110{cellPos_x + halfCellSize, cellPos_y + halfCellSize, cellPos_z - halfCellSize};
	float3 const c111{cellPos_x + halfCellSize, cellPos_y + halfCellSize, cellPos_z + halfCellSize};

	// Util to scale a vec3
	auto const mul3 = [](auto const& v, float m) -> float3 { return {v[0]*m, v[1]*m, v[2]*m}; };

	// Util to add two vec3
	auto const add3 = [](auto const& a, auto const& b) -> float3 { return {a[0]+b[0], a[1]+b[1], a[2]+b[2]}; };

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

auto const dotProd = 
[](auto const& v1, auto const& v2)
-> double 
{ return {v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]}; };

auto const scale = 
[](auto const& v1, double s)
-> std::array<float,3>
{ return {v1[0] * s + v1[1] * s + v1[2] * s}; };

auto const project =
[](auto const& v1, auto const& v2)
-> std::array<float,3> 
{ return scale(v1, dotProd(v1, v2) / dotProd(v1, v1)) };

*/