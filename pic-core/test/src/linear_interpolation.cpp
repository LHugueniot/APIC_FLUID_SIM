#include <catch2.hpp>
#include <cstdlib>

#include "picCore.h"

using namespace pic;

TEST_CASE("Linear interpolation", "[linear][interpolation]") 
{
    REQUIRE(linearInterpolation(0.f, 1.f, 0.5f) == 0.5f);
}


TEST_CASE("Trilinear interpolation", "[trilinear][interpolation]") 
{
    REQUIRE(trilinearInterpolation({1, 1, 1, 1, 1, 1, 1, 1},
    								0.5f, 0.5f, 0.5f) == 1.0f);

    REQUIRE(trilinearInterpolation({0.5f, 0.5f, 0.5f, 0.5f, 
    								0.5f, 0.5f, 0.5f, 0.5f},
    								0.5f, 0.5f, 0.5f) == 0.5f);
}

TEST_CASE("Difference", "[get][diff]")
{
	REQUIRE(getDiff(0, 10, 3) == 0.3f);
	REQUIRE(getDiff(5, 25, 10) == 0.25f);
}

TEST_CASE("Transfer attribute", "[transfer][attribute]")
{
	GridAttributes grid(10, 10, 10, 1);

	float const particle_x = 2.5f;
	float const particle_y = 2.5f;
	float const particle_z = 2.5f;

	float const particle_mass = 1.f;

	float const cellSize = grid.constants.cellSize;

	//Get grid index from particle position
	auto [cell_i, cell_j, cell_k] = findGridIndex(particle_x, particle_y, particle_z,
												  cellSize, cellSize, cellSize);
	//std::cout<<"cell_i:"<<cell_i<<"cell_j:"<<cell_j<<"cell_k:"<<cell_k<<std::endl;
	//DEBUG();

	int c000_index = grid.flatten3DIndex(cell_i, cell_j, cell_k);
	int c100_index = grid.flatten3DIndex(cell_i + 1, cell_j, cell_k);
	int c001_index = grid.flatten3DIndex(cell_i, cell_j, cell_k + 1);
	int c101_index = grid.flatten3DIndex(cell_i + 1, cell_j, cell_k + 1);
	int c010_index = grid.flatten3DIndex(cell_i, cell_j + 1, cell_k);
	int c110_index = grid.flatten3DIndex(cell_i + 1, cell_j + 1, cell_k);
	int c011_index = grid.flatten3DIndex(cell_i, cell_j + 1, cell_k + 1);
	int c111_index = grid.flatten3DIndex(cell_i + 1, cell_j + 1, cell_k + 1);

	//DEBUG();

	//Get cell space position
	auto [cellSpace_x, cellSpace_y, cellSpace_z] = cellSpacePos(particle_x, particle_y, particle_z, 
				 												cellSize, cellSize, cellSize,
				 												cell_i, cell_j, cell_k);

	//Get interpolated weights
	auto [w000, w100, w001, w101, w010, w110, w011, w111] = getWeights(cellSpace_x, cellSpace_y, cellSpace_z);

	//Store the interpolated mass.
	//std::cout<<"c000_index "<<c000_index<<std::endl;
	//std::cout<<"c100_index "<<c100_index<<std::endl;
	//std::cout<<"c001_index "<<c001_index<<std::endl;
	//std::cout<<"c101_index "<<c101_index<<std::endl;
	//std::cout<<"c010_index "<<c010_index<<std::endl;
	//std::cout<<"c110_index "<<c110_index<<std::endl;
	//std::cout<<"c011_index "<<c011_index<<std::endl;
	//std::cout<<"c111_index "<<c111_index<<std::endl;
	//DEBUG();
	grid.masses[c000_index] += w000 * particle_mass;
	grid.masses[c100_index] += w100 * particle_mass;
	grid.masses[c001_index] += w001 * particle_mass;
	grid.masses[c101_index] += w101 * particle_mass;
	grid.masses[c010_index] += w010 * particle_mass;
	grid.masses[c110_index] += w110 * particle_mass;
	grid.masses[c011_index] += w011 * particle_mass;
	grid.masses[c111_index] += w111 * particle_mass;
	//DEBUG();

	//std::cout<<grid.masses[c000_index]<<std::endl;
	//std::cout<<grid.masses[c100_index]<<std::endl;
	//std::cout<<grid.masses[c001_index]<<std::endl;
	//std::cout<<grid.masses[c101_index]<<std::endl;
	//std::cout<<grid.masses[c010_index]<<std::endl;
	//std::cout<<grid.masses[c110_index]<<std::endl;
	//std::cout<<grid.masses[c011_index]<<std::endl;
	//std::cout<<grid.masses[c111_index]<<std::endl;
}

TEST_CASE("Randomised particle bounding box", "[Randomised][Particle][Bounding][Box]")
{
	int nParticles = 400;
	ParticleAttributes particleSim;
	randomisedParticleBB(particleSim, nParticles, 
						 0, 0, 0,
						 10, 10, 10);

	REQUIRE(particleSim.positions_x.size() == nParticles);
	REQUIRE(particleSim.positions_y.size() == nParticles);
	REQUIRE(particleSim.positions_z.size() == nParticles);
	REQUIRE(particleSim.velocities_x.size() == nParticles);
	REQUIRE(particleSim.velocities_y.size() == nParticles);
	REQUIRE(particleSim.velocities_z.size() == nParticles);
	REQUIRE(particleSim.masses.size() == nParticles);

	for (int i = 0 ; i < nParticles ; i ++)
	{
		REQUIRE(particleSim.positions_x[i] > 0);
		REQUIRE(particleSim.positions_x[i] < 10);
		REQUIRE(particleSim.positions_y[i] > 0);
		REQUIRE(particleSim.positions_y[i] < 10);
		REQUIRE(particleSim.positions_z[i] > 0);
		REQUIRE(particleSim.positions_z[i] < 10);
	}
}

TEST_CASE("Transfer Attributes particles 2 grid", "[Transfer][Attributes]")
{
	int nParticles = 5000;
	ParticleAttributes particleSim;
	randomisedParticleBB(particleSim, nParticles, 
						 0, 0, 0,
						 10, 10, 10);

	GridAttributes grid(10, 10, 10, 1);

	for (int i = 0 ; i < nParticles ; i ++)
	{
		REQUIRE(particleSim.positions_x[i] > 0);
		REQUIRE(particleSim.positions_x[i] < 10);
		REQUIRE(particleSim.positions_y[i] > 0);
		REQUIRE(particleSim.positions_y[i] < 10);
		REQUIRE(particleSim.positions_z[i] > 0);
		REQUIRE(particleSim.positions_z[i] < 10);
	}

	auto old_grid_masses = grid.masses;

	transferAttributes(particleSim, grid);

	for (int i = 0 ; i < grid.masses.size() ; i ++)
	{
		if(grid.masses[i] == old_grid_masses[i])
			std::cout<<i<<std::endl;
		REQUIRE(grid.masses[i] != old_grid_masses[i]);
	}
}

TEST_CASE("Transfer Attributes grid 2 particles", "[Transfer][Attributes]")
{
	int nParticles = 5000;
	ParticleAttributes particleSim;
	randomisedParticleBB(particleSim, nParticles, 
						 0, 0, 0,
						 10, 10, 10);

	GridAttributes grid(10, 10, 10, 1);

	for (int i = 0 ; i < nParticles ; i ++)
	{
		REQUIRE(particleSim.positions_x[i] > 0);
		REQUIRE(particleSim.positions_x[i] < 10);
		REQUIRE(particleSim.positions_y[i] > 0);
		REQUIRE(particleSim.positions_y[i] < 10);
		REQUIRE(particleSim.positions_z[i] > 0);
		REQUIRE(particleSim.positions_z[i] < 10);
	}

	auto old_particle_masses = particleSim.masses;

	transferAttributes(particleSim, grid);

	transferAttributes(grid, particleSim);

	for (int i = 0 ; i < particleSim.masses.size() ; i ++)
	{
		if(particleSim.masses[i] == old_particle_masses[i])
			std::cout<<i<<std::endl;
		REQUIRE(particleSim.masses[i] != old_particle_masses[i]);
	}
}

TEST_CASE("timeStep", "[Transfer][Attributes]")
{
	int nParticles = 1;
	ParticleAttributes particleSim;
	randomisedParticleBB(particleSim, nParticles, 
						 0, 0, 0,
						 10, 10, 10);

	GridAttributes grid(10, 10, 10, 1);

	auto old_particle_masses = particleSim.masses;

		DEBUG();

	auto sphereCollision = 
	[](double pPos_x, double pPos_y, double pPos_z,
	   double pVel_x, double pVel_y, double pVel_z)
	-> std::array<double,3>
	{
		//=========================================Vector Funcs=========================================

		//==============================================================================================

		double sCenter_x = 5;
		double sCenter_y = 5;
		double sCenter_z = 5;

		double sphereRadius = 4;

		double t;

		std::array<double,3> newPos = add({pPos_x, pPos_y, pPos_z}, {pVel_x, pVel_y, pVel_z});

		auto intersected = intersectRaySphere({pPos_x, pPos_y, pPos_z}, {pVel_x, pVel_y, pVel_z},
			{sCenter_x, sCenter_y, sCenter_z}, sphereRadius, t, newPos);

		//std::cout<<"Intersection result: "<<intersected<<std::endl;
//
		//auto sRef_pPos_x = pPos_x - sCenter_x;
		//auto sRef_pPos_y = pPos_y - sCenter_y;
		//auto sRef_pPos_z = pPos_z - sCenter_z;
//
		//auto sRef_pNewPos_x = sRef_pPos_x + pVel_x;
		//auto sRef_pNewPos_y = sRef_pPos_y + pVel_y;
		//auto sRef_pNewPos_z = sRef_pPos_z + pVel_z;
//
		////Particle pos to circle center length
		//double p2cLen = length({sRef_pPos_x, sRef_pPos_y, sRef_pPos_z});
//
		////New particle pos to circle center length
		//double np2cLen = length({sRef_pNewPos_x, sRef_pNewPos_y, sRef_pNewPos_z});
//
		//double newPvel_x, newPvel_y, newPvel_z = 0;
//
		////Deal with case where particle is leaving circle
		//if(p2cLen < cRadius && cRadius < np2cLen)
		//{
		//	auto dx = pVel_x;
		//	auto dy = pVel_y;
		//	auto dz = pVel_z;
//
//
		//}
		return newPos;
	};

	for (int i = 0 ; i < 100 ; i ++)
	{
		DEBUG();
		transferAttributes(particleSim, grid);
		DEBUG();
		transferAttributes(grid, particleSim);
		DEBUG();
		timeStep(particleSim, sphereCollision, 1);
		DEBUG();
	}
}
