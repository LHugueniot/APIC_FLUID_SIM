#include "common.h"
#include "PicCore.h"

using namespace pic;

TEST_CASE("PicCore getClosestNBRCellIdcs", "[MacGrid][getClosestNBRCellIdcs]"){
	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 1.f);
	tuple3i pCell{0,0,0};
	Vector3d particlePos(0.75, 0.75, 0.75);
	tuple3i pNBRCell = getClosestNBRCellIdcs(my_grid, particlePos, pCell[0], pCell[1], pCell[2]);
	REQUIRE(pNBRCell == tuple3i{1,1,1});
}

TEST_CASE("PicCore getOrderedNBRIdcs", "[MacGrid][getOrderedNBRIdcs]"){
	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 1.f);
	tuple3i cellCoord{0,0,0};
	Vector3d particlePos(.75f, .75f, .75f);
	tuple6i orderedNBRIdcs = getOrderedNBRIdcs(my_grid, particlePos, cellCoord[0], cellCoord[1], cellCoord[2]);
	REQUIRE(orderedNBRIdcs == tuple6i{0,0,0,1,1,1});

	particlePos = {.25f, .25f, .25f};
	orderedNBRIdcs = getOrderedNBRIdcs(my_grid, particlePos, cellCoord[0], cellCoord[1], cellCoord[2]);
	REQUIRE(orderedNBRIdcs == tuple6i{-1,-1,-1,0,0,0});

	particlePos = {.5f, .5f, .75f};
	orderedNBRIdcs = getOrderedNBRIdcs(my_grid, particlePos, cellCoord[0], cellCoord[1], cellCoord[2]);
	REQUIRE(orderedNBRIdcs == tuple6i{0,0,0,1,1,1});

	particlePos = {.4999999f, .4999999f, .75f};
	orderedNBRIdcs = getOrderedNBRIdcs(my_grid, particlePos, cellCoord[0], cellCoord[1], cellCoord[2]);
	REQUIRE(orderedNBRIdcs == tuple6i{-1,-1,0,0,0,1});

	cellCoord={1,0,0};
	particlePos = {1.1f, .4999999f, .75f};
	orderedNBRIdcs = getOrderedNBRIdcs(my_grid, particlePos, cellCoord[0], cellCoord[1], cellCoord[2]);
	REQUIRE(orderedNBRIdcs == tuple6i{0,-1,0,1,0,1});

	cellCoord={1,1,1};
	particlePos = {1.1f, 1.5f, 1.99999f};
	orderedNBRIdcs = getOrderedNBRIdcs(my_grid, particlePos, cellCoord[0], cellCoord[1], cellCoord[2]);
	REQUIRE(orderedNBRIdcs == tuple6i{0,1,1,1,2,2});
}

std::array<Vector3d, 8> getStageredCellFacePos_u(MacGrid const & grid, int min_i, int min_j, int min_k, int max_i, int max_j, int max_k){
	return {grid.cellMinFacePos_u(min_i, min_j, min_k), grid.cellMinFacePos_u(max_i, min_j, min_k),
			grid.cellMinFacePos_u(min_i, min_j, max_k), grid.cellMinFacePos_u(max_i, min_j, max_k),
			grid.cellMinFacePos_u(min_i, max_j, min_k), grid.cellMinFacePos_u(max_i, max_j, min_k),
			grid.cellMinFacePos_u(min_i, max_j, max_k), grid.cellMinFacePos_u(max_i, max_j, max_k)};
}

std::array<Vector3d, 8> getStageredCellFacePos_v(MacGrid const & grid, int min_i, int min_j, int min_k, int max_i, int max_j, int max_k){
	return {grid.cellMinFacePos_v(min_i, min_j, min_k), grid.cellMinFacePos_v(max_i, min_j, min_k),
			grid.cellMinFacePos_v(min_i, min_j, max_k), grid.cellMinFacePos_v(max_i, min_j, max_k),
			grid.cellMinFacePos_v(min_i, max_j, min_k), grid.cellMinFacePos_v(max_i, max_j, min_k),
			grid.cellMinFacePos_v(min_i, max_j, max_k), grid.cellMinFacePos_v(max_i, max_j, max_k)};
}

std::array<Vector3d, 8> getStageredCellFacePos_w(MacGrid const & grid, int min_i, int min_j, int min_k, int max_i, int max_j, int max_k){
	return {grid.cellMinFacePos_w(min_i, min_j, min_k), grid.cellMinFacePos_w(max_i, min_j, min_k),
			grid.cellMinFacePos_w(min_i, min_j, max_k), grid.cellMinFacePos_w(max_i, min_j, max_k),
			grid.cellMinFacePos_w(min_i, max_j, min_k), grid.cellMinFacePos_w(max_i, max_j, min_k),
			grid.cellMinFacePos_w(min_i, max_j, max_k), grid.cellMinFacePos_w(max_i, max_j, max_k)};
}

TEST_CASE("PicCore getStageredCellFaceNBRIdcs_u", "[MacGrid][getStageredCellFaceNBRIdcs_u]"){
	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 1.f);

	//c000, c100, c001, c101, c010, c110, c011, c111
	int c000Idx = my_grid.cellMinFaceIdx_u(0,0,0);//c000
	int c100Idx = my_grid.cellMaxFaceIdx_u(0,0,0);//c100

	int c001Idx = my_grid.cellMinFaceIdx_u(0,0,1);//c001
	int c101Idx = my_grid.cellMaxFaceIdx_u(0,0,1);//c101

	int c010Idx = my_grid.cellMinFaceIdx_u(0,1,0);//c010
	int c110Idx = my_grid.cellMaxFaceIdx_u(0,1,0);//c110

	int c011Idx = my_grid.cellMinFaceIdx_u(0,1,1);//c011
	int c111Idx = my_grid.cellMaxFaceIdx_u(0,1,1);//c111

	
	my_grid.cellFaceVel_u[c000Idx] = 10;
	my_grid.cellFaceVel_u[c100Idx] = 10;
	my_grid.cellFaceVel_u[c001Idx] = 10;
	my_grid.cellFaceVel_u[c101Idx] = 10;
	my_grid.cellFaceVel_u[c010Idx] = 10;
	my_grid.cellFaceVel_u[c110Idx] = 10;
	my_grid.cellFaceVel_u[c011Idx] = 10;
	my_grid.cellFaceVel_u[c111Idx] = 10;

	tuple3i cellCoord{0,0,0};
	Vector3d particlePos(.75f, .75f, .75f);
	auto [min_i, min_j, min_k, max_i, max_j, max_k] = getOrderedNBRIdcs(my_grid, particlePos, cellCoord[0], cellCoord[1], cellCoord[2]);
	auto facePoses = getStageredCellFacePos_u(my_grid, min_i, min_j, min_k, max_i, max_j, max_k);

	
	REQUIRE(tuple6i{min_i, min_j, min_k, max_i, max_j, max_k} == tuple6i{0,0,0,1,1,1});

	auto faceIdcs_u = getStageredCellFaceNBRIdcs_u(my_grid, min_i, min_j, min_k, max_i, max_j, max_k);
	REQUIRE(faceIdcs_u == tuple8i{c000Idx, c100Idx, c001Idx, c101Idx, c010Idx, c110Idx, c011Idx, c111Idx});


}

TEST_CASE("PicCore getStageredCellFaceNBRIdcs_v", "[MacGrid][getStageredCellFaceNBRIdcs_v]"){
	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 1.f);

	//c000, c100, c001, c101, c010, c110, c011, c111
	int c000Idx = my_grid.cellMinFaceIdx_v(0,0,0);//c000
	int c100Idx = my_grid.cellMinFaceIdx_v(1,0,0);//c100

	int c001Idx = my_grid.cellMinFaceIdx_v(0,0,1);//c001
	int c101Idx = my_grid.cellMinFaceIdx_v(1,0,1);//c101

	int c010Idx = my_grid.cellMaxFaceIdx_v(0,0,0);//c010
	int c110Idx = my_grid.cellMaxFaceIdx_v(1,0,0);//c110

	int c011Idx = my_grid.cellMaxFaceIdx_v(0,0,1);//c011
	int c111Idx = my_grid.cellMaxFaceIdx_v(1,0,1);//c111

	
	my_grid.cellFaceVel_v[c000Idx] = 10;
	my_grid.cellFaceVel_v[c100Idx] = 10;
	my_grid.cellFaceVel_v[c001Idx] = 10;
	my_grid.cellFaceVel_v[c101Idx] = 10;
	my_grid.cellFaceVel_v[c010Idx] = 10;
	my_grid.cellFaceVel_v[c110Idx] = 10;
	my_grid.cellFaceVel_v[c011Idx] = 10;
	my_grid.cellFaceVel_v[c111Idx] = 10;

	tuple3i cellCoord{0,0,0};
	Vector3d particlePos(.75f, .75f, .75f);
	auto [min_i, min_j, min_k, max_i, max_j, max_k] = getOrderedNBRIdcs(my_grid, particlePos, cellCoord[0], cellCoord[1], cellCoord[2]);
	REQUIRE(tuple6i{min_i, min_j, min_k, max_i, max_j, max_k} == tuple6i{0,0,0,1,1,1});

	auto faceIdcs_v = getStageredCellFaceNBRIdcs_v(my_grid, min_i, min_j, min_k, max_i, max_j, max_k);
	REQUIRE(faceIdcs_v == tuple8i{c000Idx, c100Idx, c001Idx, c101Idx, c010Idx, c110Idx, c011Idx, c111Idx});


}

TEST_CASE("PicCore getStageredCellFaceNBRIdcs_w", "[MacGrid][getStageredCellFaceNBRIdcs_w]"){
	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 1.f);

	//c000, c100, c001, c101, c010, c110, c011, c111
	int c000Idx = my_grid.cellMinFaceIdx_w(0,0,0);//c000
	int c100Idx = my_grid.cellMinFaceIdx_w(1,0,0);//c100

	int c001Idx = my_grid.cellMinFaceIdx_w(0,0,1);//c001
	int c101Idx = my_grid.cellMinFaceIdx_w(1,0,1);//c101

	int c010Idx = my_grid.cellMinFaceIdx_w(0,1,0);//c010
	int c110Idx = my_grid.cellMinFaceIdx_w(1,1,0);//c110

	int c011Idx = my_grid.cellMinFaceIdx_w(0,1,1);//c011
	int c111Idx = my_grid.cellMinFaceIdx_w(1,1,1);//c111

	
	my_grid.cellFaceVel_w[c000Idx] = 10;
	my_grid.cellFaceVel_w[c100Idx] = 10;
	my_grid.cellFaceVel_w[c001Idx] = 10;
	my_grid.cellFaceVel_w[c101Idx] = 10;
	my_grid.cellFaceVel_w[c010Idx] = 10;
	my_grid.cellFaceVel_w[c110Idx] = 10;
	my_grid.cellFaceVel_w[c011Idx] = 10;
	my_grid.cellFaceVel_w[c111Idx] = 10;

	tuple3i cellCoord{0,0,0};
	Vector3d particlePos(.75f, .75f, .75f);
	auto [min_i, min_j, min_k, max_i, max_j, max_k] = getOrderedNBRIdcs(my_grid, particlePos, cellCoord[0], cellCoord[1], cellCoord[2]);
	REQUIRE(tuple6i{min_i, min_j, min_k, max_i, max_j, max_k} == tuple6i{0,0,0,1,1,1});

	auto faceIdcs_v = getStageredCellFaceNBRIdcs_w(my_grid, min_i, min_j, min_k, max_i, max_j, max_k);
	REQUIRE(faceIdcs_v == tuple8i{c000Idx, c100Idx, c001Idx, c101Idx, c010Idx, c110Idx, c011Idx, c111Idx});
}

TEST_CASE("Particle to Grid attribute transfer", "[MacGrid][MacParticles][transferAttributes][p2g][g2p]"){
	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 1.f);
	std::vector<Vector3d> particlePoses{Vector3d(.51f, .76f, 1.f)};//, Vector3d(.75f, .75f, .75f)};
	std::vector<Vector3d> particleVels{Vector3d(.75f, 100.f, .75f)};//, Vector3d(.75f, .75f, .75f)};
	MacParticles my_particles(1.f, particlePoses, particleVels);

	transferAttributes(my_particles, my_grid);
	transferAttributes(my_grid, my_particles);

	for(int vi = 0 ; vi < my_particles.velocities.size(); vi++ ){
		std::cout<<my_particles.velocities[0]<<std::endl;
		std::cout<<particleVels[0]<<std::endl;
		REQUIRE(my_particles.velocities[vi] == particleVels[vi]);
	}

}

TEST_CASE("Particle to Grid attribute transfer Random", "[MacGrid][MacParticles][transferAttributes][p2g][g2p][random]"){

    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_real_distribution<> x_dist(0, 10);
    std::uniform_real_distribution<> y_dist(0, 10);
    std::uniform_real_distribution<> z_dist(0, 10);


	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 1.f);
	MacParticles my_particles(1.f, {}, {});
	int num_iterations = 100;

	for(int i=0 ; i<num_iterations ; i++){
		std::vector<Vector3d> particlePoses{ Vector3d(x_dist(eng), y_dist(eng), z_dist(eng)) };
		std::vector<Vector3d> particleVels{ Vector3d(x_dist(eng), y_dist(eng), z_dist(eng)) };

		my_particles.positions = particlePoses;
		my_particles.velocities = particleVels;

		transferAttributes(my_particles, my_grid);
		transferAttributes(my_grid, my_particles);

		//std::cout<<my_particles.velocities[0]<<std::endl;
		//std::cout<<particleVels[0]<<std::endl;
		REQUIRE(my_particles.velocities[0] == particleVels[0]);
	}

}


TEST_CASE("Calculate time sub step", "[MacGrid][MacParticles][calculateSubStep]"){
	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 2.f);
	std::vector<Vector3d> particlePoses{Vector3d(3.f, 3.f, 4.f)};//, Vector3d(.75f, .75f, .75f)};
	std::vector<Vector3d> particleVels{Vector3d(.75f, 100.f, .75f)};//, Vector3d(.75f, .75f, .75f)};
	MacParticles my_particles(1.f, particlePoses, particleVels);

	double timeStep = 1.f / 60.f;

	transferAttributes(my_particles, my_grid);
	double subStep = calculateSubStep(my_grid, timeStep);

	REQUIRE(particleVels[0][0] * subStep <= my_grid.cellSize);
	REQUIRE(particleVels[0][1] * subStep <= my_grid.cellSize);
	REQUIRE(particleVels[0][2] * subStep <= my_grid.cellSize);
}

TEST_CASE("Grid applyExternalForces", "[MacGrid][MacParticles][applyExternalForces]"){
	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 2.f);
	std::vector<Vector3d> particlePoses{Vector3d(3.f, 3.f, 4.f)};//, Vector3d(.75f, .75f, .75f)};
	std::vector<Vector3d> particleVels{Vector3d(.75f, 100.f, .75f)};//, Vector3d(.75f, .75f, .75f)};
	MacParticles my_particles(1.f, particlePoses, particleVels);

	double subStep = 1.f / 60.f;
	transferAttributes(my_particles, my_grid);
	applyExternalForces(my_grid, subStep);
	transferAttributes(my_grid, my_particles);

	REQUIRE(my_particles.velocities[0][0] == particleVels[0][0]);
	REQUIRE(my_particles.velocities[0][1] == particleVels[0][1] - 9.8f * subStep);
	REQUIRE(my_particles.velocities[0][2] == particleVels[0][2]);
}

TEST_CASE("Grid initializeLaplacianNBRMat", "[MacGrid][initializeLaplacianNBRMat]"){
	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 2.f);
	initializeLaplacianNBRMat(my_grid);

	//std::cout<<my_grid.laplacianSparseNBR<<std::endl;
}

TEST_CASE("Grid applyPressureForces", "[MacGrid][applyPressureForces]"){
	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 2.f);
	initializeLaplacianNBRMat(my_grid);

	std::vector<Vector3d> particlePoses{Vector3d(3.f, 3.f, 4.f)};//, Vector3d(.75f, .75f, .75f)};
	std::vector<Vector3d> particleVels{Vector3d(.75f, 100.f, .75f)};//, Vector3d(.75f, .75f, .75f)};
	MacParticles my_particles(1.f, particlePoses, particleVels);

	double subStep = 1.f / 60.f;
	transferAttributes(my_particles, my_grid);
	applyExternalForces(my_grid, subStep);
	applyPressureForces(my_grid, subStep);
	transferAttributes(my_grid, my_particles);
}

//TEST_CASE("Calculate time sub step", "[MacGrid][calculateSubStep]"){
//	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 2.f);
//	std::vector<Vector3d> particlePoses{Vector3d(.51f, .76f, 1.f)};//, Vector3d(.75f, .75f, .75f)};
//	std::vector<Vector3d> particleVels{Vector3d(.75f, 100.f, .75f)};//, Vector3d(.75f, .75f, .75f)};
//	MacParticles my_particles(1.f, particlePoses, particleVels);
//	
//	transferAttributes(my_particles, my_grid);
//
//
//	REQUIRE(particleVels[0][0] * subStep <= my_grid.cellSize);
//	REQUIRE(particleVels[0][1] * subStep <= my_grid.cellSize);
//	REQUIRE(particleVels[0][2] * subStep <= my_grid.cellSize);
//
//}
