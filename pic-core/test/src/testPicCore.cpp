#include "common.h"

#include "picCore.h"

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
	std::cout<<"c000{\n"<<facePoses[0]<<"\n}"<<std::endl;
	std::cout<<"c100{\n"<<facePoses[1]<<"\n}"<<std::endl;
	std::cout<<"c001{\n"<<facePoses[2]<<"\n}"<<std::endl;
	std::cout<<"c101{\n"<<facePoses[3]<<"\n}"<<std::endl;
	std::cout<<"c010{\n"<<facePoses[4]<<"\n}"<<std::endl;
	std::cout<<"c110{\n"<<facePoses[5]<<"\n}"<<std::endl;
	std::cout<<"c011{\n"<<facePoses[6]<<"\n}"<<std::endl;
	std::cout<<"c111{\n"<<facePoses[7]<<"\n}"<<std::endl;

	auto cellFacePosc000_u = my_grid.cellMinFacePos_u(min_i, min_j, min_k);
	auto cellFacePosc111_u = my_grid.cellMinFacePos_u(max_i, max_j, max_k);
	std::cout<<"c000{\n"<<cellFacePosc000_u<<"\n}"<<std::endl;
	std::cout<<"c111{\n"<<cellFacePosc111_u<<"\n}"<<std::endl;

	
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

TEST_CASE("Particle to Grid attribute transfer", "[MacGrid][MacParticles][transferAttributes][p2g]"){
	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 1.f);
	std::vector<Vector3d> particlePoses{Vector3d(.51f, .76f, 1.f)};//, Vector3d(.75f, .75f, .75f)};
	std::vector<Vector3d> particleVels{Vector3d(.75f, 100.f, .75f)};//, Vector3d(.75f, .75f, .75f)};
	MacParticles my_particles(1.f, particlePoses, particleVels);

	transferAttributes(my_particles, my_grid);
	transferAttributes(my_grid, my_particles);

	for(int vi = 0 ; vi < my_particles.vel.size(); vi++ ){
		std::cout<<my_particles.vel[0]<<std::endl;
		std::cout<<particleVels[0]<<std::endl;
		REQUIRE(my_particles.vel[vi] == particleVels[vi]);
	}

}

TEST_CASE("Particle to Grid attribute transfer", "[MacGrid][MacParticles][transferAttributes][p2g]"){

	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 1.f);
	MacParticles my_particles(1.f, {}, {});
	int num_iterations = 100;
	for(int i=0 ; i<num_iterations ; i++){
		
		my_particles.positions = 
	}
	std::vector<Vector3d> particlePoses{Vector3d(.51f, .76f, 1.f)};//, Vector3d(.75f, .75f, .75f)};
	std::vector<Vector3d> particleVels{Vector3d(.75f, 100.f, .75f)};//, Vector3d(.75f, .75f, .75f)};
	


	transferAttributes(my_particles, my_grid);
	transferAttributes(my_grid, my_particles);

	for(int vi = 0 ; vi < my_particles.vel.size(); vi++ ){
		std::cout<<my_particles.vel[0]<<std::endl;
		std::cout<<particleVels[0]<<std::endl;
		REQUIRE(my_particles.vel[vi] == particleVels[vi]);
	}

}



//transferAttributes(MacParticles const & particles, MacGrid & grid)
//transferAttributes(MacGrid const & grid, MacParticles & particles)

//TEST_CASE("Staggered grid interpolation", "[MacGrid][getClosestNBRCellIdcs]")
//{
//	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 1.f);
//	tuple3i pCell{0,0,0};
//	Vector3d particlePos(0.75, 0.75, 0.75);
//	tuple3i pNBRCell = getClosestNBRCellIdcs(my_grid, particlePos, pCell);
//
//	REQUIRE(pNBRCell == tuple3i{1,1,1});
//
//	auto faceIdcs_u = getStageredCellFaceNBRIdcs_u(my_grid, pCell, pNBRCell);
//	auto faceIdcs_v = getStageredCellFaceNBRIdcs_v(my_grid, pCell, pNBRCell);
//	auto faceIdcs_w = getStageredCellFaceNBRIdcs_w(my_grid, pCell, pNBRCell);
//
//}

//TEST_CASE("computeAllCellCenterVel small", "[MacGrid][computeAllCellCenterVel]")
//{
//	MacGrid my_grid(Vector3d(0,0,0), 1, 1, 1, 1.f);
//
//	my_grid.cellFaceVel_u[my_grid.cellMinFaceIdx_u(0, 0, 0)] = -0.5;
//	my_grid.cellFaceVel_v[my_grid.cellMinFaceIdx_v(0, 0, 0)] = -0.5;
//	my_grid.cellFaceVel_w[my_grid.cellMinFaceIdx_w(0, 0, 0)] = -0.5;
//	my_grid.cellFaceVel_u[my_grid.cellMaxFaceIdx_u(0, 0, 0)] = 0.5;
//	my_grid.cellFaceVel_v[my_grid.cellMaxFaceIdx_v(0, 0, 0)] = 0.5;
//	my_grid.cellFaceVel_w[my_grid.cellMaxFaceIdx_w(0, 0, 0)] = 0.5;
//	DEBUG();
//	computeAllCellCenterVel(my_grid);
//	DEBUG();
//	std::vector<Vector3d> expectedCellCenterVel = {
//		Vector3d{0.f, 0.f, 0.f}
//	};
//	DEBUG();
//	for(int i = 0; i < my_grid.cellCenterVel.size() ; i++ ) 
//		REQUIRE(my_grid.cellCenterVel[i] == expectedCellCenterVel[i]);
//	DEBUG();
//	my_grid.cellFaceVel_u[my_grid.cellMinFaceIdx_u(0, 0, 0)] = 1;
//	my_grid.cellFaceVel_v[my_grid.cellMinFaceIdx_v(0, 0, 0)] = 1;
//	my_grid.cellFaceVel_w[my_grid.cellMinFaceIdx_w(0, 0, 0)] = 1;
//	my_grid.cellFaceVel_u[my_grid.cellMaxFaceIdx_u(0, 0, 0)] = 0.5;
//	my_grid.cellFaceVel_v[my_grid.cellMaxFaceIdx_v(0, 0, 0)] = 0.5;
//	my_grid.cellFaceVel_w[my_grid.cellMaxFaceIdx_w(0, 0, 0)] = 0.5;
//	computeAllCellCenterVel(my_grid);
//
//	expectedCellCenterVel = {
//		Vector3d{.75f, .75f, .75f}
//	};
//	REQUIRE(my_grid.cellCenterVel == expectedCellCenterVel);
//}
//
//TEST_CASE("computeAllCellCenterVel large", "[MacGrid][computeAllCellCenterVel]")
//{
//	MacGrid my_grid(Vector3d(0,0,0), 2, 1, 1, 1.f);
//
//	my_grid.cellFaceVel_u[my_grid.cellMinFaceIdx_u(0, 0, 0)] = -0.5;
//	my_grid.cellFaceVel_v[my_grid.cellMinFaceIdx_v(0, 0, 0)] = -0.5;
//	my_grid.cellFaceVel_w[my_grid.cellMinFaceIdx_w(0, 0, 0)] = -0.5;
//
//	my_grid.cellFaceVel_u[my_grid.cellMaxFaceIdx_u(0, 0, 0)] = 2.f;
//	my_grid.cellFaceVel_v[my_grid.cellMaxFaceIdx_v(0, 0, 0)] = 2.f;
//	my_grid.cellFaceVel_w[my_grid.cellMaxFaceIdx_w(0, 0, 0)] = 2.f;
//
//	my_grid.cellFaceVel_v[my_grid.cellMinFaceIdx_v(1, 0, 0)] = 1;
//	my_grid.cellFaceVel_w[my_grid.cellMinFaceIdx_w(1, 0, 0)] = 1;
//
//	my_grid.cellFaceVel_u[my_grid.cellMaxFaceIdx_u(1, 0, 0)] = 3.f;
//	my_grid.cellFaceVel_v[my_grid.cellMaxFaceIdx_v(1, 0, 0)] = 3.f;
//	my_grid.cellFaceVel_w[my_grid.cellMaxFaceIdx_w(1, 0, 0)] = 3.f;
//
//	computeAllCellCenterVel(my_grid);
//
//	std::vector<Vector3d> expectedCellCenterVel = {
//		Vector3d{.75f, .75f, .75f},
//		Vector3d{2.f, 2.f, 2.f}
//	};
//
//	REQUIRE(my_grid.cellCenterVel == expectedCellCenterVel);
//}


//void computeAllCellCenterVel(MacGrid & grid);
//void cellFaceVelTransfer_u(Vector3d const & pos, Vector3d const & vel, tuple8i const & cellUFaceNBRIdcs, MacGrid & grid);
//void cellFaceVelTransfer_v(Vector3d const & pos, Vector3d const & vel, tuple8i const & cellVFaceNBRIdcs, MacGrid & grid);
//void cellFaceVelTransfer_w(Vector3d const & pos, Vector3d const & vel, tuple8i const & cellWFaceNBRIdcs, MacGrid & grid);
//tuple3i getClosestNBRCellIdcs(MacGrid const & grid, Vector3d const & worldSpacePos);
//tuple8i getStageredCellFaceNBRIdcs_u(MacGrid const & grid, int min_i, int min_j, int min_k, int max_u, int max_j, int max_k);
//tuple8i getStageredCellFaceNBRIdcs_v(MacGrid const & grid, int min_i, int min_j, int min_k, int max_u, int max_j, int max_k);
//tuple8i getStageredCellFaceNBRIdcs_w(MacGrid const & grid, int min_i, int min_j, int min_k, int max_u, int max_j, int max_k);
//void transferAttributes(MacParticles const & particles, MacGrid & grid);
//void transferAttributes(MacGrid const & grid, MacParticles & particles);