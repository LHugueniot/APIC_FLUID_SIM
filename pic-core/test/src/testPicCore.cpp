#include "common.h"
#include "PicCore.h"

using namespace pic;

/**
void setDefaultCellStates(MacGrid & grid){

void setCellStates(MacGrid & grid, std::vector<tuple3i> & cellCoords, std::vector<MacGrid::CellState> & cellCenterStates){

void transferAttributes(Particles const & particles, MacGrid & grid){

void applyExternalForces(MacGrid & grid, double timeStep){

void initializeLaplacianNBRMat(MacGrid & grid){

void initializeCellCenterDivergence(MacGrid & grid){

inline double projectPressureOnFace(MacGrid & grid, int i, int j, int k, double cCP, double invCellSize){

void applyPressureForces(MacGrid & grid, double timeStep){

void cIsBound(MacGrid & grid, int i, int j, int k){

void enforceBoundaryVelocities(MacGrid & grid){

void extrapolateBoundaryVelocities(MacGrid & grid){

void enforceBoundaryPressure(MacGrid & grid){

void enforceBoundary(MacGrid & grid){

void collisionBasedParticleUpdate(Particles & particles, MacGrid & grid, double subStep){

void advanceStep(Particles & particles, MacGrid & grid, double timeStep){

void initializeLaplacianNBRMat(MacGrid & grid){

void applyPressureForces(MacGrid & grid, double subStep){
**/

TEST_CASE("PicCore getClosestNBRCellIdcs", "[MacGrid][getClosestNBRCellIdcs]"){
	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 1.f);
	tuple3i pCell{0,0,0};
	Vector3d particlePos(0.75, 0.75, 0.75);
	tuple3i pNBRCell = getClosestNBRCellIdcs(my_grid, particlePos, pCell[0], pCell[1], pCell[2]);
	REQUIRE(pNBRCell == tuple3i{1,1,1});
}

TEST_CASE("PicCore getLocalNBRIdcs", "[MacGrid][getLocalNBRIdcs]"){
	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 1.f);
	tuple3i cellCoord{0,0,0};
	Vector3d particlePos(.75f, .75f, .75f);
	tuple6i orderedNBRIdcs = getLocalNBRIdcs(my_grid, particlePos, cellCoord[0], cellCoord[1], cellCoord[2]);
	REQUIRE(orderedNBRIdcs == tuple6i{0,0,0,1,1,1});

	particlePos = {.25f, .25f, .25f};
	orderedNBRIdcs = getLocalNBRIdcs(my_grid, particlePos, cellCoord[0], cellCoord[1], cellCoord[2]);
	REQUIRE(orderedNBRIdcs == tuple6i{-1,-1,-1,0,0,0});

	particlePos = {.5f, .5f, .75f};
	orderedNBRIdcs = getLocalNBRIdcs(my_grid, particlePos, cellCoord[0], cellCoord[1], cellCoord[2]);
	REQUIRE(orderedNBRIdcs == tuple6i{0,0,0,1,1,1});

	particlePos = {.4999999f, .4999999f, .75f};
	orderedNBRIdcs = getLocalNBRIdcs(my_grid, particlePos, cellCoord[0], cellCoord[1], cellCoord[2]);
	REQUIRE(orderedNBRIdcs == tuple6i{-1,-1,0,0,0,1});

	cellCoord={1,0,0};
	particlePos = {1.1f, .4999999f, .75f};
	orderedNBRIdcs = getLocalNBRIdcs(my_grid, particlePos, cellCoord[0], cellCoord[1], cellCoord[2]);
	REQUIRE(orderedNBRIdcs == tuple6i{0,-1,0,1,0,1});

	cellCoord={1,1,1};
	particlePos = {1.1f, 1.5f, 1.99999f};
	orderedNBRIdcs = getLocalNBRIdcs(my_grid, particlePos, cellCoord[0], cellCoord[1], cellCoord[2]);
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

TEST_CASE("PicCore getStageredCellFaceNBRIdcs_u", "[MacGrid][getStageredCellFaceNBRIdcs_u][getStageredCellFaceNBRIdcs]"){
	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 1.f);

	//c000, c100, c001, c101, c010, c110, c011, c111
	int c000Idx = my_grid.cellFaceIdx<minU>(1,1,1);//c000
	int c100Idx = my_grid.cellFaceIdx<maxU>(1,1,1);//c100

	int c001Idx = my_grid.cellFaceIdx<minU>(1,1,2);//c001
	int c101Idx = my_grid.cellFaceIdx<maxU>(1,1,2);//c101

	int c010Idx = my_grid.cellFaceIdx<minU>(1,2,1);//c010
	int c110Idx = my_grid.cellFaceIdx<maxU>(1,2,1);//c110

	int c011Idx = my_grid.cellFaceIdx<minU>(1,2,2);//c011
	int c111Idx = my_grid.cellFaceIdx<maxU>(1,2,2);//c111

	
	my_grid.cellFaceVel_u[c000Idx] = 10;
	my_grid.cellFaceVel_u[c100Idx] = 10;
	my_grid.cellFaceVel_u[c001Idx] = 10;
	my_grid.cellFaceVel_u[c101Idx] = 10;
	my_grid.cellFaceVel_u[c010Idx] = 10;
	my_grid.cellFaceVel_u[c110Idx] = 10;
	my_grid.cellFaceVel_u[c011Idx] = 10;
	my_grid.cellFaceVel_u[c111Idx] = 10;

	Vector3d particlePos(1.75f, 1.75f, 1.75f);
	auto [i, j, k] = my_grid.gridCoord(particlePos);
	auto [min_i, min_j, min_k, max_i, max_j, max_k] = getLocalNBRIdcs(my_grid, particlePos, i, j, k);
	auto facePoses = getStageredCellFacePos_u(my_grid, min_i, min_j, min_k, max_i, max_j, max_k);

	REQUIRE(tuple6i{min_i, min_j, min_k, max_i, max_j, max_k} == tuple6i{1,1,1,2,2,2});

	auto faceIdcs_u = getStageredCellFaceNBRIdcs_u(my_grid, i, min_j, min_k, max_j, max_k);

	REQUIRE(faceIdcs_u == tuple8i{c000Idx, c100Idx, c001Idx, c101Idx, c010Idx, c110Idx, c011Idx, c111Idx});

	auto templateFaceIdcs_u = getStageredCellFaceNBRIdcs<minU>(my_grid, min_i, min_j, min_k, max_i, max_j, max_k);
	DEBUG_VAR(templateFaceIdcs_u);
	REQUIRE(faceIdcs_u == templateFaceIdcs_u);
}

TEST_CASE("PicCore getStageredCellFaceNBRIdcs_v", "[MacGrid][getStageredCellFaceNBRIdcs_v][getStageredCellFaceNBRIdcs]"){
	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 1.f);

	//c000, c100, c001, c101, c010, c110, c011, c111
	int c000Idx = my_grid.cellFaceIdx<minV>(1,1,1);//c000
	int c100Idx = my_grid.cellFaceIdx<minV>(2,1,1);//c100

	int c001Idx = my_grid.cellFaceIdx<minV>(1,1,2);//c001
	int c101Idx = my_grid.cellFaceIdx<minV>(2,1,2);//c101

	int c010Idx = my_grid.cellFaceIdx<maxV>(1,1,1);//c010
	int c110Idx = my_grid.cellFaceIdx<maxV>(2,1,1);//c110

	int c011Idx = my_grid.cellFaceIdx<maxV>(1,1,2);//c011
	int c111Idx = my_grid.cellFaceIdx<maxV>(2,1,2);//c111

	
	my_grid.cellFaceVel_v[c000Idx] = 10;
	my_grid.cellFaceVel_v[c100Idx] = 10;
	my_grid.cellFaceVel_v[c001Idx] = 10;
	my_grid.cellFaceVel_v[c101Idx] = 10;
	my_grid.cellFaceVel_v[c010Idx] = 10;
	my_grid.cellFaceVel_v[c110Idx] = 10;
	my_grid.cellFaceVel_v[c011Idx] = 10;
	my_grid.cellFaceVel_v[c111Idx] = 10;

	Vector3d particlePos(1.75f, 1.75f, 1.75f);
	auto [i, j, k] = my_grid.gridCoord(particlePos);
	auto [min_i, min_j, min_k, max_i, max_j, max_k] = getLocalNBRIdcs(my_grid, particlePos, i, j, k);

	REQUIRE(tuple6i{min_i, min_j, min_k, max_i, max_j, max_k} == tuple6i{1,1,1,2,2,2});

	auto faceIdcs_v = getStageredCellFaceNBRIdcs_v(my_grid, j, min_i, min_k, max_i, max_k);
	REQUIRE(faceIdcs_v == tuple8i{c000Idx, c100Idx, c001Idx, c101Idx, c010Idx, c110Idx, c011Idx, c111Idx});
	auto templateFaceIdcs_v = getStageredCellFaceNBRIdcs<minV>(my_grid, min_i, min_j, min_k, max_i, max_j, max_k);
	REQUIRE(faceIdcs_v == templateFaceIdcs_v);


}

TEST_CASE("PicCore getStageredCellFaceNBRIdcs_w", "[MacGrid][getStageredCellFaceNBRIdcs_w][getStageredCellFaceNBRIdcs]"){
	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 1.f);

	//c000, c100, c001, c101, c010, c110, c011, c111
	int c000Idx = my_grid.cellFaceIdx<minW>(1,1,1);//c000
	int c100Idx = my_grid.cellFaceIdx<minW>(2,1,1);//c100

	int c001Idx = my_grid.cellFaceIdx<maxW>(1,1,1);//c001
	int c101Idx = my_grid.cellFaceIdx<maxW>(2,1,1);//c101

	int c010Idx = my_grid.cellFaceIdx<minW>(1,2,1);//c010
	int c110Idx = my_grid.cellFaceIdx<minW>(2,2,1);//c110

	int c011Idx = my_grid.cellFaceIdx<maxW>(1,2,1);//c011
	int c111Idx = my_grid.cellFaceIdx<maxW>(2,2,1);//c111

	
	my_grid.cellFaceVel_w[c000Idx] = 10;
	my_grid.cellFaceVel_w[c100Idx] = 10;
	my_grid.cellFaceVel_w[c001Idx] = 10;
	my_grid.cellFaceVel_w[c101Idx] = 10;
	my_grid.cellFaceVel_w[c010Idx] = 10;
	my_grid.cellFaceVel_w[c110Idx] = 10;
	my_grid.cellFaceVel_w[c011Idx] = 10;
	my_grid.cellFaceVel_w[c111Idx] = 10;

	Vector3d particlePos(1.75f, 1.75f, 1.75f);
	auto [i, j, k] = my_grid.gridCoord(particlePos);
	auto [min_i, min_j, min_k, max_i, max_j, max_k] = getLocalNBRIdcs(my_grid, particlePos, i, j, k);

	REQUIRE(tuple6i{min_i, min_j, min_k, max_i, max_j, max_k} == tuple6i{1,1,1,2,2,2});

	auto faceIdcs_v = getStageredCellFaceNBRIdcs_w(my_grid, k, min_i, min_j, max_i, max_j);
	REQUIRE(faceIdcs_v == tuple8i{c000Idx, c100Idx, c001Idx, c101Idx, c010Idx, c110Idx, c011Idx, c111Idx});
	auto templateFaceIdcs_v = getStageredCellFaceNBRIdcs<minW>(my_grid, min_i, min_j, min_k, max_i, max_j, max_k);
	REQUIRE(faceIdcs_v == templateFaceIdcs_v);
}

TEST_CASE("PicCore cellFacePos", "[MacGrid][cellFacePos]"){
	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 1.f);

	REQUIRE(my_grid.cellFacePos<minU>(0,0,0) == Vector3d(0, .5f, .5f));
	REQUIRE(my_grid.cellFacePos<minV>(0,0,0) == Vector3d(.5f, 0, .5f));
	REQUIRE(my_grid.cellFacePos<minW>(0,0,0) == Vector3d(.5f, .5f, 0));

	REQUIRE(my_grid.cellFacePos<maxU>(0,0,0) == Vector3d(1, .5f, .5f));
	REQUIRE(my_grid.cellFacePos<maxV>(0,0,0) == Vector3d(.5f, 1, .5f));
	REQUIRE(my_grid.cellFacePos<maxW>(0,0,0) == Vector3d(.5f, .5f, 1));

	REQUIRE(my_grid.cellFacePos<minU>(1,0,0) == Vector3d(1, .5f, .5f));
	REQUIRE(my_grid.cellFacePos<minV>(1,0,0) == Vector3d(1.5f, 0, .5f));
	REQUIRE(my_grid.cellFacePos<minW>(1,0,0) == Vector3d(1.5f, .5f, 0));

	REQUIRE(my_grid.cellFacePos<maxU>(1,0,0) == Vector3d(2, .5f, .5f));
	REQUIRE(my_grid.cellFacePos<maxV>(1,0,0) == Vector3d(1.5f, 1, .5f));
	REQUIRE(my_grid.cellFacePos<maxW>(1,0,0) == Vector3d(1.5f, .5f, 1));

	REQUIRE(my_grid.cellFacePos<minU>(2,4,6) == Vector3d(2, 4.5f, 6.5f));
	REQUIRE(my_grid.cellFacePos<minV>(2,4,6) == Vector3d(2.5f, 4, 6.5f));
	REQUIRE(my_grid.cellFacePos<minW>(2,4,6) == Vector3d(2.5f, 4.5f, 6));

	REQUIRE(my_grid.cellFacePos<maxU>(2,4,6) == Vector3d(3, 4.5f, 6.5f));
	REQUIRE(my_grid.cellFacePos<maxV>(2,4,6) == Vector3d(2.5f, 5, 6.5f));
	REQUIRE(my_grid.cellFacePos<maxW>(2,4,6) == Vector3d(2.5f, 4.5f, 7));

	MacGrid my_grid2(Vector3d(0,0,0), 10, 10, 10, 1.5f);

	REQUIRE(my_grid2.cellFacePos<minU>(0,0,0) == Vector3d(0, .75f, .75f));
	REQUIRE(my_grid2.cellFacePos<minV>(0,0,0) == Vector3d(.75f, 0, .75f));
	REQUIRE(my_grid2.cellFacePos<minW>(0,0,0) == Vector3d(.75f, .75f, 0));

	REQUIRE(my_grid2.cellFacePos<maxU>(0,0,0) == Vector3d(1.5f, .75f, .75f));
	REQUIRE(my_grid2.cellFacePos<maxV>(0,0,0) == Vector3d(.75f, 1.5f, .75f));
	REQUIRE(my_grid2.cellFacePos<maxW>(0,0,0) == Vector3d(.75f, .75f, 1.5f));
}

TEST_CASE("PicCore getStageredCellFaceNBRIdcs", "[MacGrid][cellFacePos][getStageredCellFaceNBRIdcs]"){
	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 1.f);
	Vector3d pPos(1.94f, 7.47f, 1.94f);

	auto [i, j, k] = my_grid.gridCoord(pPos);
	REQUIRE(i == 1);
	REQUIRE(j == 7);
	REQUIRE(k == 1);

	auto localNBRIdcs = getLocalNBRIdcs(my_grid, pPos, i, j, k);

	REQUIRE(localNBRIdcs == tuple6i{1,6,1,2,7,2});

	auto cellFacePosc000_u = my_grid.cellFacePos<minU>(i,6,1);
	auto cellFacePosc111_u = my_grid.cellFacePos<maxU>(i,7,2);

	REQUIRE(cellFacePosc000_u == Vector3d(1.f, 6.5f, 1.5f));
	REQUIRE(cellFacePosc111_u == Vector3d(2.f, 7.5f, 2.5f));

	DEBUG_VAR(cellFacePosc000_u);
	DEBUG_VAR(cellFacePosc111_u);

	auto cellFacePosc000_v = my_grid.cellFacePos<minV>(1,j,1);
	auto cellFacePosc111_v = my_grid.cellFacePos<maxV>(2,j,2);

	REQUIRE(cellFacePosc000_v == Vector3d(1.5f, 7.f, 1.5f));
	REQUIRE(cellFacePosc111_v == Vector3d(2.5f, 8.f, 2.5f));

	DEBUG_VAR(cellFacePosc000_v);
	DEBUG_VAR(cellFacePosc111_v);

	auto cellFacePosc000_w = my_grid.cellFacePos<minW>(1,6,k);
	auto cellFacePosc111_w = my_grid.cellFacePos<maxW>(2,7,k);

	REQUIRE(cellFacePosc000_w == Vector3d(1.5f, 6.5f, 1.f));
	REQUIRE(cellFacePosc111_w == Vector3d(2.5f, 7.5f, 2.f));

	DEBUG_VAR(cellFacePosc000_w);
	DEBUG_VAR(cellFacePosc111_w);

}

TEST_CASE("Particle to Grid attribute transfer", "[MacGrid][Particles][transferAttributes][p2g][g2p]"){
	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 1.f);
	std::vector<double> particlePoses{1.51f, 1.76f, 2.f};//, Vector3d(.75f, .75f, .75f)};
	std::vector<double> particleVels{.75f, 100.f, .75f};//, Vector3d(.75f, .75f, .75f)};
	Particles my_particles(1.f, particlePoses, particleVels);

	transferAttributes(my_particles, my_grid);
	REQUIRE(my_particles.velocities == particleVels);
	transferAttributes(my_grid, my_particles);

	REQUIRE(my_particles.velocities.size() == particleVels.size());

	for(int vi = 0 ; vi < my_particles.velocities.size(); vi++ ){
		std::cout<<my_particles.velocities[vi]<<std::endl;
		std::cout<<particleVels[vi]<<std::endl;
		REQUIRE(my_particles.velocities[vi] == particleVels[vi]);
	}

}

TEST_CASE("Particle to Grid attribute transfer Random", "[MacGrid][Particles][transferAttributes][p2g][g2p][random]"){

    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_real_distribution<> x_dist(0, 10);
    std::uniform_real_distribution<> y_dist(0, 10);
    std::uniform_real_distribution<> z_dist(0, 10);


	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 1.f);
	Particles my_particles(1.f, std::vector<double>{}, std::vector<double>{});
	int num_iterations = 100;

	for(int i=0 ; i<num_iterations ; i++){
		std::vector<double> particlePoses{ x_dist(eng), y_dist(eng), z_dist(eng) };
		std::vector<double> particleVels{ x_dist(eng), y_dist(eng), z_dist(eng) };

		my_particles.positions = particlePoses;
		my_particles.velocities = particleVels;

		transferAttributes(my_particles, my_grid);
		transferAttributes(my_grid, my_particles);

		//std::cout<<my_particles.velocities[0]<<std::endl;
		//std::cout<<particleVels[0]<<std::endl;
		REQUIRE(my_particles.velocities == particleVels);
	}

}

TEST_CASE("Calculate time sub step", "[MacGrid][Particles][calculateSubStep]"){
	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 2.f);
	std::vector<double> particlePoses{3.f, 3.f, 4.f};
	std::vector<double> particleVels{.75f, 100.f, .75f};
	Particles my_particles(1.f, particlePoses, particleVels);

	double timeStep = 1.f / 60.f;

	transferAttributes(my_particles, my_grid);
	double subStep = calculateSubStep(my_grid, timeStep);

	REQUIRE(particleVels[0] * subStep <= my_grid.cellSize);
	REQUIRE(particleVels[1] * subStep <= my_grid.cellSize);
	REQUIRE(particleVels[2] * subStep <= my_grid.cellSize);
}

TEST_CASE("Grid applyExternalForces", "[MacGrid][Particles][applyExternalForces]"){
	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 2.f);
	std::vector<double> particlePoses{3.f, 3.f, 4.f};
	std::vector<double> particleVels{.75f, 100.f, .75f};
	Particles my_particles(1.f, particlePoses, particleVels);

	double timeStep = 1.f / 60.f;
	transferAttributes(my_particles, my_grid);
	applyExternalForces(my_grid, timeStep);
	transferAttributes(my_grid, my_particles);

	REQUIRE(my_particles.velocities[0] == particleVels[0]);
	REQUIRE(my_particles.velocities[1] == particleVels[1] - 9.8f);
	REQUIRE(my_particles.velocities[2] == particleVels[2]);
}

TEST_CASE("Grid initializeLaplacianNBRMat", "[MacGrid][initializeLaplacianNBRMat]"){
	MacGrid my_grid(Vector3d(0,0,0), 5, 5, 5, 2.f);
	auto rows = my_grid.laplacianSparseNBR.rows();
	auto cols = my_grid.laplacianSparseNBR.cols();
	DEBUG_VAR(rows);
	DEBUG_VAR(cols);
	initializeLaplacianNBRMat(my_grid);
	for(int i=0; i<5 ; i++)
		for(int j=0; j<5 ; j++)
			for(int k=0; k<5 ; k++){
				auto idx = my_grid.cellCenterIdx(i,j,k);
				if(my_grid.cellIsInBounds(i+1,j,k) && !my_grid.isSolidCell(i+1,j,k)){
					REQUIRE(my_grid.laplacianSparseNBR.coeffRef(idx, my_grid.cellCenterIdx(i+1,j,k)) == 1);
					REQUIRE(my_grid.laplacianSparseNBR.coeffRef(my_grid.cellCenterIdx(i+1,j,k), idx) == 1);
				}
				if(my_grid.cellIsInBounds(i-1,j,k) && !my_grid.isSolidCell(i-1,j,k)){
					REQUIRE(my_grid.laplacianSparseNBR.coeffRef(idx, my_grid.cellCenterIdx(i-1,j,k)) == 1);
					REQUIRE(my_grid.laplacianSparseNBR.coeffRef(my_grid.cellCenterIdx(i-1,j,k), idx) == 1);
				}
				if(my_grid.cellIsInBounds(i,j+1,k) && !my_grid.isSolidCell(i,j+1,k)){
					REQUIRE(my_grid.laplacianSparseNBR.coeffRef(idx, my_grid.cellCenterIdx(i,j+1,k)) == 1);
					REQUIRE(my_grid.laplacianSparseNBR.coeffRef(my_grid.cellCenterIdx(i,j+1,k), idx) == 1);
				}
				if(my_grid.cellIsInBounds(i,j-1,k) && !my_grid.isSolidCell(i,j-1,k)){
					REQUIRE(my_grid.laplacianSparseNBR.coeffRef(idx, my_grid.cellCenterIdx(i,j-1,k)) == 1);
					REQUIRE(my_grid.laplacianSparseNBR.coeffRef(my_grid.cellCenterIdx(i,j-1,k), idx) == 1);
				}
				if(my_grid.cellIsInBounds(i,j,k+1) && !my_grid.isSolidCell(i,j,k+1)){
					REQUIRE(my_grid.laplacianSparseNBR.coeffRef(idx, my_grid.cellCenterIdx(i,j,k+1)) == 1);
					REQUIRE(my_grid.laplacianSparseNBR.coeffRef(my_grid.cellCenterIdx(i,j,k+1), idx) == 1);
				}
				if(my_grid.cellIsInBounds(i,j,k-1) && !my_grid.isSolidCell(i,j,k-1)){
					REQUIRE(my_grid.laplacianSparseNBR.coeffRef(idx, my_grid.cellCenterIdx(i,j,k-1)) == 1);
					REQUIRE(my_grid.laplacianSparseNBR.coeffRef(my_grid.cellCenterIdx(i,j,k-1), idx) == 1);
				}
				REQUIRE(my_grid.laplacianSparseNBR.coeffRef(idx, idx) <= 0);
			}
	//DEBUG_VAR(my_grid.laplacianSparseNBR);
	REQUIRE(&my_grid.cg != nullptr);
}


TEST_CASE("Grid initializeCellCenterDivergence",
	"[MacGrid][Particles][initializeCellCenterDivergence]"){
	MacGrid my_grid(Vector3d(0,0,0), 5, 5, 5, 1.f);
	//my_grid.printOutMaxFaceVels();
	//initializeLaplacianNBRMat(my_grid);

	std::vector<double> particlePositions = AABBRandomParticles(Vector3d(2,2,2), Vector3d(2,4,2), 10);
	std::vector<double> particleVelocities(particlePositions.size(), 0);

	Particles my_particles(1.f, particlePositions, particleVelocities);

	double timeStep = 1.f / 24.f;

	transferAttributes(my_particles, my_grid);
	applyExternalForces(my_grid, timeStep);
	enforceBoundaryVelocities(my_grid);
	initializeCellCenterDivergence(my_grid,timeStep);
	DEBUG();
	//i++;

	DEBUG_VAR(my_grid.cellFaceVel_u);
	DEBUG_VAR(my_grid.cellFaceVel_w);
	DEBUG();
}

/*
TEST_CASE("Grid initializeCellCenterDivergence",
	"[MacGrid][Particles][initializeCellCenterDivergence]"){
	MacGrid my_grid(Vector3d(0,0,0), 5, 5, 5, 1.f);
	my_grid.printOutMaxFaceVels();
	initializeLaplacianNBRMat(my_grid);

	std::vector<double> particlePoses{3.f, 3.f, 3.f};
	std::vector<double> particleVels{0.f, 0.f, 0.f};
	Particles my_particles(1.f, particlePoses, particleVels);
	
	double timeStep = 1.f / 24.f;
	
	int i = 0;
	int max = 1;
	while( i < max){
	
		transferAttributes(my_particles, my_grid);
		applyExternalForces(my_grid, timeStep);
		//enforceBoundaryVelocities(my_grid);
		//initializeCellCenterDivergence(my_grid);
		DEBUG();
		i++;
	}
	DEBUG_VAR(my_grid.cellFaceVel_u);
	DEBUG_VAR(my_grid.cellFaceVel_w);
	DEBUG();
}
*/

TEST_CASE("Grid applyPressureForces", "[MacGrid][Particles][applyPressureForces]"){
	MacGrid my_grid(Vector3d(0,0,0), 5, 5, 5, 1.f);
	my_grid.printOutMaxFaceVels();
	initializeLaplacianNBRMat(my_grid);

	std::vector<double> particlePoses{3.5f, 3.5f, 3.5f};
	std::vector<double> particleVels{0.f, 0.f, 0.f};
	Particles my_particles(1.f, particlePoses, particleVels);

	double timeStep = 1.f / 24.f;

	int i = 0;
	int max = 1;
	while( i < max){

		transferAttributes(my_particles, my_grid);
		applyExternalForces(my_grid, timeStep);
		my_grid.printOutMaxFaceVels();
		applyPressureForces(my_grid, timeStep);
		transferAttributes(my_grid, my_particles);
		i++;
	}
	//DEBUG_VAR(my_grid.cellFaceVel_u);
	//DEBUG_VAR(my_grid.cellFaceVel_w);
	DEBUG();
}

TEST_CASE("Generate random particles in bounding box particles", "[Particles][AABBRandomParticles]"){

	std::vector<double> particlePositions = AABBRandomParticles(Vector3d(0,0,0), Vector3d(10,10,10), 10);

	REQUIRE(particlePositions.size() == 30);
	for(auto &pPos : particlePositions){
		REQUIRE(0 < pPos);
		REQUIRE(pPos < 10);
	}
}

TEST_CASE("Set boundary cells to solid", "[MacGrid][setDefaultCellStates]"){

    pic::MacGrid my_grid(Vector3d(0,0,0), 100, 100, 100, 1.f);
    setDefaultCellStates(my_grid);

	for(int i=0; i<my_grid.cellNum_i ; i++)
		for(int j=0; j<my_grid.cellNum_j ; j++)
			for(int k=0; k<my_grid.cellNum_k ; k++){
				if(	(i == 0) || (i == my_grid.cellNum_i - 1) 
				 || (j == 0) || (j == my_grid.cellNum_j - 1)
				 || (k == 0) || (k == my_grid.cellNum_k - 1))
				REQUIRE(my_grid.isSolidCell(i,j,k));
			}
}

TEST_CASE("Set enforce boundary pressure based on neighbourhood.", "[MacGrid][enforceBoundaryPressure]"){

    pic::MacGrid my_grid(Vector3d(0,0,0), 5, 5, 5, 1.f);
    //setDefaultCellStates(my_grid);
	for (uint i = 0 ; i < my_grid.cellCenterPressure.rows() ; i++ )
		my_grid.cellCenterPressure.coeffRef(i, 0) = 2;
	enforceBoundaryPressure(my_grid);

	auto gridCCP = [&my_grid](int i, int j, int k){
		auto idx = my_grid.cellCenterIdx(i,j,k);
		return my_grid.cellCenterPressure.coeffRef(idx, 0);
	};
	REQUIRE(gridCCP(1,1,1) == 4);
	REQUIRE(gridCCP(3,1,1) == 4);
	REQUIRE(gridCCP(3,3,1) == 4);
	REQUIRE(gridCCP(3,3,3) == 4);
	REQUIRE(gridCCP(1,3,3) == 4);
	REQUIRE(gridCCP(1,1,3) == 4);
	REQUIRE(gridCCP(3,1,3) == 4);
	REQUIRE(gridCCP(1,3,1) == 4);

	REQUIRE(gridCCP(1,1,2) == 3);
	REQUIRE(gridCCP(2,1,1) == 3);
	REQUIRE(gridCCP(2,1,3) == 3);
	REQUIRE(gridCCP(3,1,2) == 3);
	REQUIRE(gridCCP(1,3,2) == 3);
	REQUIRE(gridCCP(2,3,1) == 3);
	REQUIRE(gridCCP(2,3,3) == 3);
	REQUIRE(gridCCP(3,3,2) == 3);


	REQUIRE(gridCCP(2,1,2) == 2.4);

	REQUIRE(gridCCP(2,2,1) == 2.4);
	REQUIRE(gridCCP(1,2,2) == 2.4);
	REQUIRE(gridCCP(2,2,3) == 2.4);
	REQUIRE(gridCCP(3,2,2) == 2.4);

	REQUIRE(gridCCP(2,3,2) == 2.4);

}


TEST_CASE("Update Particles multiple times","[Particles][updateParticles]"){

	std::vector<double> particlePoses{
		0,0,0,
		10,10,10,
		1,2,3
	};
	std::vector<double> particleVels{
		20.f, 13.5f, 0.f,
		12.f, 6.f, 6.f,
		5.f, 5.f, 4.f
	};

	Particles my_particles(1.f, particlePoses, particleVels);

	REQUIRE(my_particles.num == 3);
	updateParticles(my_particles, 1.f/2.f);

	std::vector<double> expectedPositions{
		10.f, 6.75f, 0.f,
		16.f,13.f,13.f,
		3.5f,4.5f,5.f
		};
	REQUIRE(my_particles.positions == expectedPositions);
}

TEST_CASE("update particles", "[MacGrid][updateParticles][1000it]"){

	std::cout<<"===================================Update \
Particles==================================="<<std::endl;


    pic::MacGrid my_grid(Vector3d(0,0,0), 5, 5, 5, 1.f);
    initializeLaplacianNBRMat(my_grid);

	std::vector<double> particlePoses = pic::AABBRandomParticles(
		Vector3d(2, 3 ,2), Vector3d(2, 4 ,2), 1);
    std::vector<double> particleVels(particlePoses.size(), 0.f);

	double timeStep = 1.f/48.f;

	Particles my_particles(1.f, particlePoses, particleVels);
	int iterations = 100;

	for (int i = 0 ; i < iterations ; i++){
		advanceStep(my_particles, my_grid, timeStep);
	}
}

TEST_CASE("Enforce the Boundary Velocities", "[MacGrid][enforceBoundaryVelocities]"){

    pic::MacGrid my_grid(Vector3d(0,0,0), 5, 5, 5, 1.f);
    applyExternalForces(my_grid, 1.f);
    for(auto & vel_u : my_grid.cellFaceVel_u)
    	vel_u = 10;
    my_grid.printOutMaxFaceVels();

    enforceBoundaryVelocities(my_grid);
    DEBUG();

    my_grid.printOutMaxFaceVels();


	for(int i=0; i<my_grid.cellNum_i ; i++)
		for(int j=0; j<my_grid.cellNum_j ; j++)
			for(int k=0; k<my_grid.cellNum_k ; k++){
	
				if(i == 0)
					REQUIRE(my_grid.cellFaceVel<maxU>(i,j,k) >= 0);
				else if(i == my_grid.cellNum_i - 1)
					REQUIRE(my_grid.cellFaceVel<minU>(i,j,k) <= 0);
				if(j == 0)
					REQUIRE(my_grid.cellFaceVel<maxV>(i,j,k) >= 0);
				else if(j == my_grid.cellNum_j - 1)
					REQUIRE(my_grid.cellFaceVel<minV>(i,j,k) <= 0);
				if(k == 0)
					REQUIRE(my_grid.cellFaceVel<maxW>(i,j,k) >= 0);
				else if(k == my_grid.cellNum_k - 1)
					REQUIRE(my_grid.cellFaceVel<minW>(i,j,k) <= 0);
			}

}


TEST_CASE("Extrapolate the Boundary Velocities", "[MacGrid][extrapolateBoundaryVelocities]"){

    pic::MacGrid my_grid(Vector3d(0,0,0), 5, 5, 5, 1.f);
    //applyExternalForces(my_grid, 1.f);

    for(int i = 0 ; i < my_grid.cellFaceVel_v.size() ; i++)
    	my_grid.cellFaceVel_v[i] = -i;

    my_grid.printOutMaxFaceVels();

    enforceBoundaryVelocities(my_grid);
    my_grid.printOutMaxFaceVels();
    extrapolateBoundaryVelocities(my_grid);
    my_grid.printOutMaxFaceVels();

}


TEST_CASE("Test iterating over grid boundary", "[MacGrid][iterateBoundary]"){

    pic::MacGrid my_grid(Vector3d(0,0,0), 5, 5, 5, 1.f);
	std::fill(my_grid.cellCenterState.begin(), my_grid.cellCenterState.end(), SOLID);

	std::cout<<"cell state"<<std::endl;
	for ( uint j = 0 ; j < my_grid.cellNum_j ; j++){
		for ( uint i = 0 ; i < my_grid.cellNum_i ; i++){
			for ( uint k = 0 ; k < my_grid.cellNum_k ; k++){
				std::cout<<my_grid.getCellState(i,j,k)<<" ";
			}
			std::cout<<std::endl;
		}
		std::cout<<std::endl;
	}

	for (uint i = 0 ; i < my_grid.cellNum_i ; i += (my_grid.cellNum_i - 1))
		for (uint j = 1 ; j < my_grid.cellNum_j - 1 ; j ++)
			for (uint k = 1 ; k < my_grid.cellNum_k - 1 ; k++)
				my_grid.cellCenterState.at(my_grid.cellCenterIdx(i, j, k)) = AIR;
	
	for (uint i = 1 ; i < my_grid.cellNum_i - 1 ; i ++)
		for (uint j = 0 ; j < my_grid.cellNum_j ; j += (my_grid.cellNum_j - 1))
			for (uint k = 1 ; k < my_grid.cellNum_k - 1 ; k ++)
				my_grid.cellCenterState.at(my_grid.cellCenterIdx(i, j, k)) = AIR;
	
	for (uint i = 1 ; i < my_grid.cellNum_i - 1 ; i ++)
		for (uint j = 1 ; j < my_grid.cellNum_j - 1 ; j++)
			for (uint k = 0 ; k < my_grid.cellNum_k ; k += (my_grid.cellNum_k - 1))
				my_grid.cellCenterState.at(my_grid.cellCenterIdx(i, j, k)) = AIR;

/** Edge Cell states

	for (uint i = 0 ; i < my_grid.cellNum_i ; i += (my_grid.cellNum_i - 1))
		for (uint j = 0 ; j < my_grid.cellNum_j ; j += (my_grid.cellNum_j - 1))
			for (uint k = 0 ; k < my_grid.cellNum_k ; k++)
				my_grid.cellCenterState.at(my_grid.cellCenterIdx(i, j, k)) = AIR;
	
	for (uint i = 0 ; i < my_grid.cellNum_i ; i ++)
		for (uint j = 0 ; j < my_grid.cellNum_j ; j += (my_grid.cellNum_j - 1))
			for (uint k = 0 ; k < my_grid.cellNum_k ; k += (my_grid.cellNum_k - 1))
				my_grid.cellCenterState.at(my_grid.cellCenterIdx(i, j, k)) = AIR;
	
	for (uint i = 0 ; i < my_grid.cellNum_i ; i += (my_grid.cellNum_i - 1))
		for (uint j = 0 ; j < my_grid.cellNum_j ; j++)
			for (uint k = 0 ; k < my_grid.cellNum_k ; k += (my_grid.cellNum_k - 1))
				my_grid.cellCenterState.at(my_grid.cellCenterIdx(i, j, k)) = AIR;
**/

	std::cout<<"cell state"<<std::endl;
	for ( uint j = 0 ; j < my_grid.cellNum_j ; j++){
		for ( uint i = 0 ; i < my_grid.cellNum_i ; i++){
			for ( uint k = 0 ; k < my_grid.cellNum_k ; k++){
				std::cout<<my_grid.getCellState(i,j,k)<<" ";
			}
			std::cout<<std::endl;
		}
		std::cout<<std::endl;
	}
}

TEST_CASE("Generate uniform particles in a cubic bounding box", "[Particles][AABCubeUniformParticles]"){

    std::vector<double> particlePositions = pic::AABCubeUniformParticles(Vector3d(2.5 , 8 ,2.5 ), 2, 125);

	REQUIRE(particlePositions.size() == 125 * 3);
	for (uint i = 0 ; i < particlePositions.size() ; i ++){
		DEBUG_VAR(i);
		REQUIRE(0 < particlePositions[i]);
		REQUIRE(particlePositions[i] < 2.f);
	}

}


/**
To test:
void setCellStates(MacGrid & grid, std::vector<tuple3i> & cellCoords, std::vector<MacGrid::CellState> & cellCenterStates);

**/

//TEST_CASE("Calculate time sub step", "[MacGrid][calculateSubStep]"){
//	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 2.f);
//	std::vector<Vector3d> particlePoses{Vector3d(.51f, .76f, 1.f)};//, Vector3d(.75f, .75f, .75f)};
//	std::vector<Vector3d> particleVels{Vector3d(.75f, 100.f, .75f)};//, Vector3d(.75f, .75f, .75f)};
//	Particles my_particles(1.f, particlePoses, particleVels);
//	
//	transferAttributes(my_particles, my_grid);
//
//
//	REQUIRE(particleVels[0][0] * subStep <= my_grid.cellSize);
//	REQUIRE(particleVels[0][1] * subStep <= my_grid.cellSize);
//	REQUIRE(particleVels[0][2] * subStep <= my_grid.cellSize);
//
//}
