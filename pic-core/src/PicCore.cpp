#include "PicCore.h"

namespace pic{

//=============================================================================================
//=====================================HELPER FUNCTIONS========================================
//=============================================================================================

//genereate randome particles in an Axis aligned bounding box
std::vector<double> AABBRandomParticles(Vector3d const & fieldC000, 
	Vector3d const & fieldC111, uint particleSystemSize)
{
    std::vector<double> particleSystemVertexData(particleSystemSize * 3);

    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_real_distribution<> dist_x(fieldC000[0], fieldC111[0]);
    std::uniform_real_distribution<> dist_y(fieldC000[1], fieldC111[1]);
    std::uniform_real_distribution<> dist_z(fieldC000[2], fieldC111[2]);

    for (uint i = 0 ; i < particleSystemSize ; i ++)
    {
        particleSystemVertexData[i * 3] = dist_x(eng);
        particleSystemVertexData[i * 3 + 1] = dist_y(eng);
        particleSystemVertexData[i * 3 + 2] = dist_z(eng);
    }
    return particleSystemVertexData;
}

std::vector<double> AABCubeUniformParticles(Vector3d const & fieldC000, 
	Vector3d const & fieldC111, double interParticleDistance){
	assert(interParticleDistance != 0);
	uint particleSystemSize = 
		((fieldC111[0] - fieldC000[0])/interParticleDistance + 1) *
		((fieldC111[1] - fieldC000[1])/interParticleDistance + 1) *
		((fieldC111[2] - fieldC000[2])/interParticleDistance + 1);

    std::vector<double> particleSystemVertexData(particleSystemSize * 3);

	uint idx = 0;

    for (double i = fieldC000[0] ; i < fieldC111[0] ; i += interParticleDistance)
    	for (double j = fieldC000[1] ; j < fieldC111[1] ; j += interParticleDistance)
    		for (double k = fieldC000[2] ; k < fieldC111[2] ; k += interParticleDistance){

				particleSystemVertexData.at(idx) = i;
				particleSystemVertexData.at(idx + 1) = j;
				particleSystemVertexData.at(idx + 2) = k;
				idx+=3;
    		}
    if (idx != particleSystemSize)
    	particleSystemVertexData.resize(idx);
    return particleSystemVertexData;
}

tuple3i getClosestNBRCellIdcs(MacGrid const & grid,
	Vector3d const & worldSpacePos, int i, int j, int k){

	Vector3d const cellSpacePos = grid.cellSpacePos(worldSpacePos, i, j, k);
	return {
		cellSpacePos[0] < .5f ? i - 1 : i + 1,
		cellSpacePos[1] < .5f ? j - 1 : j + 1,
		cellSpacePos[2] < .5f ? k - 1 : k + 1
	};
}

tuple6i getCellNBRIdcs(int i, int j, int k){
	return {i - 1, j - 1, k - 1, i + 1, j + 1, k + 1};
}

tuple6i getLocalNBRIdcs(MacGrid const & grid, Vector3d const & worldSpacePos, int i, int j, int k){
	Vector3d const cellSpacePos = grid.cellSpacePos(worldSpacePos, i, j, k);

	DEBUG_VAR("cellSpacePos");
	DEBUG_VAR(i);
	DEBUG_VAR(j);
	DEBUG_VAR(k);

	DEBUG_VAR(cellSpacePos);
	int nbr_i = cellSpacePos[0] < .5f ? i - 1 : i + 1;
	int nbr_j = cellSpacePos[1] < .5f ? j - 1 : j + 1;
	int nbr_k = cellSpacePos[2] < .5f ? k - 1 : k + 1;
	return {std::min(i, nbr_i), std::min(j, nbr_j), std::min(k, nbr_k), 
			std::max(i, nbr_i), std::max(j, nbr_j), std::max(k, nbr_k)};
}

////Get cell face Indices
//template<FaceDim D>
//tuple8i getStageredCellFaceNBRIdcs(MacGrid const & grid, 
//	int min_i, int min_j, int min_k, 
//	int max_i, int max_j, int max_k){
//
//	return {grid.cellFaceIdx<D>(min_i, min_j, min_k), grid.cellFaceIdx<D>(max_i, min_j, min_k),
//			grid.cellFaceIdx<D>(min_i, min_j, max_k), grid.cellFaceIdx<D>(max_i, min_j, max_k),
//			grid.cellFaceIdx<D>(min_i, max_j, min_k), grid.cellFaceIdx<D>(max_i, max_j, min_k),
//			grid.cellFaceIdx<D>(min_i, max_j, max_k), grid.cellFaceIdx<D>(max_i, max_j, max_k)};
//}


//c000, c100, c001, c101, c010, c110, c011, c111

tuple8i getStageredCellFaceNBRIdcs_u(MacGrid const & grid,
	int i, int min_j, int min_k, int max_j, int max_k){

	return {grid.cellFaceIdx<minU>(i, min_j, min_k), grid.cellFaceIdx<maxU>(i, min_j, min_k),
			grid.cellFaceIdx<minU>(i, min_j, max_k), grid.cellFaceIdx<maxU>(i, min_j, max_k),
			grid.cellFaceIdx<minU>(i, max_j, min_k), grid.cellFaceIdx<maxU>(i, max_j, min_k),
			grid.cellFaceIdx<minU>(i, max_j, max_k), grid.cellFaceIdx<maxU>(i, max_j, max_k)};
}

tuple8i getStageredCellFaceNBRIdcs_v(MacGrid const & grid,
	int j, int min_i, int min_k, int max_i, int max_k){

	return {grid.cellFaceIdx<minV>(min_i, j, min_k), grid.cellFaceIdx<minV>(max_i, j, min_k),
			grid.cellFaceIdx<minV>(min_i, j, max_k), grid.cellFaceIdx<minV>(max_i, j, max_k),
			grid.cellFaceIdx<maxV>(min_i, j, min_k), grid.cellFaceIdx<maxV>(max_i, j, min_k),
			grid.cellFaceIdx<maxV>(min_i, j, max_k), grid.cellFaceIdx<maxV>(max_i, j, max_k)};
}

tuple8i getStageredCellFaceNBRIdcs_w(MacGrid const & grid,
	int k, int min_i, int min_j, int max_i, int max_j){

	return {grid.cellFaceIdx<minW>(min_i, min_j, k), grid.cellFaceIdx<minW>(max_i, min_j, k), 
			grid.cellFaceIdx<maxW>(min_i, min_j, k), grid.cellFaceIdx<maxW>(max_i, min_j, k),
			grid.cellFaceIdx<minW>(min_i, max_j, k), grid.cellFaceIdx<minW>(max_i, max_j, k),
			grid.cellFaceIdx<maxW>(min_i, max_j, k), grid.cellFaceIdx<maxW>(max_i, max_j, k)};
}

inline std::array<Vector3d, 8> getStageredCellFacePos_u(MacGrid const & grid,
	int i, int min_j, int min_k, int max_j, int max_k){

	return {grid.cellFacePos<minU>(i, min_j, min_k), grid.cellFacePos<maxU>(i, min_j, min_k),
			grid.cellFacePos<minU>(i, min_j, max_k), grid.cellFacePos<maxU>(i, min_j, max_k),
			grid.cellFacePos<minU>(i, max_j, min_k), grid.cellFacePos<maxU>(i, max_j, min_k),
			grid.cellFacePos<minU>(i, max_j, max_k), grid.cellFacePos<maxU>(i, max_j, max_k)};
}

inline std::array<Vector3d, 8> getStageredCellFacePos_v(MacGrid const & grid,
	int j, int min_i, int min_k, int max_i, int max_k){

	return {grid.cellFacePos<minV>(min_i, j, min_k), grid.cellFacePos<minV>(max_i, j, min_k),
			grid.cellFacePos<minV>(min_i, j, max_k), grid.cellFacePos<minV>(max_i, j, max_k),
			grid.cellFacePos<maxV>(min_i, j, min_k), grid.cellFacePos<maxV>(max_i, j, min_k),
			grid.cellFacePos<maxV>(min_i, j, max_k), grid.cellFacePos<maxV>(max_i, j, max_k)};
}

inline std::array<Vector3d, 8> getStageredCellFacePos_w(MacGrid const & grid,
	int k, int min_i, int min_j, int max_i, int max_j){

	return {grid.cellFacePos<minW>(min_i, min_j, k), grid.cellFacePos<minW>(max_i, min_j, k), 
			grid.cellFacePos<maxW>(min_i, min_j, k), grid.cellFacePos<maxW>(max_i, min_j, k),
			grid.cellFacePos<minW>(min_i, max_j, k), grid.cellFacePos<minW>(max_i, max_j, k),
			grid.cellFacePos<maxW>(min_i, max_j, k), grid.cellFacePos<maxW>(max_i, max_j, k)};
}

void setDefaultCellStates(MacGrid & grid){

	for (uint i = 0 ; i < grid.cellNum_i ; i += grid.cellNum_i - 1)
		for (uint j = 0 ; j < grid.cellNum_j ; j++)
			for (uint k = 0 ; k < grid.cellNum_k ; k++)
				grid.cellCenterState.at(grid.cellCenterIdx(i, j, k)) = SOLID;

	for (uint i = 0 ; i < grid.cellNum_i ; i ++)
		for (uint j = 0 ; j < grid.cellNum_j ; j += grid.cellNum_j - 1)
			for (uint k = 0 ; k < grid.cellNum_k ; k++)
				grid.cellCenterState.at(grid.cellCenterIdx(i, j, k)) = SOLID;

	for (uint i = 0 ; i < grid.cellNum_i ; i ++)
		for (uint j = 0 ; j < grid.cellNum_j ; j++)
			for (uint k = 0 ; k < grid.cellNum_k ; k += grid.cellNum_k - 1)
				grid.cellCenterState.at(grid.cellCenterIdx(i, j, k)) = SOLID;
}

void setCellStates(MacGrid & grid, std::vector<tuple3i> & cellCoords,
	std::vector<CellState> & cellCenterStates){

	assert(cellCoords.size() == cellCenterStates.size());
	for(uint i=0 ; i<cellCoords.size() ; i++)
		grid.cellCenterState[
			grid.cellCenterIdx(cellCoords[i][0], cellCoords[i][1], cellCoords[i][2])] = cellCenterStates[i];
}

//=============================================================================================
//=====================================FLUID SIM FUNCTIONS=====================================
//=============================================================================================

void applyExternalForces(MacGrid & grid, double timeStep){

	for ( uint i = 0 ; i < grid.cellNum_i; i++)
		for ( uint j = 0 ; j < grid.cellNum_j ; j++)
			for ( uint k = 0 ; k < grid.cellNum_k ; k++)
				if(grid.isValidFace<maxV>(i, j, k))
					grid.cellFaceVel<maxV>(i, j, k) += GRAV_Y * timeStep;
			//}
			//{
				//if( (grid.getCellState(i,j,k) == LIQUID) && grid.isValidFace<minV>(i, j, k))
}

inline double projectPressureOnFace(MacGrid & grid, int idxNBR, double cCP){

	double nbrCCP = grid.cellCenterPressure.at(idxNBR);
	return (cCP - nbrCCP);
}

void applyPressureForces(MacGrid & grid, double timeStep){

	std::fill(grid.cellCenterPressure.begin(), grid.cellCenterPressure.end(), 0);

	SparseMatrix<double> laplacianSparseNBR(grid.cellNum_i * grid.cellNum_j * grid.cellNum_k, grid.cellNum_i * grid.cellNum_j * grid.cellNum_k);
	SparseVector<double> cellCenterPressure(grid.cellNum_i * grid.cellNum_j * grid.cellNum_k);
	SparseMatrix<double> cellCenterDivergence(grid.cellNum_i * grid.cellNum_j * grid.cellNum_k, 1);

	laplacianSparseNBR.setZero();
	cellCenterDivergence.setZero();

	double scale = timeStep/ (grid.cellSize * grid.cellSize);
	//std::cout<<"SCALE: "<<scale<<std::endl;
	for ( uint i = 0 ; i < grid.cellNum_i - 0 ; i++)
		for ( uint j = 0 ; j < grid.cellNum_j - 0 ; j++)
			for ( uint k = 0 ; k < grid.cellNum_k - 0 ; k++){
				DEBUG();

				int totalValidNeighbours = 0;

				auto cellCenterIdx = grid.cellCenterIdx(i, j, k);
				if(grid.cellCenterState.at(cellCenterIdx) != LIQUID)
					continue;

				//Calculate Divergence of each cell
				auto [minFaceIdx_u, minFaceIdx_v, minFaceIdx_w,
					maxFaceIdx_u, maxFaceIdx_v, maxFaceIdx_w] =
					grid.cellFaceIdcs(i, j, k);

				DEBUG();


					if(!grid.isSolidCell(i - 1, j, k) && grid.isCellInGrid_i(i - 1)){
						totalValidNeighbours += 1;
						if(grid.getCellState(i - 1, j, k) == LIQUID){
							auto cellCenterMinIdx_i = grid.cellCenterIdx(i - 1, j, k);
							laplacianSparseNBR.coeffRef(cellCenterIdx, cellCenterMinIdx_i) = -scale;
						}
					}
				DEBUG();

					if(!grid.isSolidCell(i + 1, j, k) && grid.isCellInGrid_i(i + 1)){
						totalValidNeighbours += 1;
						if(grid.getCellState(i + 1, j, k) == LIQUID){
							auto cellCenterMaxIdx_i = grid.cellCenterIdx(i + 1, j, k);
							laplacianSparseNBR.coeffRef(cellCenterIdx, cellCenterMaxIdx_i) = -scale;
						}
					}

				DEBUG();

					if(!grid.isSolidCell(i, j - 1, k) && grid.isCellInGrid_j(j - 1)){
						totalValidNeighbours += 1;
						if(grid.getCellState(i, j - 1, k) == LIQUID){
							auto cellCenterMinIdx_j = grid.cellCenterIdx(i, j - 1, k);
							laplacianSparseNBR.coeffRef(cellCenterIdx, cellCenterMinIdx_j) = -scale;
						}
					}
				DEBUG();

					if(!grid.isSolidCell(i, j + 1, k) && grid.isCellInGrid_j(j + 1)){
						totalValidNeighbours += 1;
						if(grid.getCellState(i, j + 1, k) == LIQUID){
							auto cellCenterMaxIdx_j = grid.cellCenterIdx(i, j + 1, k);
							laplacianSparseNBR.coeffRef(cellCenterIdx, cellCenterMaxIdx_j) = -scale;
						}
					}
				DEBUG();

					if(!grid.isSolidCell(i, j, k - 1) && grid.isCellInGrid_k(k - 1)){
						totalValidNeighbours += 1;
						if(grid.getCellState(i, j, k - 1) == LIQUID){
							auto cellCenterMinIdx_k = grid.cellCenterIdx(i, j, k - 1);
							laplacianSparseNBR.coeffRef(cellCenterIdx, cellCenterMinIdx_k) = -scale;
						}
					}
				DEBUG();

					if(!grid.isSolidCell(i, j, k + 1) && grid.isCellInGrid_j(k + 1)){
						totalValidNeighbours += 1;
						if(grid.getCellState(i, j, k + 1) == LIQUID){
							auto cellCenterMaxIdx_k = grid.cellCenterIdx(i, j, k + 1);
							laplacianSparseNBR.coeffRef(cellCenterIdx, cellCenterMaxIdx_k) = -scale;
						}
					}
				//DEBUG();

				laplacianSparseNBR.coeffRef(cellCenterIdx, cellCenterIdx) = totalValidNeighbours * scale;


				double div = 
					-(grid.cellFaceVel_u.at(maxFaceIdx_u) - grid.cellFaceVel_u.at(minFaceIdx_u))/grid.cellSize +
					-(grid.cellFaceVel_v.at(maxFaceIdx_v) - grid.cellFaceVel_v.at(minFaceIdx_v))/grid.cellSize +
					-(grid.cellFaceVel_w.at(maxFaceIdx_w) - grid.cellFaceVel_w.at(minFaceIdx_w))/grid.cellSize;

				if(i == 0 || grid.isSolidCell(i - 1, j, k))
					div -= grid.cellFaceVel_u.at(minFaceIdx_u) /grid.cellSize;

				if(i == grid.cellNum_i - 1 || grid.isSolidCell(i + 1, j, k))
					div += grid.cellFaceVel_u.at(maxFaceIdx_u) /grid.cellSize;

				if(j == 0 || grid.isSolidCell(i, j - 1, k))
					div -= grid.cellFaceVel_v.at(minFaceIdx_v) /grid.cellSize;

				if(j == grid.cellNum_j - 1 || grid.isSolidCell(i, j + 1, k))
					div += grid.cellFaceVel_v.at(maxFaceIdx_v) /grid.cellSize;

				if(k == 0 || grid.isSolidCell(i, j, k - 1))
					div -= grid.cellFaceVel_w.at(minFaceIdx_w) /grid.cellSize;

				if(k == grid.cellNum_k - 1 || grid.isSolidCell(i, j, k + 1))
					div += grid.cellFaceVel_w.at(maxFaceIdx_w) /grid.cellSize;

				cellCenterDivergence.coeffRef(cellCenterIdx, 0) = div;
			}
	//checkForNan(grid.cellCenterDivergence);

	grid.cg.compute(laplacianSparseNBR);
	//std::cout<<cellCenterDivergence<<std::endl;

	cellCenterPressure.setZero();
	
	//initializeCellCenterDivergence(grid, timeStep);
	cellCenterPressure = grid.cg.solve(cellCenterDivergence);

	if(ENFORCE_INCOMPRESSIBILITY){
		for(Eigen::SparseVector<double>::InnerIterator it(cellCenterPressure); it; ++it)
    	    if (grid.cellCenterState.at(it.index()) == LIQUID)
    	    	grid.cellCenterPressure.at(it.index()) = it.value() * (double)grid.cellCenterDensity.at(it.index()) / WATER_PARTICLE_DENSITY;
	}
    else{
		for(Eigen::SparseVector<double>::InnerIterator it(cellCenterPressure); it; ++it)
    	    if (grid.cellCenterState.at(it.index()) == LIQUID)
    	    	grid.cellCenterPressure.at(it.index()) = it.value();
    }
    //enforceBoundaryPressure(grid);

	scale = timeStep/grid.cellSize;

	for ( uint i = 0 ; i < grid.cellNum_i - 0 ; i++)
		for ( uint j = 0 ; j < grid.cellNum_j - 0 ; j++)
			for ( uint k = 0 ; k < grid.cellNum_k - 0 ; k++){

				auto cellCenterIdx = grid.cellCenterIdx(i, j, k);
				auto idxMin_i = grid.cellCenterIdx(i - 1, j, k);
				auto idxMin_j = grid.cellCenterIdx(i, j - 1, k);
				auto idxMin_k = grid.cellCenterIdx(i, j, k - 1);

				auto cellState = grid.cellCenterState.at(cellCenterIdx);

				double & cCP = grid.cellCenterPressure.at(cellCenterIdx);

    			if(grid.isCellInGrid_i(i - 1)){

					auto cellStateMin_i = grid.cellCenterState.at(idxMin_i);

        			bool leftExists = i > 0;
        			bool rightExists = i < grid.cellNum_i;
        			bool leftFluid = leftExists && cellStateMin_i == LIQUID;
        			bool rightFluid = rightExists && cellState == LIQUID;

 					if (leftFluid || rightFluid)
    					//grid.cellFaceVel_u.at(grid.cellFaceIdx<minU>(i, j, k)) -= projectPressureOnFace(grid, idxMin_i, cCP) * scale;// *grid.cellSize;
    					grid.cellFaceVel<minU>(i, j, k) -= projectPressureOnFace(grid, idxMin_i, cCP) * scale;// *grid.cellSize;
    			}
			
    			if(grid.isCellInGrid_j(j - 1)){ 
					auto cellStateMin_j = grid.cellCenterState.at(idxMin_j);

        			bool leftExists = j > 0;
        			bool rightExists = j < grid.cellNum_j;
        			bool leftFluid = leftExists && cellStateMin_j == LIQUID;
        			bool rightFluid = rightExists && cellState == LIQUID;

 					if (leftFluid || rightFluid)
    					//grid.cellFaceVel_v.at(grid.cellFaceIdx<minV>(i, j, k)) -= projectPressureOnFace(grid, idxMin_j, cCP) * scale;// *grid.cellSize;
    					grid.cellFaceVel<minV>(i, j, k) -= projectPressureOnFace(grid, idxMin_j, cCP) * scale;// *grid.cellSize;
				}
    			if(grid.isCellInGrid_k(k - 1)){
					auto cellStateMin_k = grid.cellCenterState.at(idxMin_k);

        			bool leftExists = k > 0;
        			bool rightExists = k < grid.cellNum_k;
        			bool leftFluid = leftExists && cellStateMin_k == LIQUID;
        			bool rightFluid = rightExists && cellState == LIQUID;

 					if (leftFluid || rightFluid)
    					//grid.cellFaceVel_w.at(grid.cellFaceIdx<minW>(i, j, k)) -= projectPressureOnFace(grid, idxMin_k, cCP) * scale;// *grid.cellSize;
    					grid.cellFaceVel<minW>(i, j, k) -= projectPressureOnFace(grid, idxMin_k, cCP) * scale;// *grid.cellSize;
    			}
			}
}

void enforceBoundaryVelocities(MacGrid & grid){

	for ( uint i = 0 ; i < grid.cellNum_i ; i++)
		for ( uint j = 0 ; j < grid.cellNum_j ; j++)
			for ( uint k = 0 ; k < grid.cellNum_k ; k++){
				
				//std::cout<<"enforceBoundaryVelocities start"<<std::endl;
				if(grid.getCellState(i,j,k) != LIQUID)
					continue;
				//std::cout<<"grid.getCellState(i,j,k) != SOLID"<<std::endl;

				if(grid.isValidFace<minU>(i, j, k) && grid.isSolidCell(i - 1, j, k))
					grid.cellFaceVel<minU>(i, j, k) = std::max<double>(0.f, grid.cellFaceVel<minU>(i, j, k));

				if(grid.isValidFace<maxU>(i, j, k) && grid.isSolidCell(i + 1, j, k))
					grid.cellFaceVel<maxU>(i, j, k) = std::min<double>(0.f, grid.cellFaceVel<maxU>(i, j, k));

				if(grid.isValidFace<minV>(i, j, k) && grid.isSolidCell(i, j - 1, k))
					grid.cellFaceVel<minV>(i, j, k) = std::max<double>(0.f, grid.cellFaceVel<minV>(i, j, k));

				if(grid.isValidFace<maxV>(i, j, k) && grid.isSolidCell(i, j + 1, k))
					grid.cellFaceVel<maxV>(i, j, k) = std::min<double>(0.f, grid.cellFaceVel<maxV>(i, j, k));

				if(grid.isValidFace<minW>(i, j, k) && grid.isSolidCell(i, j, k - 1))
					grid.cellFaceVel<minW>(i, j, k) = std::max<double>(0.f, grid.cellFaceVel<minW>(i, j, k));

				if(grid.isValidFace<maxW>(i, j, k) && grid.isSolidCell(i, j, k + 1))
					grid.cellFaceVel<maxW>(i, j, k) = std::min<double>(0.f, grid.cellFaceVel<maxW>(i, j, k));
			}
}

//void enforceBoundaryVelocities(MacGrid & grid){
//
//	for ( uint i = 0 ; i < grid.cellNum_i ; i++)
//		for ( uint j = 0 ; j < grid.cellNum_j ; j++)
//			for ( uint k = 0 ; k < grid.cellNum_k ; k++){
//				
//				//std::cout<<"enforceBoundaryVelocities start"<<std::endl;
//				if(grid.getCellState(i,j,k) != LIQUID)
//					continue;
//				//std::cout<<"grid.getCellState(i,j,k) != SOLID"<<std::endl;
//
//				if(grid.isValidFace<minU>(i, j, k) && grid.isSolidCell(i - 1, j, k))
//					grid.cellFaceVel<minU>(i, j, k) = std::max<double>(0.f, grid.cellFaceVel<minU>(i, j, k));
//
//				if(grid.isValidFace<maxU>(i, j, k) && grid.isSolidCell(i + 1, j, k))
//					grid.cellFaceVel<maxU>(i, j, k) = std::min<double>(0.f, grid.cellFaceVel<maxU>(i, j, k));
//
//				if(grid.isValidFace<minV>(i, j, k) && grid.isSolidCell(i, j - 1, k))
//					grid.cellFaceVel<minV>(i, j, k) = std::max<double>(0.f, grid.cellFaceVel<minV>(i, j, k));
//
//				if(grid.isValidFace<maxV>(i, j, k) && grid.isSolidCell(i, j + 1, k))
//					grid.cellFaceVel<maxV>(i, j, k) = std::min<double>(0.f, grid.cellFaceVel<maxV>(i, j, k));
//
//				if(grid.isValidFace<minW>(i, j, k) && grid.isSolidCell(i, j, k - 1))
//					grid.cellFaceVel<minW>(i, j, k) = std::max<double>(0.f, grid.cellFaceVel<minW>(i, j, k));
//
//				if(grid.isValidFace<maxW>(i, j, k) && grid.isSolidCell(i, j, k + 1))
//					grid.cellFaceVel<maxW>(i, j, k) = std::min<double>(0.f, grid.cellFaceVel<maxW>(i, j, k));
//			}
//}

void extrapolateBoundaryVelocities(MacGrid & grid){

    for (uint i = 0 ; i < grid.cellNum_i ; i ++)
        for (uint j = 0 ; j < grid.cellNum_j ; j ++)
            for (uint k = 0 ; k < grid.cellNum_k ; k++){
                bool currentCellIsValid = grid.cellIsInBounds(i,j,k) && grid.getCellState(i,j,k) != LIQUID;
                //U Faces
                if (currentCellIsValid || (grid.cellIsInBounds(i+1,j,k) && grid.getCellState(i+1,j,k) != LIQUID)) {

                    bool fromUp    = (grid.cellIsInBounds(i  ,j+1,k  ) && grid.getCellState(i  ,j+1,k  ) == LIQUID) || (grid.cellIsInBounds(i+1,j+1,k) && grid.getCellState(i+1,j+1,k  ) == LIQUID);
                    bool fromDown  = (grid.cellIsInBounds(i  ,j-1,k  ) && grid.getCellState(i  ,j-1,k  ) == LIQUID) || (grid.cellIsInBounds(i+1,j-1,k) && grid.getCellState(i+1,j-1,k  ) == LIQUID);
                    bool fromFront = (grid.cellIsInBounds(i  ,j  ,k+1) && grid.getCellState(i  ,j  ,k+1) == LIQUID) || (grid.cellIsInBounds(i+1,j,k+1) && grid.getCellState(i+1,j  ,k+1) == LIQUID);
                    bool fromBack  = (grid.cellIsInBounds(i  ,j  ,k-1) && grid.getCellState(i  ,j  ,k-1) == LIQUID) || (grid.cellIsInBounds(i+1,j,k-1) && grid.getCellState(i+1,j  ,k-1) == LIQUID);

                    double val = 0;
                    int count = fromUp + fromDown + fromFront + fromBack;

                    if (fromUp) val += grid.cellFaceVel<maxU>(i,j+1,k);
                    if (fromDown) val += grid.cellFaceVel<maxU>(i,j-1,k);
                    if (fromFront) val += grid.cellFaceVel<maxU>(i,j,k+1);
                    if (fromBack) val += grid.cellFaceVel<maxU>(i,j,k-1);
                
                    if (count > 0) {
                        grid.cellFaceVel<maxU>(i,j,k) = val / count;
                    }
                }
                //V Faces
                if (currentCellIsValid || (grid.cellIsInBounds(i,j+1,k) && grid.getCellState(i,j+1,k) != LIQUID)) {

                    bool fromUp    = (grid.cellIsInBounds(i+1,j,k) && grid.getCellState(i+1,j,k) == LIQUID) || (grid.cellIsInBounds(i+1,j+1,k) && grid.getCellState(i+1,j+1,k) == LIQUID);
                    bool fromDown  = (grid.cellIsInBounds(i-1,j,k) && grid.getCellState(i-1,j,k) == LIQUID) || (grid.cellIsInBounds(i-1,j+1,k) && grid.getCellState(i-1,j+1,k) == LIQUID);
                    bool fromFront = (grid.cellIsInBounds(i,j,k+1) && grid.getCellState(i,j,k+1) == LIQUID) || (grid.cellIsInBounds(i,j+1,k+1) && grid.getCellState(i,j+1,k+1) == LIQUID);
                    bool fromBack  = (grid.cellIsInBounds(i,j,k-1) && grid.getCellState(i,j,k-1) == LIQUID) || (grid.cellIsInBounds(i,j+1,k-1) && grid.getCellState(i,j+1,k-1) == LIQUID);
                
                    double val = 0;
                    int count = fromUp + fromDown + fromFront + fromBack;
                    if (fromUp) val += grid.cellFaceVel<maxV>(i+1,j,k);
                    if (fromDown) val += grid.cellFaceVel<maxV>(i-1,j,k);
                    if (fromFront) val += grid.cellFaceVel<maxV>(i,j,k+1);
                    if (fromBack) val += grid.cellFaceVel<maxV>(i,j,k-1);
                
                    if (count > 0) {
                        grid.cellFaceVel<maxV>(i,j,k) = val / count;
                    }
                }
                //W Faces
                if (currentCellIsValid || (grid.cellIsInBounds(i,j,k+1) && grid.getCellState(i,j,k+1) != LIQUID)) {

                    bool fromUp    = (grid.cellIsInBounds(i+1,j,k) && grid.getCellState(i+1,j,k) == LIQUID) || (grid.cellIsInBounds(i+1,j,k+1) && grid.getCellState(i+1,j,k+1) == LIQUID);
                    bool fromDown  = (grid.cellIsInBounds(i-1,j,k) && grid.getCellState(i-1,j,k) == LIQUID) || (grid.cellIsInBounds(i-1,j,k+1) && grid.getCellState(i-1,j,k+1) == LIQUID);
                    bool fromFront = (grid.cellIsInBounds(i,j+1,k) && grid.getCellState(i,j+1,k) == LIQUID) || (grid.cellIsInBounds(i,j+1,k+1) && grid.getCellState(i,j+1,k+1) == LIQUID);
                    bool fromBack  = (grid.cellIsInBounds(i,j-1,k) && grid.getCellState(i,j-1,k) == LIQUID) || (grid.cellIsInBounds(i,j-1,k+1) && grid.getCellState(i,j-1,k+1) == LIQUID);
                
                    double val = 0;
                    int count = fromUp + fromDown + fromFront + fromBack;
                    if (fromUp) val += grid.cellFaceVel<maxW>(i+1,j,k);
                    if (fromDown) val += grid.cellFaceVel<maxW>(i-1,j,k);
                    if (fromFront) val += grid.cellFaceVel<maxW>(i,j+1,k);
                    if (fromBack) val += grid.cellFaceVel<maxW>(i,j-1,k);
                
                    if (count > 0) {
                        grid.cellFaceVel<maxW>(i,j,k) = val / count;
                    }
                }
            }
}

void enforceBoundaryPressure(MacGrid & grid){
	uint jump_j = 0;
	uint jump_k = 0;
	for (uint i = 1 ; i < grid.cellNum_i - 1 ; i++ )
		for (uint j = 1 ; j < grid.cellNum_j - 1 ; j += jump_j)
			for (uint k = 1 ; k < grid.cellNum_k - 1 ; k += jump_k){

				if(2 < i && i < grid.cellNum_i - 2) jump_j = grid.cellNum_j - 2;
				else jump_j = 1;

				if(2 < j && j < grid.cellNum_j - 2) jump_k = grid.cellNum_k - 2;
				else jump_k = 1;

				//bool isSolid = grid.isSolidCell(i, j, k);
				if(grid.getCellState(i, j, k) != LIQUID)
					continue;

				uint neighbourNum = 0;
				auto idx = grid.cellCenterIdx(i, j, k);

				if(!grid.isSolidCell(i - 1, j, k))
					neighbourNum += 1;
				if(!grid.isSolidCell(i + 1, j, k))
					neighbourNum += 1;
				if(!grid.isSolidCell(i, j - 1, k))
					neighbourNum += 1;
				if(!grid.isSolidCell(i, j + 1, k))
					neighbourNum += 1;
				if(!grid.isSolidCell(i, j, k - 1))
					neighbourNum += 1;
				if(!grid.isSolidCell(i, j, k + 1))
					neighbourNum += 1;

				if(neighbourNum != 6){
					assert(neighbourNum >= 3);
					grid.cellCenterPressure.at(idx) *= 6.f/(double)neighbourNum;

				}
			}
}

void updateParticles(Particles & particles, double subStep){
	for(uint pi = 0 ; pi < particles.positions.size() ; pi++){
		auto newVel = particles.velocities[pi] * subStep;
		DEBUG_VAR(newVel);
		particles.positions[pi] += newVel;
	}
}

void collisionBasedParticleUpdate(Particles & particles, MacGrid & grid, double subStep){

	for(uint pi = 0 ; pi < particles.num ; pi++){
		auto pVel = particles.getVel(pi);
		auto pPos = particles.getPos(pi);

		//Vector3d planeNormal;
		//Vector3d planePos;

		Vector3d tempVel = pVel * subStep;
		Vector3d tempPos = pPos;

		//while(boxCollision(grid.origin, grid.end, tempPos, tempVel)){}
		if(boxCollision(grid.CBStart, grid.CBEnd, tempPos, tempVel)){
			pVel = tempVel / subStep;
		}
		else{
			pPos += tempVel;
		}
		assert(pPos[0] > grid.CBStart[0] && grid.CBEnd[0] > pPos[0]);
		assert(pPos[1] > grid.CBStart[1] && grid.CBEnd[1] > pPos[1]);
		assert(pPos[2] > grid.CBStart[2] && grid.CBEnd[2] > pPos[2]);
	}
}


double calculateSubStep(MacGrid const & grid, double timeStep){

	double maxVelComponent = grid.cellSize;

	int u,v,w;
	for (uint i = 0 ; i < grid.cellFaceVel_u.size(); i++){
		auto projected_u = std::fabs(grid.cellFaceVel_u[i]);
		if(projected_u > maxVelComponent){
			maxVelComponent = projected_u;
			u = 1; v = 0; w = 0;
		}
	}

	for (uint j = 0 ; j < grid.cellFaceVel_v.size(); j++){
		auto projected_v = std::fabs(grid.cellFaceVel_v[j]);
		if(projected_v > maxVelComponent){
			maxVelComponent = projected_v;
			u = 0; v = 1; w = 0;
		}
	}

	for (uint k = 0 ; k < grid.cellFaceVel_w.size(); k++){
		auto projected_w = std::fabs(grid.cellFaceVel_w[k]);
		if(projected_w > maxVelComponent){
			maxVelComponent = projected_w;
			u = 0; v = 0; w = 1;
		}
	}
	double subStep =  (grid.cellSize / maxVelComponent) * timeStep;
	//if(u)
	//	std::cout<<"maxVelComponent was in u dir: "<<maxVelComponent<<std::endl;
	//else if(v)
	//	std::cout<<"maxVelComponent was in v dir: "<<maxVelComponent<<std::endl;
	//else if(w)
	//	std::cout<<"maxVelComponent was in w dir: "<<maxVelComponent<<std::endl;

	return subStep;
}

//=============================================================================================
//=====================================PIC FUNCTIONS===========================================
//=============================================================================================

void transferAttributes(Particles const & particles, MacGrid & grid){

	std::fill(grid.cellCenterDensity.begin(), grid.cellCenterDensity.end(), 0);
	
	std::fill(grid.cellFaceVel_u.begin(), grid.cellFaceVel_u.end(), 0);
	std::fill(grid.cellFaceVel_v.begin(), grid.cellFaceVel_v.end(), 0);
	std::fill(grid.cellFaceVel_w.begin(), grid.cellFaceVel_w.end(), 0);

	std::fill(grid.cellFaceWeightSum_u.begin(), grid.cellFaceWeightSum_u.end(), 0);
	std::fill(grid.cellFaceWeightSum_v.begin(), grid.cellFaceWeightSum_v.end(), 0);
	std::fill(grid.cellFaceWeightSum_w.begin(), grid.cellFaceWeightSum_w.end(), 0);

	//Set inner bound cells to AIR
	for ( uint i=1 ; i<grid.cellNum_i-1 ; i++)
		for ( uint j=1 ; j<grid.cellNum_j-1 ; j++)
			for ( uint k=1 ; k<grid.cellNum_k-1 ; k++)
				grid.getCellState(i,j,k) = AIR;

	DEBUG();

	//pi is particle index
	for(uint pi = 0 ; pi < particles.num ; pi++){
		auto pPos = particles.getPos(pi);
		auto pVel = particles.getVel(pi);
		auto [i, j, k] = grid.gridCoord(pPos);
		auto localNBRIdcs = getLocalNBRIdcs(grid, pPos, i, j, k);
		auto [min_i, min_j, min_k, max_i, max_j, max_k] = localNBRIdcs;

		tuple8i cellUFaceNBRIdcs = getStageredCellFaceNBRIdcs_u(grid, i, min_j, min_k, max_j, max_k);
		auto cellFacePosc000_u = grid.cellFacePos<minU>(i, min_j, min_k);
		auto cellFacePosc111_u = grid.cellFacePos<maxU>(i, max_j, max_k);

		auto weights_u = getWeights(pPos, cellFacePosc000_u, cellFacePosc111_u);
		//auto weights_u = getDistWeights(pPos, cellFacePosc000_u, cellFacePosc111_u);

		for(uint ii=0 ; ii<8 ; ii++){
			auto fi_u = cellUFaceNBRIdcs[ii];
			grid.cellFaceVel_u.at(fi_u) += weights_u[ii] * pVel[0];
			grid.cellFaceWeightSum_u.at(fi_u) += weights_u.at(ii);
		}

		tuple8i cellVFaceNBRIdcs = getStageredCellFaceNBRIdcs_v(grid, j, min_i, min_k, max_i, max_k);
		auto cellFacePosc000_v = grid.cellFacePos<minV>(min_i, j, min_k);
		auto cellFacePosc111_v = grid.cellFacePos<maxV>(max_i, j, max_k);

		auto weights_v = getWeights(pPos, cellFacePosc000_v, cellFacePosc111_v);
		//auto weights_v = getDistWeights(pPos, cellFacePosc000_v, cellFacePosc111_v);

		for(uint ii=0 ; ii<8 ; ii++){
			auto fi_v = cellVFaceNBRIdcs[ii];
			grid.cellFaceVel_v.at(fi_v) += weights_v[ii] * pVel[1];
			grid.cellFaceWeightSum_v.at(fi_v) += weights_v[ii];
		}

		tuple8i cellWFaceNBRIdcs = getStageredCellFaceNBRIdcs_w(grid, k, min_i, min_j, max_i, max_j);
		auto cellFacePosc000_w = grid.cellFacePos<minW>(min_i, min_j, k);
		auto cellFacePosc111_w = grid.cellFacePos<maxW>(max_i, max_j, k);

		auto weights_w = getWeights(pPos, cellFacePosc000_w, cellFacePosc111_w);
		//auto weights_w = getDistWeights(pPos, cellFacePosc000_w, cellFacePosc111_w);

		for(uint ii=0 ; ii<8 ; ii++){
			auto fi_w = cellWFaceNBRIdcs[ii];
			grid.cellFaceVel_w.at(fi_w) += weights_w[ii] * pVel[2];
			grid.cellFaceWeightSum_w.at(fi_w) += weights_w[ii];
		}

		//Set Cell to liquid
		auto idx = grid.cellCenterIdx(i,j,k);
		if(ENFORCE_INCOMPRESSIBILITY)
			grid.cellCenterDensity.at(idx) += 1;
		grid.cellCenterState.at(idx) = LIQUID;
		//grid.getCellState(i,j,k) = LIQUID;
	}

	DEBUG();

	for (uint i_u = 0; i_u < grid.cellFaceVel_u.size() ; i_u++)
		if(grid.cellFaceWeightSum_u[i_u])
			grid.cellFaceVel_u[i_u] /= grid.cellFaceWeightSum_u[i_u];
		else
			grid.cellFaceVel_u[i_u] = 0;

	for (uint i_v = 0; i_v < grid.cellFaceVel_v.size() ; i_v++)
		if(grid.cellFaceWeightSum_v[i_v])
			grid.cellFaceVel_v[i_v] /= grid.cellFaceWeightSum_v[i_v];
		else
			grid.cellFaceVel_v[i_v] = 0;


	for (uint i_w = 0; i_w < grid.cellFaceVel_w.size() ; i_w++)
		if(grid.cellFaceWeightSum_w[i_w])
			grid.cellFaceVel_w[i_w] /= grid.cellFaceWeightSum_w[i_w];
		else
			grid.cellFaceVel_w[i_w] = 0;

	setDefaultCellStates(grid);
		
}

void transferAttributes(MacGrid const & grid, Particles & particles){

	std::cout<<"Pic"<<std::endl;
	//pi is particle index
	for(uint pi = 0 ; pi < particles.num ; pi++){
		Vector3dRef pPos = particles.getPos(pi);
		Vector3dRef pVel = particles.getVel(pi);
		auto [i, j, k] = grid.gridCoord(pPos);
		auto [min_i, min_j, min_k, max_i, max_j, max_k] = getLocalNBRIdcs(grid, pPos, i, j, k);

		tuple8i cellUFaceNBRIdcs = getStageredCellFaceNBRIdcs_u(grid, i, min_j, min_k, max_j, max_k);
		auto cellFacePosc000_u = grid.cellFacePos<minU>(i, min_j, min_k);
		auto cellFacePosc111_u = grid.cellFacePos<maxU>(i, max_j, max_k);
		auto cellFaceNBRVels_u = getFromIdcs(grid.cellFaceVel_u, cellUFaceNBRIdcs);
		pVel(0) = trilinearInterpolation(cellFaceNBRVels_u, getDiff(pPos, cellFacePosc000_u, cellFacePosc111_u));

		tuple8i cellVFaceNBRIdcs = getStageredCellFaceNBRIdcs_v(grid, j, min_i, min_k, max_i, max_k);
		auto cellFacePosc000_v = grid.cellFacePos<minV>(min_i, j, min_k);
		auto cellFacePosc111_v = grid.cellFacePos<maxV>(max_i, j, max_k);
		auto cellFaceNBRVels_v = getFromIdcs(grid.cellFaceVel_v, cellVFaceNBRIdcs);
		pVel(1) = trilinearInterpolation(cellFaceNBRVels_v, getDiff(pPos, cellFacePosc000_v, cellFacePosc111_v));

		tuple8i cellWFaceNBRIdcs = getStageredCellFaceNBRIdcs_w(grid, k, min_i, min_j, max_i, max_j);
		auto cellFacePosc000_w = grid.cellFacePos<minW>(min_i, min_j, k);
		auto cellFacePosc111_w = grid.cellFacePos<maxW>(max_i, max_j, k);
		auto cellFaceNBRVels_w = getFromIdcs(grid.cellFaceVel_w, cellWFaceNBRIdcs);
		pVel(2) = trilinearInterpolation(cellFaceNBRVels_w, getDiff(pPos, cellFacePosc000_w, cellFacePosc111_w));
	}
}


void advanceStep(Particles & particles, MacGrid & grid, double timeStep){
	DEBUG();
	double step = 0;
	int iterations = 0;
	while((step < timeStep)){//&& (iterations < 5)){

		//std::cout<<"===================================STEP: "<<
		//step<<"/"<<timeStep<<"==================================="<<std::endl;

		transferAttributes(particles, grid);

			//std::cout<<"applyExternalForces(grid, timeStep)"<<std::endl;
			applyExternalForces(grid, timeStep);
			//std::cout<<"enforceBoundaryVelocities(grid)"<<std::endl;
			//grid.printOutMaxFaceVels();
			enforceBoundaryVelocities(grid);
			//std::cout<<"applyPressureForces(grid, timeStep)"<<std::endl;
			applyPressureForces(grid, timeStep);
			//std::cout<<"enforceBoundaryVelocities(grid)"<<std::endl;
			//grid.printOutMaxFaceVels();
			enforceBoundaryVelocities(grid);
			//std::cout<<"extrapolateBoundaryVelocities(grid"<<std::endl;
			extrapolateBoundaryVelocities(grid);
			//std::cout<<"transferAttributes(grid, particles)"<<std::endl;
			enforceBoundaryVelocities(grid);

		transferAttributes(grid, particles);

		//double subStep = calculateSubStep(grid, timeStep);
		double subStep = timeStep;

		collisionBasedParticleUpdate(particles, grid, subStep);
		//updateParticles(particles, timeStep);

		step += subStep;
		iterations++;
	}
	//std::cout<<"advanceStep done!"<<std::endl;
}

//=============================================================================================
//=====================================FLIP FUNCTIONS==========================================
//=============================================================================================

void transferAttributes(Particles const & particles, FlipMacGrid & grid){

	grid.cellFaceVelOld_u.swap(grid.cellFaceVel_u);
	grid.cellFaceVelOld_v.swap(grid.cellFaceVel_v);
	grid.cellFaceVelOld_w.swap(grid.cellFaceVel_w);

	std::fill(grid.cellCenterDensity.begin(), grid.cellCenterDensity.end(), 0);

	std::fill(grid.cellFaceVel_u.begin(), grid.cellFaceVel_u.end(), 0);
	std::fill(grid.cellFaceVel_v.begin(), grid.cellFaceVel_v.end(), 0);
	std::fill(grid.cellFaceVel_w.begin(), grid.cellFaceVel_w.end(), 0);

	std::fill(grid.cellFaceWeightSum_u.begin(), grid.cellFaceWeightSum_u.end(), 0);
	std::fill(grid.cellFaceWeightSum_v.begin(), grid.cellFaceWeightSum_v.end(), 0);
	std::fill(grid.cellFaceWeightSum_w.begin(), grid.cellFaceWeightSum_w.end(), 0);

	//Set inner bound cells to AIR
	for ( uint i=1 ; i<grid.cellNum_i-1 ; i++)
		for ( uint j=1 ; j<grid.cellNum_j-1 ; j++)
			for ( uint k=1 ; k<grid.cellNum_k-1 ; k++)
				grid.getCellState(i,j,k) = AIR;

	DEBUG();

	//pi is particle index
	for(uint pi = 0 ; pi < particles.num ; pi++){
		auto pPos = particles.getPos(pi);
		auto pVel = particles.getVel(pi);
		auto [i, j, k] = grid.gridCoord(pPos);
		auto localNBRIdcs = getLocalNBRIdcs(grid, pPos, i, j, k);
		auto [min_i, min_j, min_k, max_i, max_j, max_k] = localNBRIdcs;

		{
			tuple8i cellUFaceNBRIdcs = getStageredCellFaceNBRIdcs_u(grid, i, min_j, min_k, max_j, max_k);
			auto cellFacePosc000_u = grid.cellFacePos<minU>(i, min_j, min_k);
			auto cellFacePosc111_u = grid.cellFacePos<maxU>(i, max_j, max_k);
	
			auto weights_u = getWeights(pPos, cellFacePosc000_u, cellFacePosc111_u);
			//auto weights_u = getDistWeights(pPos, cellFacePosc000_u, cellFacePosc111_u);
	
			for(uint ii=0 ; ii<8 ; ii++){
				auto fi_u = cellUFaceNBRIdcs[ii];
				grid.cellFaceVel_u.at(fi_u) += weights_u[ii] * pVel[0];
				grid.cellFaceWeightSum_u.at(fi_u) += weights_u.at(ii);
			}
		}

		tuple8i cellVFaceNBRIdcs = getStageredCellFaceNBRIdcs_v(grid, j, min_i, min_k, max_i, max_k);
		auto cellFacePosc000_v = grid.cellFacePos<minV>(min_i, j, min_k);
		auto cellFacePosc111_v = grid.cellFacePos<maxV>(max_i, j, max_k);

		auto weights_v = getWeights(pPos, cellFacePosc000_v, cellFacePosc111_v);
		//auto weights_v = getDistWeights(pPos, cellFacePosc000_v, cellFacePosc111_v);

		for(uint ii=0 ; ii<8 ; ii++){
			auto fi_v = cellVFaceNBRIdcs[ii];
			grid.cellFaceVel_v.at(fi_v) += weights_v[ii] * pVel[1];
			grid.cellFaceWeightSum_v.at(fi_v) += weights_v[ii];
		}

		tuple8i cellWFaceNBRIdcs = getStageredCellFaceNBRIdcs_w(grid, k, min_i, min_j, max_i, max_j);
		auto cellFacePosc000_w = grid.cellFacePos<minW>(min_i, min_j, k);
		auto cellFacePosc111_w = grid.cellFacePos<maxW>(max_i, max_j, k);

		auto weights_w = getWeights(pPos, cellFacePosc000_w, cellFacePosc111_w);
		//auto weights_w = getDistWeights(pPos, cellFacePosc000_w, cellFacePosc111_w);

		for(uint ii=0 ; ii<8 ; ii++){
			auto fi_w = cellWFaceNBRIdcs[ii];
			grid.cellFaceVel_w.at(fi_w) += weights_w[ii] * pVel[2];
			grid.cellFaceWeightSum_w.at(fi_w) += weights_w[ii];
		}

		//Set Cell to liquid
		auto idx = grid.cellCenterIdx(i,j,k);
		if(ENFORCE_INCOMPRESSIBILITY)
			grid.cellCenterDensity.at(idx) += 1;
		grid.cellCenterState.at(idx) = LIQUID;
		//grid.getCellState(i,j,k) = LIQUID;
	}

	DEBUG();

	for (uint i_u = 0; i_u < grid.cellFaceVel_u.size() ; i_u++)
		if(grid.cellFaceWeightSum_u[i_u])
			grid.cellFaceVel_u[i_u] /= grid.cellFaceWeightSum_u[i_u];
		else
			grid.cellFaceVel_u[i_u] = 0;

	for (uint i_v = 0; i_v < grid.cellFaceVel_v.size() ; i_v++)
		if(grid.cellFaceWeightSum_v[i_v])
			grid.cellFaceVel_v[i_v] /= grid.cellFaceWeightSum_v[i_v];
		else
			grid.cellFaceVel_v[i_v] = 0;


	for (uint i_w = 0; i_w < grid.cellFaceVel_w.size() ; i_w++)
		if(grid.cellFaceWeightSum_w[i_w])
			grid.cellFaceVel_w[i_w] /= grid.cellFaceWeightSum_w[i_w];
		else
			grid.cellFaceVel_w[i_w] = 0;

	setDefaultCellStates(grid);
		
}

void transferAttributes(FlipMacGrid const & grid, Particles & particles){

	std::cout<<"Flip"<<std::endl;
	//pi is particle index
	for(uint pi = 0 ; pi < particles.num ; pi++){
		Vector3dRef pPos = particles.getPos(pi);
		Vector3dRef pVel = particles.getVel(pi);
		Vector3d vel;
		Vector3d velOld;
		auto [i, j, k] = grid.gridCoord(pPos);
		auto [min_i, min_j, min_k, max_i, max_j, max_k] = getLocalNBRIdcs(grid, pPos, i, j, k);

		tuple8i cellUFaceNBRIdcs = getStageredCellFaceNBRIdcs_u(grid, i, min_j, min_k, max_j, max_k);
		auto cellFacePosc000_u = grid.cellFacePos<minU>(i, min_j, min_k);
		auto cellFacePosc111_u = grid.cellFacePos<maxU>(i, max_j, max_k);
		auto cellFaceNBRVels_u = getFromIdcs(grid.cellFaceVel_u, cellUFaceNBRIdcs);
		auto cellFaceNBRVelsOld_u = getFromIdcs(grid.cellFaceVelOld_u, cellUFaceNBRIdcs);
		vel[0] = trilinearInterpolation(cellFaceNBRVels_u, getDiff(pPos, cellFacePosc000_u, cellFacePosc111_u));
		velOld[0] = trilinearInterpolation(cellFaceNBRVelsOld_u, getDiff(pPos, cellFacePosc000_u, cellFacePosc111_u));

		tuple8i cellVFaceNBRIdcs = getStageredCellFaceNBRIdcs_v(grid, j, min_i, min_k, max_i, max_k);
		auto cellFacePosc000_v = grid.cellFacePos<minV>(min_i, j, min_k);
		auto cellFacePosc111_v = grid.cellFacePos<maxV>(max_i, j, max_k);
		auto cellFaceNBRVels_v = getFromIdcs(grid.cellFaceVel_v, cellVFaceNBRIdcs);
		auto cellFaceNBRVelsOld_v = getFromIdcs(grid.cellFaceVelOld_v, cellVFaceNBRIdcs);
		vel[1] = trilinearInterpolation(cellFaceNBRVels_v, getDiff(pPos, cellFacePosc000_v, cellFacePosc111_v));
		velOld[1] = trilinearInterpolation(cellFaceNBRVelsOld_v, getDiff(pPos, cellFacePosc000_v, cellFacePosc111_v));

		tuple8i cellWFaceNBRIdcs = getStageredCellFaceNBRIdcs_w(grid, k, min_i, min_j, max_i, max_j);
		auto cellFacePosc000_w = grid.cellFacePos<minW>(min_i, min_j, k);
		auto cellFacePosc111_w = grid.cellFacePos<maxW>(max_i, max_j, k);
		auto cellFaceNBRVels_w = getFromIdcs(grid.cellFaceVel_w, cellWFaceNBRIdcs);
		auto cellFaceNBRVelsOld_w = getFromIdcs(grid.cellFaceVelOld_w, cellWFaceNBRIdcs);
		vel[2] = trilinearInterpolation(cellFaceNBRVels_w, getDiff(pPos, cellFacePosc000_w, cellFacePosc111_w));
		velOld[2] = trilinearInterpolation(cellFaceNBRVelsOld_w, getDiff(pPos, cellFacePosc000_w, cellFacePosc111_w));

		//     Pic Equation       		  Flip Equation
		//pVel = vel * (1.f - FPIC_RATIO) + (pVel + (vel - velOld)) * FPIC_RATIO;
		pVel = vel + (pVel + (vel - velOld)) * FPIC_RATIO;
	}
}


void advanceStep(Particles & particles, FlipMacGrid & grid, double timeStep){
	DEBUG();
	double step = 0;
	int iterations = 0;
	while((step < timeStep)){//&& (iterations < 5)){

		//std::cout<<"===================================STEP: "<<
		//step<<"/"<<timeStep<<"==================================="<<std::endl;

		transferAttributes(particles, grid);

			//std::cout<<"applyExternalForces(grid, timeStep)"<<std::endl;
			applyExternalForces(grid, timeStep);
			//std::cout<<"enforceBoundaryVelocities(grid)"<<std::endl;
			//grid.printOutMaxFaceVels();
			enforceBoundaryVelocities(grid);
			//std::cout<<"applyPressureForces(grid, timeStep)"<<std::endl;
			applyPressureForces(grid, timeStep);
			//std::cout<<"enforceBoundaryVelocities(grid)"<<std::endl;
			//grid.printOutMaxFaceVels();
			enforceBoundaryVelocities(grid);
			//std::cout<<"extrapolateBoundaryVelocities(grid"<<std::endl;
			extrapolateBoundaryVelocities(grid);
			//std::cout<<"transferAttributes(grid, particles)"<<std::endl;
			enforceBoundaryVelocities(grid);

		transferAttributes(grid, particles);

		//double subStep = calculateSubStep(grid, timeStep);
		double subStep = timeStep;

		collisionBasedParticleUpdate(particles, grid, subStep);
		//updateParticles(particles, timeStep);

		step += subStep;
		//iterations++;
	}
	//std::cout<<"advanceStep done!"<<std::endl;
}


//=============================================================================================
//=====================================APIC FUNCTIONS==========================================
//=============================================================================================


void transferAttributes(AffineParticles const & particles, MacGrid & grid){

	std::fill(grid.cellCenterDensity.begin(), grid.cellCenterDensity.end(), 0);

	std::fill(grid.cellFaceVel_u.begin(), grid.cellFaceVel_u.end(), 0);
	std::fill(grid.cellFaceVel_v.begin(), grid.cellFaceVel_v.end(), 0);
	std::fill(grid.cellFaceVel_w.begin(), grid.cellFaceVel_w.end(), 0);

	std::fill(grid.cellFaceWeightSum_u.begin(), grid.cellFaceWeightSum_u.end(), 0);
	std::fill(grid.cellFaceWeightSum_v.begin(), grid.cellFaceWeightSum_v.end(), 0);
	std::fill(grid.cellFaceWeightSum_w.begin(), grid.cellFaceWeightSum_w.end(), 0);

	//Set inner bound cells to AIR
	for ( uint i=1 ; i<grid.cellNum_i-1 ; i++)
		for ( uint j=1 ; j<grid.cellNum_j-1 ; j++)
			for ( uint k=1 ; k<grid.cellNum_k-1 ; k++)
				grid.getCellState(i,j,k) = AIR;

	DEBUG();

	//pi is particle index
	for(uint pi = 0 ; pi < particles.num ; pi++){
		auto pPos = particles.getPos(pi);
		auto pVel = particles.getVel(pi);
		auto [i, j, k] = grid.gridCoord(pPos);
		auto localNBRIdcs = getLocalNBRIdcs(grid, pPos, i, j, k);
		auto [min_i, min_j, min_k, max_i, max_j, max_k] = localNBRIdcs;

		auto pAffine_x = particles.getAffine<X>(pi);
		tuple8i cellUFaceNBRIdcs = getStageredCellFaceNBRIdcs_u(grid, i, min_j, min_k, max_j, max_k);
		std::array<Vector3d, 8> cellFaceNBRPoses_u = getStageredCellFacePos_u(grid, i, min_j, min_k, max_j, max_k);
		auto cellFacePosc000_u = cellFaceNBRPoses_u[0];
		auto cellFacePosc111_u = cellFaceNBRPoses_u[7];

		auto weights_u = getWeights(pPos, cellFacePosc000_u, cellFacePosc111_u);
		//auto weights_u = getDistWeights(pPos, cellFacePosc000_u, cellFacePosc111_u);
		for(uint ii=0 ; ii<8 ; ii++){
			auto fi_u = cellUFaceNBRIdcs[ii];
			grid.cellFaceVel_u.at(fi_u) += weights_u[ii] * (pVel[0] + pAffine_x.dot(cellFaceNBRPoses_u[ii] - pPos) * APIC_RATIO);
			grid.cellFaceWeightSum_u.at(fi_u) += weights_u.at(ii);
		}

		auto pAffine_y = particles.getAffine<Y>(pi);
		tuple8i cellVFaceNBRIdcs = getStageredCellFaceNBRIdcs_v(grid, j, min_i, min_k, max_i, max_k);
		std::array<Vector3d, 8> cellFaceNBRPoses_v = getStageredCellFacePos_v(grid, j, min_i, min_k, max_i, max_k);
		auto cellFacePosc000_v = cellFaceNBRPoses_v[0];
		auto cellFacePosc111_v = cellFaceNBRPoses_v[7];

		auto weights_v = getWeights(pPos, cellFacePosc000_v, cellFacePosc111_v);
		//auto weights_v = getDistWeights(pPos, cellFacePosc000_v, cellFacePosc111_v);

		for(uint ii=0 ; ii<8 ; ii++){
			auto fi_v = cellVFaceNBRIdcs[ii];
			grid.cellFaceVel_v.at(fi_v) += weights_v[ii] * (pVel[1] + pAffine_y.dot(cellFaceNBRPoses_v[ii] - pPos) * APIC_RATIO);
			grid.cellFaceWeightSum_v.at(fi_v) += weights_v[ii];
		}

		auto pAffine_z = particles.getAffine<Z>(pi);
		tuple8i cellWFaceNBRIdcs = getStageredCellFaceNBRIdcs_w(grid, k, min_i, min_j, max_i, max_j);
		std::array<Vector3d, 8> cellFaceNBRPoses_w = getStageredCellFacePos_w(grid, k, min_i, min_j, max_i, max_j);
		auto cellFacePosc000_w = grid.cellFacePos<minW>(min_i, min_j, k);
		auto cellFacePosc111_w = grid.cellFacePos<maxW>(max_i, max_j, k);

		auto weights_w = getWeights(pPos, cellFacePosc000_w, cellFacePosc111_w);
		//auto weights_w = getDistWeights(pPos, cellFacePosc000_w, cellFacePosc111_w);

		for(uint ii=0 ; ii<8 ; ii++){
			auto fi_w = cellWFaceNBRIdcs[ii];
			grid.cellFaceVel_w.at(fi_w) += weights_w[ii] * (pVel[2] + pAffine_z.dot(cellFaceNBRPoses_w[ii] - pPos) * APIC_RATIO);
			grid.cellFaceWeightSum_w.at(fi_w) += weights_w[ii];
		}

		//Set Cell to liquid
		auto idx = grid.cellCenterIdx(i,j,k);
		if(ENFORCE_INCOMPRESSIBILITY)
			grid.cellCenterDensity.at(idx) += 1;
		grid.cellCenterState.at(idx) = LIQUID;
		//grid.getCellState(i,j,k) = LIQUID;
	}

	DEBUG();

	for (uint i_u = 0; i_u < grid.cellFaceVel_u.size() ; i_u++)
		if(grid.cellFaceWeightSum_u[i_u])
			grid.cellFaceVel_u[i_u] /= grid.cellFaceWeightSum_u[i_u];
		else
			grid.cellFaceVel_u[i_u] = 0;

	for (uint i_v = 0; i_v < grid.cellFaceVel_v.size() ; i_v++)
		if(grid.cellFaceWeightSum_v[i_v])
			grid.cellFaceVel_v[i_v] /= grid.cellFaceWeightSum_v[i_v];
		else
			grid.cellFaceVel_v[i_v] = 0;


	for (uint i_w = 0; i_w < grid.cellFaceVel_w.size() ; i_w++)
		if(grid.cellFaceWeightSum_w[i_w])
			grid.cellFaceVel_w[i_w] /= grid.cellFaceWeightSum_w[i_w];
		else
			grid.cellFaceVel_w[i_w] = 0;

	setDefaultCellStates(grid);
		
}

void transferAttributes(MacGrid const & grid, AffineParticles & particles){

	std::cout<<"Apic"<<std::endl;
	//pi is particle index
	for(uint pi = 0 ; pi < particles.num ; pi++){
		Vector3dRef pPos = particles.getPos(pi);
		Vector3dRef pVel = particles.getVel(pi);

		auto [i, j, k] = grid.gridCoord(pPos);
		auto [min_i, min_j, min_k, max_i, max_j, max_k] = getLocalNBRIdcs(grid, pPos, i, j, k);

		auto pAffine_x = particles.getAffine<X>(pi);
		tuple8i cellUFaceNBRIdcs = getStageredCellFaceNBRIdcs_u(grid, i, min_j, min_k, max_j, max_k);
		auto cellFacePosc000_u = grid.cellFacePos<minU>(i, min_j, min_k);
		auto cellFacePosc111_u = grid.cellFacePos<maxU>(i, max_j, max_k);
		auto cellFaceNBRVels_u = getFromIdcs(grid.cellFaceVel_u, cellUFaceNBRIdcs);
		pVel[0] = trilinearInterpolation(cellFaceNBRVels_u, getDiff(pPos, cellFacePosc000_u, cellFacePosc111_u));
		pAffine_x = gradientTrilinearInterpolation(cellFaceNBRVels_u, getDiff(pPos, cellFacePosc000_u, cellFacePosc111_u));// * grid.invCellSize;
		//std::cout<<pAffine_x<<std::endl;

		auto pAffine_y = particles.getAffine<Y>(pi);
		tuple8i cellVFaceNBRIdcs = getStageredCellFaceNBRIdcs_v(grid, j, min_i, min_k, max_i, max_k);
		auto cellFacePosc000_v = grid.cellFacePos<minV>(min_i, j, min_k);
		auto cellFacePosc111_v = grid.cellFacePos<maxV>(max_i, j, max_k);
		auto cellFaceNBRVels_v = getFromIdcs(grid.cellFaceVel_v, cellVFaceNBRIdcs);
		pVel[1] = trilinearInterpolation(cellFaceNBRVels_v, getDiff(pPos, cellFacePosc000_v, cellFacePosc111_v));
		pAffine_y = gradientTrilinearInterpolation(cellFaceNBRVels_v, getDiff(pPos, cellFacePosc000_v, cellFacePosc111_v));// * grid.invCellSize;
		//std::cout<<pAffine_y<<std::endl;

		auto pAffine_z = particles.getAffine<Z>(pi);
		tuple8i cellWFaceNBRIdcs = getStageredCellFaceNBRIdcs_w(grid, k, min_i, min_j, max_i, max_j);
		auto cellFacePosc000_w = grid.cellFacePos<minW>(min_i, min_j, k);
		auto cellFacePosc111_w = grid.cellFacePos<maxW>(max_i, max_j, k);
		auto cellFaceNBRVels_w = getFromIdcs(grid.cellFaceVel_w, cellWFaceNBRIdcs);
		pVel[2] = trilinearInterpolation(cellFaceNBRVels_w, getDiff(pPos, cellFacePosc000_w, cellFacePosc111_w));
		pAffine_z = gradientTrilinearInterpolation(cellFaceNBRVels_w, getDiff(pPos, cellFacePosc000_w, cellFacePosc111_w));// * grid.invCellSize;
		//std::cout<<pAffine_z<<std::endl;
	}
}


void advanceStep(AffineParticles & particles, MacGrid & grid, double timeStep){
	DEBUG();
	double step = 0;
	int iterations = 0;
	while((step < timeStep)){//&& (iterations < 5)){

		//std::cout<<"===================================STEP: "<<
		//step<<"/"<<timeStep<<"==================================="<<std::endl;

		transferAttributes(particles, grid);

			//std::cout<<"applyExternalForces(grid, timeStep)"<<std::endl;
			applyExternalForces(grid, timeStep);
			//std::cout<<"enforceBoundaryVelocities(grid)"<<std::endl;
			//grid.printOutMaxFaceVels();
			enforceBoundaryVelocities(grid);
			//std::cout<<"applyPressureForces(grid, timeStep)"<<std::endl;
			applyPressureForces(grid, timeStep);
			//std::cout<<"enforceBoundaryVelocities(grid)"<<std::endl;
			//grid.printOutMaxFaceVels();
			enforceBoundaryVelocities(grid);
			//std::cout<<"extrapolateBoundaryVelocities(grid"<<std::endl;
			extrapolateBoundaryVelocities(grid);
			//std::cout<<"transferAttributes(grid, particles)"<<std::endl;
			enforceBoundaryVelocities(grid);

		transferAttributes(grid, particles);

		//double subStep = calculateSubStep(grid, timeStep);
		double subStep = timeStep;

		collisionBasedParticleUpdate(particles, grid, subStep);
		//updateParticles(particles, timeStep);

		step += subStep;
		//iterations++;
	}
	//std::cout<<"advanceStep done!"<<std::endl;
}

/**

void extrapolateBoundaryVelocities(MacGrid & grid){

	auto extrapolateBoundaryCenterVel_i = [&grid](uint i, uint j, uint k){
		if(grid.isMaxBoundaryCell_i(i) && grid.isValidFace<maxV>(i, j, k))
			grid.cellFaceVel<maxV>(i, j, k) = grid.cellFaceVel<maxV>(i - 1, j, k);
		else if(grid.isMinBoundaryCell_i(i) && grid.isValidFace<maxV>(i, j, k))
			grid.cellFaceVel<maxV>(i, j, k) = grid.cellFaceVel<maxV>(i + 1, j, k);

		if(grid.isMaxBoundaryCell_i(i) && grid.isValidFace<maxW>(i, j, k))
			grid.cellFaceVel<maxW>(i, j, k) = grid.cellFaceVel<maxW>(i - 1, j, k);
		else if(grid.isMinBoundaryCell_i(i) && grid.isValidFace<maxW>(i, j, k))
			grid.cellFaceVel<maxW>(i, j, k) = grid.cellFaceVel<maxW>(i + 1, j, k);
	};

	auto extrapolateBoundaryCenterVel_j = [&grid](uint i, uint j, uint k){
		if(grid.isMaxBoundaryCell_j(j) && grid.isValidFace<maxU>(i, j, k))
			grid.cellFaceVel<maxU>(i, j, k) = grid.cellFaceVel<maxU>(i, j - 1, k);
		else if(grid.isMinBoundaryCell_j(j) && grid.isValidFace<maxU>(i, j, k))
			grid.cellFaceVel<maxU>(i, j, k) = grid.cellFaceVel<maxU>(i, j + 1, k);

		if(grid.isMaxBoundaryCell_j(j) && grid.isValidFace<maxW>(i, j, k))
			grid.cellFaceVel<maxW>(i, j, k) = grid.cellFaceVel<maxW>(i, j - 1, k);
		else if(grid.isMinBoundaryCell_j(j) && grid.isValidFace<maxW>(i, j, k))
			grid.cellFaceVel<maxW>(i, j, k) = grid.cellFaceVel<maxW>(i, j + 1, k);
	};

	auto extrapolateBoundaryCenterVel_k = [&grid](uint i, uint j, uint k){
		if(grid.isMaxBoundaryCell_k(k) && grid.isValidFace<maxU>(i, j, k))
			grid.cellFaceVel<maxU>(i, j, k) = grid.cellFaceVel<maxU>(i, j, k - 1);
		else if(grid.isMinBoundaryCell_k(k) && grid.isValidFace<maxU>(i, j, k))
			grid.cellFaceVel<maxU>(i, j, k) = grid.cellFaceVel<maxU>(i, j, k + 1);

		if(grid.isMaxBoundaryCell_k(k) && grid.isValidFace<maxV>(i, j, k))
			grid.cellFaceVel<maxV>(i, j, k) = grid.cellFaceVel<maxV>(i, j, k - 1);
		else if(grid.isMinBoundaryCell_k(k) && grid.isValidFace<maxV>(i, j, k))
			grid.cellFaceVel<maxV>(i, j, k) = grid.cellFaceVel<maxV>(i, j, k + 1);
	};

	for (uint i = 0 ; i < grid.cellNum_i ; i += (grid.cellNum_i - 1))
		for (uint j = 1 ; j < grid.cellNum_j - 1 ; j ++)
			for (uint k = 1 ; k < grid.cellNum_k - 1 ; k++)
				extrapolateBoundaryCenterVel_i(i,j,k);
	for (uint i = 1 ; i < grid.cellNum_i - 1 ; i ++)
		for (uint j = 0 ; j < grid.cellNum_j ; j += (grid.cellNum_j - 1))
			for (uint k = 1 ; k < grid.cellNum_k - 1 ; k ++)
				extrapolateBoundaryCenterVel_j(i,j,k);
	for (uint i = 1 ; i < grid.cellNum_i - 1 ; i ++)
		for (uint j = 1 ; j < grid.cellNum_j - 1 ; j++)
			for (uint k = 0 ; k < grid.cellNum_k ; k += (grid.cellNum_k - 1))
				extrapolateBoundaryCenterVel_k(i,j,k);


	for (uint i = 1 ; i < grid.cellNum_i - 1 ; i ++)
		for (uint j = 0 ; j < grid.cellNum_j ; j += (grid.cellNum_j - 1))
			for (uint k = 0 ; k < grid.cellNum_k ; k += (grid.cellNum_k - 1))
				extrapolateBoundaryCenterVel_i(i,j,k);
	for (uint i = 0 ; i < grid.cellNum_i ; i += (grid.cellNum_i - 1))
		for (uint j = 1 ; j < grid.cellNum_j - 1 ; j++)
			for (uint k = 0 ; k < grid.cellNum_k ; k += (grid.cellNum_k - 1))
				extrapolateBoundaryCenterVel_j(i,j,k);
	for (uint i = 0 ; i < grid.cellNum_i ; i += (grid.cellNum_i - 1))
		for (uint j = 0 ; j < grid.cellNum_j ; j += (grid.cellNum_j - 1))
			for (uint k = 1 ; k < grid.cellNum_k - 1 ; k++)
				extrapolateBoundaryCenterVel_k(i,j,k);

	for (uint i = 0 ; i < grid.cellNum_i ; i += (grid.cellNum_i - 1))
		for (uint j = 0 ; j < grid.cellNum_j ; j += (grid.cellNum_j - 1))
			for (uint k = 0 ; k < grid.cellNum_k ; k += (grid.cellNum_k - 1)){
				extrapolateBoundaryCenterVel_i(i,j,k);
				extrapolateBoundaryCenterVel_j(i,j,k);
				extrapolateBoundaryCenterVel_k(i,j,k);
			}
}

bool collides(MacGrid const & grid, Vector3d const & prev, Vector3d const & next, Vector3d & normal){
	//x collisions
	if(	grid.boundMax_y > prev[1] && prev[1] >= grid.boundMin_y && 
		grid.boundMax_z > prev[2] && prev[2] >= grid.boundMin_z){

		if(prev[0] > grid.boundMin_x && next[0] < grid.boundMin_x){
			normal = Vector3d(1.f, 0.f, 0.f);
			return true;
		}
		else if(prev[0] < grid.boundMin_x && next[0] > grid.boundMin_x){
			normal = Vector3d(-1.f, 0.f, 0.f);
			return true;
		}
		if(prev[0] > grid.boundMax_x && next[0] < grid.boundMax_x){
			normal = Vector3d(1.f, 0.f, 0.f);
			return true;
		}
		else if(prev[0] < grid.boundMax_x && next[0] > grid.boundMax_x){
			normal = Vector3d(-1.f, 0.f, 0.f);
			return true;
		}
	}

	//y collisions
	if(	grid.boundMax_x > prev[0] && prev[0] >= grid.boundMin_x && 
		grid.boundMax_z > prev[2] && prev[2] >= grid.boundMin_z){

		if(prev[1] > grid.boundMin_y && next[1] < grid.boundMin_y){
			normal = Vector3d(0.f, 1.f, 0.f);
			return true;
		}
		else if(prev[1] < grid.boundMin_y && next[1] > grid.boundMin_y){
			normal = Vector3d(0.f, -1.f, 0.f);
			return true;
		}
		if(prev[1] > grid.boundMax_y && next[1] < grid.boundMax_y){
			normal = Vector3d(0.f, 1.f, 0.f);
			return true;
		}
		else if(prev[1] < grid.boundMax_y && next[1] > grid.boundMax_y){
			normal = Vector3d(0.f, -1.f, 0.f);
			return true;
		}
	}
	//z collisions
	if(	grid.boundMax_x > prev[0] && prev[0] >= grid.boundMin_x && 
		grid.boundMax_y > prev[1] && prev[1] >= grid.boundMin_y){

		if(prev[2] > grid.boundMin_z && next[2] < grid.boundMin_z){
			normal = Vector3d(0.f, 0.f,  1.f);
			return true;
		}
		else if(prev[2] < grid.boundMin_z && next[2] > grid.boundMin_z){
			normal = Vector3d(0.f, 0.f,  -1.f);
			return true;
		}
		if(prev[2] > grid.boundMax_z && next[2] < grid.boundMax_z){
			normal = Vector3d(0.f, 0.f,  1.f);
			return true;
		}
		else if(prev[2] < grid.boundMax_z && next[2] > grid.boundMax_z){
			normal = Vector3d(0.f, 0.f,  -1.f);
			return true;
		}
	}
	return false;
}

void collisionBasedParticleUpdate(Particles & particles, MacGrid & grid, double subStep){

	double collisionBoundary_minI = 
	double collisionBoundary_maxI = 

	double collisionBoundary_minJ = 
	double collisionBoundary_maxJ = 

	double collisionBoundary_minK = 
	double collisionBoundary_maxK = 

	for(uint pi = 0 ; pi < particles.num ; pi++){
		auto pVel = particles.getVel(pi);
		auto pPos = particles.getPos(pi);

		Vector3d projectedPos = pPos + pVel * subStep;

		if(pVel[0] > 0 && projectedPos[0] > collisionBoundary_maxI)
			projectedPos[0] = collisionBoundary_maxI - .01f;
		if(pVel[0] < 0 && projectedPos[0] < collisionBoundary_minI)
			projectedPos[0] = collisionBoundary_minI + .01f;

		if(pVel[1] > 0 && projectedPos[1] > collisionBoundary_maxJ)
			projectedPos[1] = collisionBoundary_maxJ - .01f;
		if(pVel[1] < 0 && projectedPos[1] < collisionBoundary_minJ)
			projectedPos[1] = collisionBoundary_minJ + .01f;

		if(pVel[2] > 0 && projectedPos[2] > collisionBoundary_maxK)
			projectedPos[2] = collisionBoundary_maxK - .01f;
		if(pVel[2] < 0 && projectedPos[2] < collisionBoundary_minK)
			projectedPos[2] = collisionBoundary_minK + .01f;

		pPos = projectedPos;
		assert(!std::isnan(pPos[0]));
		assert(!std::isnan(pPos[1]));
		assert(!std::isnan(pPos[2]));
	}
}

inline double projectPressureOnFace(MacGrid & grid, int idxNBR, double cCP){

	double nbrCCP = grid.cellCenterPressure(idxNBR, 0);
	return invWaterDensity * grid.invCellSize * (nbrCCP - cCP);
}

void applyPressureForces(MacGrid & grid, double timeStep){

	grid.cellCenterPressure.col(0).setZero();

	initializeCellCenterDivergence(grid, timeStep);

	grid.cellCenterPressure = grid.cg.solve(grid.cellCenterDivergence);

	//enforceBoundaryPressure(grid);

	DEBUG();
	//double scale = timeStep/grid.cellSize;
	double scale = timeStep/(double)grid.cellSize;

	for ( uint i = 0 ; i < grid.cellNum_i - 0 ; i++)
		for ( uint j = 0 ; j < grid.cellNum_j - 0 ; j++)
			for ( uint k = 0 ; k < grid.cellNum_k - 0 ; k++){

				auto cellCenterIdx = grid.cellCenterIdx(i, j, k);
				auto idxMax_i = grid.cellCenterIdx(i + 1, j, k);
				auto idxMax_j = grid.cellCenterIdx(i, j + 1, k);
				auto idxMax_k = grid.cellCenterIdx(i, j, k + 1);

				auto cellState = grid.cellCenterState.at(cellCenterIdx);
				auto cellStateMin_i = grid.cellCenterState.at(idxMax_i);
				auto cellStateMin_j = grid.cellCenterState.at(idxMax_j);
				auto cellStateMin_k = grid.cellCenterState.at(idxMax_k);

				double & cCP = grid.cellCenterPressure(cellCenterIdx, 0);

				if(cellState == AIR)
					cCP = 0;

				//double cCPMax_i = cellStateMin_i == LIQUID ? 
//
				//bool isValidNBRCell_i = grid.isCellInGrid_i(i + 1) 
				//	&& (grid.getCellState(i + 1, j, k) != SOLID);
//
				//bool isValidNBRCell_j = grid.isCellInGrid_j(j + 1) 
				//	&& (grid.getCellState(i, j + 1, k) != SOLID);
//
				//bool isValidNBRCell_k = grid.isCellInGrid_k(k + 1) 
				//	&& (grid.getCellState(i, j, k + 1) != SOLID);

				double cellFacePressure_maxU = 
				if(grid.isCellInGrid_i(i + 1))
					if(cellState == LIQUID && cellStateMin_i == LIQUID)
						grid.cellFaceVel_u.at(grid.cellFaceIdx<maxU>(i, j, k)) -= invWaterDensity * grid.invCellSize * (grid.cellCenterPressure(idxMax_i, 0) - cCP) * scale;
					else if(cellStateMin_i == AIR)
						grid.cellFaceVel_u.at(grid.cellFaceIdx<maxU>(i, j, k)) -= invWaterDensity * grid.invCellSize * (-cCP) * scale;

				if(grid.isCellInGrid_j(j + 1))
					if(grid.getCellState(i, j + 1, k) == LIQUID)
						grid.cellFaceVel_v.at(grid.cellFaceIdx<maxV>(i, j, k)) -=
							projectPressureOnFace(grid, i, j + 1, k, cCP) * scale;
					else
						grid.cellFaceVel_v.at(grid.cellFaceIdx<maxV>(i, j, k)) -= invWaterDensity * grid.invCellSize * (-cCP) * scale;

				if(grid.isCellInGrid_k(k + 1))
					if(grid.getCellState(i, j, k + 1) == LIQUID)
						grid.cellFaceVel_w.at(grid.cellFaceIdx<maxW>(i, j, k)) -= 
							projectPressureOnFace(grid, i, j, k + 1, cCP) * scale;
					else
						grid.cellFaceVel_w.at(grid.cellFaceIdx<maxW>(i, j, k)) -= invWaterDensity * grid.invCellSize * (-cCP) * scale;

			}
	//enforceBoundaryVelocities(grid);
}

void cIsBound(MacGrid & grid, int i, int j, int k){
	assert(grid.getCellState(i, j, k) == SOLID);
}




				auto maxIdx_i = grid.cellFaceIdx<maxU>(i, j, k);
				auto maxIdx_j = grid.cellFaceIdx<maxV>(i, j, k);
				auto maxIdx_k = grid.cellFaceIdx<maxW>(i, j, k);
				if(maxIdx_i > grid.cellFaceVel_u.size()
				|| maxIdx_j > grid.cellFaceVel_v.size()
				|| maxIdx_k > grid.cellFaceVel_w.size())
					continue;
				auto cellCenterIdx = grid.cellCenterIdx(i, j, k);
				double cCP = grid.cellCenterPressure(cellCenterIdx, 0);

				bool cellIsSolid = grid.isSolidCell(i, j, k);

				bool NBRIsSolid_maxI = grid.isSolidCell(i + 1, j, k);
				bool NBRIsSolid_maxJ = grid.isSolidCell(i, j + 1, k);
				bool NBRIsSolid_maxK = grid.isSolidCell(i, j, k + 1);

				double newCellFaceVel_u = grid.cellFaceVel_u.at(maxIdx_i) 
					- projectPressureOnFace(grid, i + 1, j, k, cCP, invCellSize);
				double newCellFaceVel_v = grid.cellFaceVel_v.at(maxIdx_j) 
					- projectPressureOnFace(grid, i, j + 1, k, cCP, invCellSize);
				double newCellFaceVel_w = grid.cellFaceVel_w.at(maxIdx_k) 
					- projectPressureOnFace(grid, i, j, k + 1, cCP, invCellSize);

				assert(!std::isnan(newCellFaceVel_u));
				assert(!std::isnan(newCellFaceVel_v));
				assert(!std::isnan(newCellFaceVel_w));

				grid.cellFaceVel_u.at(maxIdx_i) = newCellFaceVel_u;
				grid.cellFaceVel_v.at(maxIdx_j) = newCellFaceVel_v;
				grid.cellFaceVel_w.at(maxIdx_k) = newCellFaceVel_w;

				//if(!cellIsSolid && !NBRIsSolid_maxI && !NBRIsSolid_maxJ && !NBRIsSolid_maxK){
				//	grid.cellFaceVel_u.at(maxIdx_i) = newCellFaceVel_u;
				//	grid.cellFaceVel_v.at(maxIdx_j) = newCellFaceVel_v;
				//	grid.cellFaceVel_w.at(maxIdx_k) = newCellFaceVel_w;
				//	continue;
				//}
				//if(!cellIsSolid){
				//	if(NBRIsSolid_maxI){
				//		grid.cellFaceVel_u.at(maxIdx_i) = std::min<double>(0.f, newCellFaceVel_u);
				//	}else{
				//		grid.cellFaceVel_u.at(maxIdx_i) = std::max<double>(0.f, newCellFaceVel_u);
				//	}
				//
				//	if(NBRIsSolid_maxJ){
				//		grid.cellFaceVel_v.at(maxIdx_j) = std::min<double>(0.f, newCellFaceVel_v);
				//	}else{
				//		grid.cellFaceVel_v.at(maxIdx_j) = std::max<double>(0.f, newCellFaceVel_v);
				//	}
				//
				//	if(NBRIsSolid_maxK){
				//		grid.cellFaceVel_w.at(maxIdx_k) = std::min<double>(0.f, newCellFaceVel_w);
				//	}else{
				//		grid.cellFaceVel_w.at(maxIdx_k) = std::max<double>(0.f, newCellFaceVel_w);
				//	}
				//}
				//else{
				//	if(NBRIsSolid_maxI){
				//		grid.cellFaceVel_u.at(maxIdx_i) = std::max<double>(0.f, newCellFaceVel_u);
				//	}else{
				//		grid.cellFaceVel_u.at(maxIdx_i) = newCellFaceVel_u;
				//	}
				//
				//	if(NBRIsSolid_maxJ){
				//		grid.cellFaceVel_v.at(maxIdx_j) = std::max<double>(0.f, newCellFaceVel_v);
				//	}else{
				//		grid.cellFaceVel_v.at(maxIdx_j) = newCellFaceVel_v;
				//	}
				//
				//	if(NBRIsSolid_maxK){
				//		grid.cellFaceVel_w.at(maxIdx_k) = std::max<double>(0.f, newCellFaceVel_w);
				//	}else{
				//		grid.cellFaceVel_w.at(maxIdx_k) = newCellFaceVel_w;
				//	}
				//}
**/


//=======================================================================================================

/**
	for ( uint i = 0 ; i < grid.cellNum_i ; i++ )
		for ( uint j = 0 ; j < grid.cellNum_j ; j += grid.cellNum_j - 1)
			for ( uint k = 0 ; k < grid.cellNum_k ; k++){

			}

	for ( uint i = 0 ; i < grid.cellNum_i ; i++ )
		for ( uint j = 0 ; j < grid.cellNum_j ; j++ )
			for ( uint k = 0 ; k < grid.cellNum_k ; k += grid.cellNum_k -1){
				
			}

	for ( uint i = 0 ; i < grid.cellNum_i ; i++ )
		for ( uint j = 0 ; j < grid.cellNum_j ; j += jump_j)
			for ( uint k = 0 ; k < grid.cellNum_k ; k += jump_k){
				if(1 < i && i < grid.cellNum_i)
					jump_j = grid.cellNum_j;
				else
					jump_j = 1;
				if(1 < j && j < grid.cellNum_j)
					jump_k = grid.cellNum_k;
				else
					jump_k = 1;

				

				if(!grid.cellIsOutOfBounds(i - 1, j, k))
					if(grid.isSolidCell(i - 1, j, k)){

					}

				if(!grid.cellIsOutOfBounds(i + 1, j, k))
					if(grid.isSolidCell(i + 1, j, k)){

					}

				if(!grid.cellIsOutOfBounds(i, j - 1, k))
					if(grid.isSolidCell(i, j - 1, k)){

					}

				if(!grid.cellIsOutOfBounds(i, j + 1, k))
					if(grid.isSolidCell(i, j + 1, k)){

					}

				if(!grid.cellIsOutOfBounds(i, j, k - 1))
					if(grid.isSolidCell(i, j, k - 1)){

					}

				if(!grid.cellIsOutOfBounds(i, j, k + 1))
					if(grid.isSolidCell(i, j, k + 1)){

					}


				if(isMaxBound_i && !isMaxBound_j && !isMinBound_j && !isMaxBound_k && !isMinBound_k){
					grid.cellFaceVel<minU>(i, j, k) = std::min<double>(0.f, grid.cellFaceVel<minU>(i, j, k));
					grid.cellFaceVel<maxV>(i, j, k) = std::min<double>(0.f, grid.cellFaceVel<maxV>(i - 1, j, k));
					grid.cellFaceVel<maxW>(i, j, k) = std::min<double>(0.f, grid.cellFaceVel<maxW>(i - 1, j, k));
				}
				else if(isMinBound_i && !isMaxBound_j && !isMinBound_j && !isMaxBound_k && !isMinBound_k){
					grid.cellFaceVel<maxU>(i, j, k) = std::max<double>(0.f, grid.cellFaceVel<maxU>(i, j, k));
					grid.cellFaceVel<maxV>(i, j, k) = std::max<double>(0.f, grid.cellFaceVel<maxV>(i + 1, j, k));
					grid.cellFaceVel<maxW>(i, j, k) = std::max<double>(0.f, grid.cellFaceVel<maxW>(i + 1, j, k));
				}
				else if(!isMaxBound_i && !isMinBound_i && isMaxBound_j && !isMaxBound_k && !isMinBound_k){
					grid.cellFaceVel<maxU>(i, j, k) = std::min<double>(0.f, grid.cellFaceVel<maxU>(i, j, k));
					grid.cellFaceVel<minV>(i, j, k) = std::min<double>(0.f, grid.cellFaceVel<minV>(i, j - 1, k));
					grid.cellFaceVel<maxW>(i, j, k) = std::min<double>(0.f, grid.cellFaceVel<maxW>(i, j - 1, k));
				}
				else if(isMinBound_i && !isMaxBound_j && !isMinBound_j && !isMaxBound_k && !isMinBound_k){
					grid.cellFaceVel<maxU>(i, j, k) = std::max<double>(0.f, grid.cellFaceVel<maxU>(i, j, k));
					grid.cellFaceVel<maxV>(i, j, k) = std::max<double>(0.f, grid.cellFaceVel<maxV>(i, j + 1, k));
					grid.cellFaceVel<maxW>(i, j, k) = std::max<double>(0.f, grid.cellFaceVel<maxW>(i, j + 1, k));
				}
				uint faceIdx_u
			}

	for ( uint i = 0 ; i < grid.cellNum_i ; i++)
		for ( uint j = 0 ; j < grid.cellNum_j ; j++)
			for ( uint k = 0 ; k < grid.cellNum_k ; k++){

				bool isSolid = grid.isSolidCell(i, j, k);
				bool nbrIsSolid_maxI = grid.isSolidCell(i + 1, j, k);
				bool nbrIsSolid_maxJ = grid.isSolidCell(i, j + 1, k);
				bool nbrIsSolid_maxK = grid.isSolidCell(i, j, k + 1);

				if(!isSolid && !nbrIsSolid_maxI && !nbrIsSolid_maxJ && !nbrIsSolid_maxK)
					continue;


				if(!cellIsSolid && !NBRIsSolid_maxI && !NBRIsSolid_maxJ && !NBRIsSolid_maxK){
					grid.cellFaceVel_u.at(maxIdx_i) = newCellFaceVel_u;
					grid.cellFaceVel_v.at(maxIdx_j) = newCellFaceVel_v;
					grid.cellFaceVel_w.at(maxIdx_k) = newCellFaceVel_w;
					continue;
				}
				if(!cellIsSolid){
					if(NBRIsSolid_maxI){
						grid.cellFaceVel_u.at(maxIdx_i) = std::min<double>(0.f, newCellFaceVel_u);
					}else{
						grid.cellFaceVel_u.at(maxIdx_i) = std::max<double>(0.f, newCellFaceVel_u);
					}
				
					if(NBRIsSolid_maxJ){
						grid.cellFaceVel_v.at(maxIdx_j) = std::min<double>(0.f, newCellFaceVel_v);
					}else{
						grid.cellFaceVel_v.at(maxIdx_j) = std::max<double>(0.f, newCellFaceVel_v);
					}
				
					if(NBRIsSolid_maxK){
						grid.cellFaceVel_w.at(maxIdx_k) = std::min<double>(0.f, newCellFaceVel_w);
					}else{
						grid.cellFaceVel_w.at(maxIdx_k) = std::max<double>(0.f, newCellFaceVel_w);
					}
				}
				else{
					if(NBRIsSolid_maxI){
						grid.cellFaceVel_u.at(maxIdx_i) = std::max<double>(0.f, newCellFaceVel_u);
					}else{
						grid.cellFaceVel_u.at(maxIdx_i) = newCellFaceVel_u;
					}
				
					if(NBRIsSolid_maxJ){
						grid.cellFaceVel_v.at(maxIdx_j) = std::max<double>(0.f, newCellFaceVel_v);
					}else{
						grid.cellFaceVel_v.at(maxIdx_j) = newCellFaceVel_v;
					}
				
					if(NBRIsSolid_maxK){
						grid.cellFaceVel_w.at(maxIdx_k) = std::max<double>(0.f, newCellFaceVel_w);
					}else{
						grid.cellFaceVel_w.at(maxIdx_k) = newCellFaceVel_w;
					}
				}
			}

    _MAC->_gType.iterate([&](size_t i, size_t j, size_t k) {
        switch (_MAC->_gType(i,j,k)) {
            case EMPTY:break;
            case FLUID:break;
            case SOLID:
                if (i == 0 || _MAC->_gType(i-1,j,k) != SOLID) {
                    _MAC->_gU(i,j,k) = std::min(0.f, _MAC->_gU(i,j,k));
                }
                if (i == _MAC->_gType.countX() - 1 || _MAC->_gType(i+1,j,k) != SOLID) {
                    _MAC->_gU(i+1,j,k) = std::max(0.f, _MAC->_gU(i+1,j,k));
                }
                if (j == 0 || _MAC->_gType(i,j-1,k) != SOLID) {
                    _MAC->_gV(i,j,k) = std::min(0.f, _MAC->_gV(i,j,k));
                }
                if (j == _MAC->_gType.countY() - 1 || _MAC->_gType(i,j+1,k) != SOLID) {
                    _MAC->_gV(i,j+1,k) = std::max(0.f, _MAC->_gV(i,j+1,k));
                }
                if (k == 0 || _MAC->_gType(i,j,k-1) != SOLID) {
                    _MAC->_gW(i,j,k) = std::min(0.f, _MAC->_gW(i,j,k));
                }
                if (k == _MAC->_gType.countZ() - 1 || _MAC->_gType(i,j,k+1) != SOLID) {
                    _MAC->_gW(i,j,k+1) = std::max(0.f, _MAC->_gW(i,j,k+1));
                }
                break;
            default:break;
        }
    });
**/

/**

void initializeLaplacianNBRMat(MacGrid & grid){

	std::vector<Triplet<int>> tripletsVector;


	for ( uint i = 0 ; i < grid.cellNum_i ; i++)
		for ( uint j = 0 ; j < grid.cellNum_j ; j++)
			for ( uint k = 0 ; k < grid.cellNum_k ; k++){
				auto cellCenterIdx = grid.cellCenterIdx(i, j, k);
				int totalValidNeighbours = 0;
				//auto lapNeighbourEntry = [](int i, int j, int k) -> void
				//{
				//	if (grid.isBoundaryCell(i, j, k)){
				//		auto cellCenterMinIdx_i = grid.cellCenterIdx(i, j, k);
				//		tripletsVector.resize(tripletsVector.size() + 2);
				//		tripletsVector[tripletsVector.size() - 2 ] = Triplet<int>(cellCenterIdx, cellCenterMinIdx_i, 1);
				//		tripletsVector[tripletsVector.size() - 1 ] = Triplet<int>(cellCenterMinIdx_i,cellCenterIdx , 1);
				//		totalValidNeighbours += 1;
				//	}
				//};

				if(!grid.isBoundaryCell(i - 1, j, k)){
					auto cellCenterMinIdx_i = grid.cellCenterIdx(i - 1, j, k);
					tripletsVector.push_back(Triplet<int>(cellCenterIdx, cellCenterMinIdx_i, 1));
					tripletsVector.push_back(Triplet<int>(cellCenterMinIdx_i, cellCenterIdx, 1));
					totalValidNeighbours += 1;
				}

				if(!grid.isBoundaryCell(i + 1, j, k)){
					auto cellCenterMaxIdx_i = grid.cellCenterIdx(i + 1, j, k);
					tripletsVector.push_back(Triplet<int>(cellCenterIdx, cellCenterMaxIdx_i, 1));
					tripletsVector.push_back(Triplet<int>(cellCenterMaxIdx_i, cellCenterIdx, 1));
					totalValidNeighbours += 1;
				}

				if(!grid.isBoundaryCell(i, j - 1, k)){
					auto cellCenterMinIdx_j = grid.cellCenterIdx(i, j - 1, k);
					tripletsVector.push_back(Triplet<int>(cellCenterIdx, cellCenterMinIdx_j, 1));
					tripletsVector.push_back(Triplet<int>(cellCenterMinIdx_j, cellCenterIdx, 1));
					totalValidNeighbours += 1;
				}

				if(!grid.isBoundaryCell(i, j + 1, k)){
					auto cellCenterMaxIdx_j = grid.cellCenterIdx(i, j + 1, k);
					tripletsVector.push_back(Triplet<int>(cellCenterIdx, cellCenterMaxIdx_j, 1));
					tripletsVector.push_back(Triplet<int>(cellCenterMaxIdx_j, cellCenterIdx, 1));
					totalValidNeighbours += 1;
				}

				if(!grid.isBoundaryCell(i, j, k - 1)){
					auto cellCenterMinIdx_k = grid.cellCenterIdx(i, j, k - 1);
					tripletsVector.push_back(Triplet<int>(cellCenterIdx, cellCenterMinIdx_k, 1));
					tripletsVector.push_back(Triplet<int>(cellCenterMinIdx_k, cellCenterIdx, 1));
					totalValidNeighbours += 1;
				}

				if(!grid.isBoundaryCell(i, j, k + 1)){
					auto cellCenterMaxIdx_k = grid.cellCenterIdx(i, j, k + 1);
					tripletsVector.push_back(Triplet<int>(cellCenterIdx, cellCenterMaxIdx_k, 1));
					tripletsVector.push_back(Triplet<int>(cellCenterMaxIdx_k, cellCenterIdx, 1));
					totalValidNeighbours += 1;
				}
				tripletsVector.push_back(Triplet<int>(cellCenterIdx, cellCenterIdx, -totalValidNeighbours));
			}

	grid.laplacianSparseNBR.setFromTriplets(tripletsVector.begin(), tripletsVector.end());
	DEBUG();
	checkForNan(grid.laplacianSparseNBR);

	std::cout<<"grid.laplacianSparseNBR.rows():\n"<<grid.laplacianSparseNBR.rows()<<std::endl;
	std::cout<<"grid.laplacianSparseNBR.cols():\n"<<grid.laplacianSparseNBR.cols()<<std::endl;
	//grid.cg = new CholmodSupernodalLLT<SparseMatrix<double>>(grid.laplacianSparseNBR);
	
	std::cout<<grid.laplacianSparseNBR<<std::endl;
	grid.cg.compute(grid.laplacianSparseNBR);
	DEBUG();
}

void applyPressureForces(MacGrid & grid, double subStep){
	initializeCellCenterDivergence(grid, timeStep);
	
	grid.cellCenterPressure = grid.cg.solve(grid.cellCenterDivergence);
	
	double invWaterDensity = 1.f/waterDensity;
	double invCellSize = 1.f/grid.cellSize;


	for ( uint i = 0 ; i < grid.cellNum_i ; i++)
		for ( uint j = 0 ; j < grid.cellNum_j ; j++)
			for ( uint k = 0 ; k < grid.cellNum_k ; k++){

					auto cellCenterIdx = grid.cellCenterIdx(i, j, k);

					auto max_i = grid.cellFaceIdx<maxU>(i, j, k);
					auto min_i = grid.cellFaceIdx<minU>(i, j, k);

					auto max_j = grid.cellFaceIdx<maxV>(i, j, k);
					auto min_j = grid.cellFaceIdx<minV>(i, j, k);

					auto max_k = grid.cellFaceIdx<maxW>(i, j, k);
					auto min_k = grid.cellFaceIdx<minW>(i, j, k);

					double cCP = grid.cellCenterPressure(cellCenterIdx, 0);

					double cCP_i = grid.cellCenterPressure(grid.cellCenterIdx(i + 1, j, k), 0);
					double cCP_j = grid.cellCenterPressure(grid.cellCenterIdx(i, j + 1, k), 0);
					double cCP_k = grid.cellCenterPressure(grid.cellCenterIdx(i, j, k + 1), 0);

					grid.cellFaceVel_u.at(max_i) += subStep * invWaterDensity * (cCP_i - cCP) * invCellSize;
					grid.cellFaceVel_v.at(max_j) += subStep * invWaterDensity * (cCP_j - cCP) * invCellSize;
					grid.cellFaceVel_w.at(max_k) += subStep * invWaterDensity * (cCP_k - cCP) * invCellSize;
			}
}

void transferAttributes(MacGrid const & grid, Particles & particles){
	//pi is particle index
	for(int pi = 0 ; pi < particles.num ; pi++){
		auto & pPos = particles.positions[pi];
		auto & pVel = particles.velocities[pi];
		auto [i, j, k] = grid.gridCoord(pPos);
		auto [min_i, min_j, min_k, max_i, max_j, max_k] = getLocalNBRIdcs(grid, pPos, i, j, k);

		pVel = {0,0,0};

		double cornerCoeff = 1.f/8.f;

		auto cellFacePosc000_u = grid.cellFacePos<minU>(min_i, min_j, min_k);
		auto cellFacePosc111_u = grid.cellFacePos<minU>(max_i, max_j, max_k);
		tuple8i cellUFaceNBRIdcs = getStageredCellFaceNBRIdcs_u(grid, min_i, min_j, min_k, max_i, max_j, max_k);
		auto weights_u = getWeights(pPos, cellFacePosc000_u, cellFacePosc111_u);
		std::cout<<"weights_u:\n";
		for(int ii=0 ; ii<8 ; ii++){
			pVel[0] += grid.cellFaceVel_u[cellUFaceNBRIdcs[ii]] * cornerCoeff / weights_u[ii];
			std::cout<<weights_u[ii]<<"\n";
		}

		auto cellFacePosc000_v = grid.cellFacePos<minV>(min_i, min_j, min_k);
		auto cellFacePosc111_v = grid.cellFacePos<minV>(max_i, max_j, max_k);
		tuple8i cellVFaceNBRIdcs = getStageredCellFaceNBRIdcs_v(grid, min_i, min_j, min_k, max_i, max_j, max_k);
		auto weights_v = getWeights(pPos, cellFacePosc000_v, cellFacePosc111_v);
		std::cout<<"weights_v:\n";
		for(int ii=0 ; ii<8 ; ii++){
			pVel[1] += grid.cellFaceVel_v[cellVFaceNBRIdcs[ii]] * cornerCoeff / weights_v[ii];
			std::cout<<weights_v[ii]<<"\n";
		}

		auto cellFacePosc000_w = grid.cellFacePos<minW>(min_i, min_j, min_k);
		auto cellFacePosc111_w = grid.cellFacePos<minW>(max_i, max_j, max_k);
		tuple8i cellWFaceNBRIdcs = getStageredCellFaceNBRIdcs_w(grid, min_i, min_j, min_k, max_i, max_j, max_k);
		auto weights_w = getWeights(pPos, cellFacePosc000_w, cellFacePosc111_w);
		std::cout<<"weights_w:\n";
		for(int ii=0 ; ii<8 ; ii++){
			pVel[2] += grid.cellFaceVel_w[cellWFaceNBRIdcs[ii]] * cornerCoeff / weights_w[ii];
			std::cout<<weights_w[ii]<<"\n";
		}
	}
}

**/

}