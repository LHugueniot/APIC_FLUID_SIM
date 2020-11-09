#include "MacGrid.h"

namespace pic{

//Get position relative to grid origin and grid scale
Vector3d MacGrid::gridSpacePos(Vector3d const & worldSpacePos) const{
	return (worldSpacePos - origin) * invCellSize;
}
//Get position relative to cell origin and grid scale
Vector3d MacGrid::cellSpacePos(Vector3d const & worldSpacePos, int i, int j, int k) const{
	return (worldSpacePos - origin) * invCellSize - Vector3d(i,j,k);
}
//Get center position of a cell in worldspace
Vector3d MacGrid::cellCenterPos(int i, int j, int k) const{
	return {((double)i + .5f) * cellSize + origin[0], 
		((double)j + .5f) * cellSize + origin[1], 
		((double)k + .5f) * cellSize + origin[2]};
}
tuple3i MacGrid::gridCoord(Vector3d const & worldSpacePos) const{
	Vector3d const gridSpacePos = (worldSpacePos - origin) * invCellSize;
	return {(int)std::floor(gridSpacePos[0]), (int)std::floor(gridSpacePos[1]), (int)std::floor(gridSpacePos[2])};
}
int MacGrid::cellCenterIdx(int i, int j, int k) const{
	return (k * (cellNum_i * cellNum_j) + j * cellNum_i + i);
}
int MacGrid::cellCenterIdx(Vector3d const & worldSpacePos){
	auto [i, j, k] = gridCoord(worldSpacePos);
	return cellCenterIdx(i, j, k);
}

bool MacGrid::isCellInGrid_i(int i) const{
	if ( (i < 0) || ((int)cellNum_i <= i) )
		return false;
	return true;
}
bool MacGrid::isCellInGrid_j(int j) const{
	if ( (j < 0) || ((int)cellNum_j <= j) )
		return false;
	return true;
}
bool MacGrid::isCellInGrid_k(int k) const{
	if ( (k < 0) || ((int)cellNum_k <= k) )
		return false;
	return true;
}

bool MacGrid::cellIsOnBounds(int i, int j, int k) const{
	if ( (i == 0) || ((int)cellNum_i - 1 == i) || (j == 1) || ((int)cellNum_j - 1 == j) || (k == 1) || ((int)cellNum_k - 1 == k))
		return true;
	return false;
}

bool MacGrid::cellIsOutOfBounds(int i, int j, int k) const{
	if ( (i < 0) || ((int)cellNum_i <= i) || (j < 0) || ((int)cellNum_j <= j) || (k < 0) || ((int)cellNum_k <= k))
		return true;
	return false;
}
bool MacGrid::cellIsInBounds(int i, int j, int k) const{
	if ( (0 <= i) && (i < (int)cellNum_i) && (0 <= j) && (j < (int)cellNum_j) && (0 <= k) && (k < (int)cellNum_k))
		return true;
	return false;
}

bool MacGrid::isMaxBoundaryCell_i(int i) const{
	if(i == (int)cellNum_i - 1)
		return true;
	return false;
}
bool MacGrid::isMinBoundaryCell_i(int i) const{
	if(!i)
		return true;
	return false;
}

bool MacGrid::isMaxBoundaryCell_j(int j) const{
	if(j == (int)cellNum_j - 1)
		return true;
	return false;
}
bool MacGrid::isMinBoundaryCell_j(int j) const{
	if(!j)
		return true;
	return false;
}

bool MacGrid::isMaxBoundaryCell_k(int k) const{
	if(k == (int)cellNum_k - 1)
		return true;
	return false;
}
bool MacGrid::isMinBoundaryCell_k(int k) const{
	if(!k)
		return true;
	return false;
}

template<>
bool MacGrid::isValidFace<minU>(int i, int j, int k) const{
	return ((0 <= (i - 1)) && ((i - 1) < (int)cellFaceNum_i))
		&& ((0 <= j) && (j < (int)cellNum_j)) 
		&& ((0 <= k) && (k < (int)cellNum_k));
}
template<>
bool MacGrid::isValidFace<minV>(int i, int j, int k) const{
	return ((0 <= i) && (i < (int)cellNum_i)) 
		&& ((0 <= (j - 1)) && ((j - 1) < (int)cellFaceNum_j))
		&& ((0 <= k) && (k < (int)cellNum_k));
}
template<>
bool MacGrid::isValidFace<minW>(int i, int j, int k) const{
	return ((0 <= i) && (i < (int)cellNum_i))
		&& ((0 <= j) && (j < (int)cellNum_j))
		&& ((0 <= (k - 1)) && ((k - 1) < (int)cellFaceNum_k));
}

template<>
bool MacGrid::isValidFace<maxU>(int i, int j, int k) const{
	return ((0 <= i) && (i < (int)cellFaceNum_i))
		&& ((0 <= j) && (j < (int)cellNum_j)) 
		&& ((0 <= k) && (k < (int)cellNum_k));
}
template<>
bool MacGrid::isValidFace<maxV>(int i, int j, int k) const{
	return ((0 <= i) && (i < (int)cellNum_i))
		&& ((0 <= j) && (j < (int)cellFaceNum_j)) 
		&& ((0 <= k) && (k < (int)cellNum_k));
}
template<>
bool MacGrid::isValidFace<maxW>(int i, int j, int k) const{
	return ((0 <= i) && (i < (int)cellNum_i))
		&& ((0 <= j) && (j < (int)cellNum_j))
		&& ((0 <= k) && (k < (int)cellFaceNum_k));
}

//========================================================================


template<>
int MacGrid::cellFaceIdx<minU>(int i, int j, int k) const{
	#if EZ_MODE
		assert(isValidFace<minU>(i,j,k) == true);
	#endif
	return (k * cellFaceNum_i * cellNum_j) + (j * cellFaceNum_i) + (i - 1);
}
template<>
int MacGrid::cellFaceIdx<minV>(int i, int j, int k) const{
	#if EZ_MODE
		assert(isValidFace<minV>(i,j,k) == true);
	#endif
	return (i * cellFaceNum_j * cellNum_k) + (k * cellFaceNum_j) + (j - 1);
}
template<>
int MacGrid::cellFaceIdx<minW>(int i, int j, int k) const{
	#if EZ_MODE
		assert(isValidFace<minW>(i,j,k) == true);
	#endif
	return (j * cellFaceNum_k * cellNum_i) + (i * cellFaceNum_k) + (k - 1);
}

	//return (k * (cellNum_i * cellNum_j) + j * cellNum_i + i);

template<>
int MacGrid::cellFaceIdx<maxU>(int i, int j, int k) const{
	#if EZ_MODE
		assert(isValidFace<maxU>(i,j,k) == true);
	#endif
	return (k * cellFaceNum_i * cellNum_j) + (j * cellFaceNum_i) + i;
}
template<>
int MacGrid::cellFaceIdx<maxV>(int i, int j, int k) const{
	#if EZ_MODE
		assert(isValidFace<maxV>(i,j,k) == true);
	#endif
	return (i * cellFaceNum_j * cellNum_k) + (k * cellFaceNum_j) + j;
}

template<>
int MacGrid::cellFaceIdx<maxW>(int i, int j, int k) const{
	#if EZ_MODE
		assert(isValidFace<maxW>(i,j,k) == true);
	#endif
	return (j * cellFaceNum_k * cellNum_i) + (i * cellFaceNum_k) + k;
}

template<>
double & MacGrid::cellFaceVel<minU>(int i, int j, int k){
	return cellFaceVel_u.at(cellFaceIdx<minU>(i,j,k));
}
template<>
double & MacGrid::cellFaceVel<minV>(int i, int j, int k){
	return cellFaceVel_v.at(cellFaceIdx<minV>(i,j,k));
}
template<>
double & MacGrid::cellFaceVel<minW>(int i, int j, int k){
	return cellFaceVel_w.at(cellFaceIdx<minW>(i,j,k));
}

template<>
double & MacGrid::cellFaceVel<maxU>(int i, int j, int k){
	return cellFaceVel_u.at(cellFaceIdx<maxU>(i,j,k));
}
template<>
double & MacGrid::cellFaceVel<maxV>(int i, int j, int k){
	return cellFaceVel_v.at(cellFaceIdx<maxV>(i,j,k));
}
template<>
double & MacGrid::cellFaceVel<maxW>(int i, int j, int k){
	return cellFaceVel_w.at(cellFaceIdx<maxW>(i,j,k));
}

void MacGrid::printOutMaxFaceVels(){
	std::cout<<"cellFaceVel max"<<std::endl;
	for ( int j = 0 ; j < (int)cellNum_j ; j++){
		for ( int i = 0 ; i < (int)cellNum_i ; i++){
			for ( int k = 0 ; k < (int)cellNum_k ; k++){
				std::cout<<"{";
				if(isValidFace<maxU>(i, j, k))
					std::cout<<cellFaceVel<maxU>(i, j, k)<<", ";
				else
					std::cout<<"_, ";
				if(isValidFace<maxV>(i, j, k))
					std::cout<<cellFaceVel<maxV>(i, j, k)<<", ";
				else
					std::cout<<"_, ";
				if(isValidFace<maxW>(i, j, k))
					std::cout<<cellFaceVel<maxW>(i, j, k);
				else
					std::cout<<"_";
				std::cout<<"} ";
			}
			std::cout<<std::endl;
		}
		std::cout<<std::endl;
	}
}

void MacGrid::printOutMaxValidFaces(){
	std::cout<<"cellFaceVel max"<<std::endl;
	for ( uint j = 0 ; j < cellNum_j ; j++){
		for ( uint i = 0 ; i < cellNum_i ; i++){
			for ( uint k = 0 ; k < cellNum_k ; k++){
				std::cout<<"{"<<isValidFace<maxU>(i, j, k)<<", "<<isValidFace<maxV>(i, j, k)<<", "<<isValidFace<maxW>(i, j, k)<<"} ";
			}
			std::cout<<std::endl;
		}
		std::cout<<std::endl;
	}
}

void MacGrid::printOutCellStates(){
	std::cout<<"cell states"<<std::endl;
	for ( uint j = 0 ; j < cellNum_j ; j++){
		for ( uint i = 0 ; i < cellNum_i ; i++){
			for ( uint k = 0 ; k < cellNum_k ; k++){
				std::cout<<getCellState(i,j,k)<<" ";
			}
			std::cout<<std::endl;
		}
		std::cout<<std::endl;
	}
}

//Get Face indices as tuples
tuple3i MacGrid::cellMinFaceIdcs(int i, int j, int k) const{
	return {cellFaceIdx<minU>(i, j, k),
			cellFaceIdx<minV>(i, j, k),
			cellFaceIdx<minW>(i, j, k)};
}
tuple3i MacGrid::cellMaxFaceIdcs(int i, int j, int k) const{
	return {cellFaceIdx<maxU>(i, j, k),
			cellFaceIdx<maxV>(i, j, k),
			cellFaceIdx<maxW>(i, j, k)};
}
tuple6i MacGrid::cellFaceIdcs(int i, int j, int k) const{
	return {cellFaceIdx<minU>(i, j, k), cellFaceIdx<minV>(i, j, k), cellFaceIdx<minW>(i, j, k),
			cellFaceIdx<maxU>(i, j, k), cellFaceIdx<maxV>(i, j, k), cellFaceIdx<maxW>(i, j, k)};
}
tuple6i MacGrid::cellFaceIdcs(Vector3d const & worldSpacePos) const{
	auto [i, j, k] = gridCoord(worldSpacePos);
	return cellFaceIdcs(i, j, k);
}


template<>
Vector3d MacGrid::cellFacePos<minU>(int i, int j, int k) const{
	return {(double)i * cellSize + origin[0],
			((double)j + .5f) * cellSize + origin[1],
			((double)k + .5f) * cellSize + origin[2]};
}
template<>
Vector3d MacGrid::cellFacePos<minV>(int i, int j, int k) const{
	return {((double)i + .5f) * cellSize + origin[0], 
			(double)j * cellSize + origin[1],
			((double)k + .5f) * cellSize + origin[2]};
}
template<>
Vector3d MacGrid::cellFacePos<minW>(int i, int j, int k) const{
	return {((double)i + .5f) * cellSize + origin[0],
			((double)j + .5f) * cellSize + origin[1],
			(double)k * cellSize + origin[2]};
}

template<>
Vector3d MacGrid::cellFacePos<maxU>(int i, int j, int k) const{
	return {((double)i + 1.f) * cellSize + origin[0],
			((double)j + .5f) * cellSize + origin[1],
			((double)k + .5f) * cellSize + origin[2]};
}
template<>
Vector3d MacGrid::cellFacePos<maxV>(int i, int j, int k) const{
	return {((double)i + .5f) * cellSize + origin[0],
			((double)j + 1.f) * cellSize + origin[1],
			((double)k + .5f) * cellSize + origin[2]};
}
template<>
Vector3d MacGrid::cellFacePos<maxW>(int i, int j, int k) const{
	return {((double)i + .5f) * cellSize + origin[0],
			((double)j + .5f) * cellSize + origin[1],
			((double)k + 1.f) * cellSize + origin[2]};
}


//Get Min Faces Positions
Vector3d MacGrid::cellMinFacePos_u(int i, int j, int k) const{
	return {(double)i * cellSize + origin[0],
			((double)j + .5f) * cellSize + origin[1],
			((double)k + .5f) * cellSize + origin[2]};
}
Vector3d MacGrid::cellMinFacePos_v(int i, int j, int k) const{
	return {((double)i + .5f) * cellSize + origin[0], 
			(double)j * cellSize + origin[1],
			((double)k + .5f) * cellSize + origin[2]};
}
Vector3d MacGrid::cellMinFacePos_w(int i, int j, int k) const{
	return {((double)i + .5f) * cellSize + origin[0],
			((double)j + .5f) * cellSize + origin[1],
			(double)k * cellSize + origin[2]};
}

//Get MAX Faces Positions
Vector3d MacGrid::cellMaxFacePos_u(int i, int j, int k) const{
	return {((double)i + 1.f) * cellSize + origin[0],
			((double)j + .5f) * cellSize + origin[1],
			((double)k + .5f) * cellSize + origin[2]};
}
Vector3d MacGrid::cellMaxFacePos_v(int i, int j, int k) const{
	return {((double)i + .5f) * cellSize + origin[0],
			((double)j + 1.f) * cellSize + origin[1],
			((double)k + .5f) * cellSize + origin[2]};
}
Vector3d MacGrid::cellMaxFacePos_w(int i, int j, int k) const{
	return {((double)i + .5f) * cellSize + origin[0],
			((double)j + .5f) * cellSize + origin[1],
			((double)k + 1.f) * cellSize + origin[2]};
}

Vector3d MacGrid::cellCenterVel(int i, int j, int k) const{
	return {(double)(cellFaceVel_u[cellFaceIdx<minU>(i, j, k)] + cellFaceVel_u[cellFaceIdx<maxU>(i, j, k)] ) * .5f,
			(double)(cellFaceVel_v[cellFaceIdx<minV>(i, j, k)] + cellFaceVel_v[cellFaceIdx<maxV>(i, j, k)] ) * .5f,
			(double)(cellFaceVel_w[cellFaceIdx<minW>(i, j, k)] + cellFaceVel_w[cellFaceIdx<maxW>(i, j, k)] ) * .5f};
}

Vector3d MacGrid::cellCenterVel(Vector3d const & worldSpacePos) const{
	auto [i, j, k] = gridCoord(worldSpacePos);
	return cellCenterVel(i, j, k);
}

CellState & MacGrid::getCellState(int i, int j, int k){
	return cellCenterState.at(cellCenterIdx(i, j, k));
}

}