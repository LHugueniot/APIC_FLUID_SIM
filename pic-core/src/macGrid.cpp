#include "macGrid.h"

namespace pic
{

//Get position relative to grid origin and grid scale
Vector3d MacGrid::gridSpacePos(Vector3d const & worldSpacePos) const{
	return (worldSpacePos - origin) / cellSize;
}
//Get position relative to cell origin and grid scale
Vector3d MacGrid::cellSpacePos(Vector3d const & worldSpacePos, int i, int j, int k) const{
	return (worldSpacePos - origin) / cellSize - Vector3d(i,j,k);
}
//Get center position of a cell in worldspace
Vector3d MacGrid::cellCenterPos(int i, int j, int k) const{
	return {((double)i + .5f) * cellSize + origin[0], 
		((double)j + .5f) * cellSize + origin[1], 
		((double)k + .5f) * cellSize + origin[2]};
}
tuple3i MacGrid::gridCoord(Vector3d const & worldSpacePos) const{
	Vector3d const gridSpacePos = (worldSpacePos - origin) / cellSize;
	return {(int)gridSpacePos[0], (int)gridSpacePos[1], (int)gridSpacePos[2]};
}
int MacGrid::cellCenterIdx(int i, int j, int k) const{
	return (i * (cellNum_i * cellNum_j) + j * cellNum_i + k);
}
int MacGrid::cellCenterIdx(Vector3d const & worldSpacePos){
	auto [i, j, k] = gridCoord(worldSpacePos);
	return cellCenterIdx(i, j, k);
}
////Get cell velocity ref from position
//Vector3d& MacGrid::getCellCenterVel(Vector3d const &  worldSpacePos){
//	return cellCenterVel[cellCenterIdx(worldSpacePos)];
//}
//Get Min Faces index
int MacGrid::cellMinFaceIdx_u(int i, int j, int k) const{
	return (i * (cellFaceNum_i * cellFaceNum_j) + j * cellFaceNum_i + k);
}
int MacGrid::cellMinFaceIdx_v(int i, int j, int k) const{
	return (i * (cellFaceNum_i * cellFaceNum_j) + j * cellFaceNum_i + k);
}
int MacGrid::cellMinFaceIdx_w(int i, int j, int k) const{
	return (i * (cellFaceNum_i * cellFaceNum_j) + j * cellFaceNum_i + k);
}
//Get MAX Faces index
int MacGrid::cellMaxFaceIdx_u(int i, int j, int k) const{
	return ((i + 1) * (cellFaceNum_i * cellFaceNum_j) + j * cellFaceNum_i + k);
}
int MacGrid::cellMaxFaceIdx_v(int i, int j, int k) const{
	return (i * (cellFaceNum_i * cellFaceNum_j) + (j + 1) * cellFaceNum_i + k);									
}
int MacGrid::cellMaxFaceIdx_w(int i, int j, int k) const{
	return (i * (cellFaceNum_i * cellFaceNum_j) + j * cellFaceNum_i + (k + 1));
}
//Get Face indices as tuples
tuple3i MacGrid::cellMinFaceIdcs(int i, int j, int k) const{
	return {cellMinFaceIdx_u(i, j, k),
			cellMinFaceIdx_v(i, j, k),
			cellMinFaceIdx_w(i, j, k)};
}
tuple3i MacGrid::cellMaxFaceIdcs(int i, int j, int k) const{
	return {cellMaxFaceIdx_u(i, j, k),
			cellMaxFaceIdx_v(i, j, k),
			cellMaxFaceIdx_w(i, j, k)};
}
tuple6i MacGrid::cellFaceIdcs(int i, int j, int k) const{
	return cellMinFaceIdcs(i, j, k) + cellMaxFaceIdcs(i, j, k);
}
tuple6i MacGrid::cellFaceIdcs(Vector3d const & worldSpacePos) const{
	auto [i, j, k] = gridCoord(worldSpacePos);
	return cellFaceIdcs(i, j, k);
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
	return {(cellFaceVel_u[cellMinFaceIdx_u(i, j, k)] + cellFaceVel_u[cellMaxFaceIdx_u(i, j, k)]) * .5f,
			(cellFaceVel_v[cellMinFaceIdx_v(i, j, k)] + cellFaceVel_v[cellMaxFaceIdx_v(i, j, k)]) * .5f,
			(cellFaceVel_w[cellMinFaceIdx_w(i, j, k)] + cellFaceVel_w[cellMaxFaceIdx_w(i, j, k)]) * .5f};
}

}