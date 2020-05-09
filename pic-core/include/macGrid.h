#pragma once
#ifndef MAC_GRID_H
#define MAC_GRID_H

#include "mathUtils.h"

namespace pic{


enum FaceDim{
	U = 0,
	V = 1,
	W = 2
};

struct MacGrid{

	enum CellState
	{
		LIQUID, AIR, SOLID
	};

	MacGrid(Vector3d const & _origin, size_t _cellNum_i, size_t _cellNum_j,
		size_t _cellNum_k, double _cellSize) :
	origin(_origin),
	cellSize(_cellSize),
	cellNum_i(_cellNum_i),
	cellNum_j(_cellNum_j),
	cellNum_k(_cellNum_k),
	cellFaceNum_i(_cellNum_i + 1),
	cellFaceNum_j(_cellNum_j + 1),
	cellFaceNum_k(_cellNum_k + 1),
	cellCenterPressure(cellNum_i * cellNum_j * cellNum_k, 0),
	cellFaceVel_u(cellFaceNum_i * cellFaceNum_j * cellFaceNum_k, 0),
	cellFaceVel_v(cellFaceNum_i * cellFaceNum_j * cellFaceNum_k, 0),
	cellFaceVel_w(cellFaceNum_i * cellFaceNum_j * cellFaceNum_k, 0),
	cellFaceWeightSum_u(cellFaceNum_i * cellFaceNum_j * cellFaceNum_k, 0),
	cellFaceWeightSum_v(cellFaceNum_i * cellFaceNum_j * cellFaceNum_k, 0),
	cellFaceWeightSum_w(cellFaceNum_i * cellFaceNum_j * cellFaceNum_k, 0),
	cellState(cellNum_i * cellNum_j * cellNum_k, CellState::AIR){
	}

	//Get position relative to grid origin and grid scale
	Vector3d gridSpacePos(Vector3d const & worldSpacePos) const;

	//Get position relative to cell origin and grid scale
	Vector3d cellSpacePos(Vector3d const & worldSpacePos, int i, int j, int k) const;

	//Get center position of a cell in worldspace
	Vector3d cellCenterPos(int i, int j, int k) const;

	//Gets the coord of a cell that can be used in the form i, j, k
	tuple3i gridCoord(Vector3d const & worldSpacePos) const;

	//From a cell coord, get the index for cell center data vectors
	int cellCenterIdx(int i, int j, int k) const;

	//From a world space position, get the index for cell center data vectors
	int cellCenterIdx(Vector3d const & worldSpacePos);

	//Get cell velocity ref from position
	Vector3d& getCellCenterVel(Vector3d const &  worldSpacePos);

	//Get indices for individual 
	//Get Min Faces index
	int cellMinFaceIdx_u(int i, int j, int k) const;
	int cellMinFaceIdx_v(int i, int j, int k) const;
	int cellMinFaceIdx_w(int i, int j, int k) const;

	//Get MAX Faces index
	int cellMaxFaceIdx_u(int i, int j, int k) const;
	int cellMaxFaceIdx_v(int i, int j, int k) const;
	int cellMaxFaceIdx_w(int i, int j, int k) const;

	//template<FaceDim D>
	//int cellMaxFaceIdx(int i, int j, int k) const{
	//	switch(D){
	//		case U:
	//		return ((i + 1) * (cellFaceNum_i * cellFaceNum_j) + j * cellFaceNum_i + k);
	//		case V:
	//		return (i * (cellFaceNum_i * cellFaceNum_j) + (j + 1) * cellFaceNum_i + k);
	//		case W:
	//		return (i * (cellFaceNum_i * cellFaceNum_j) + j * cellFaceNum_i + (k + 1));
	//	}
	//}

	//Get Face indices as tuples
	tuple3i cellMinFaceIdcs(int i, int j, int k) const;
	tuple3i cellMaxFaceIdcs(int i, int j, int k) const;
	tuple6i cellFaceIdcs(int i, int j, int k) const;
	tuple6i cellFaceIdcs(Vector3d const & worldSpacePos) const;

	//No need to cache as get function involves 1 array retrieval, one function and one computation
	//These on the other hand involve one computation which is cpu parallelised anyway.

	//Get Min Faces Positions
	Vector3d cellMinFacePos_u(int i, int j, int k) const;
	Vector3d cellMinFacePos_v(int i, int j, int k) const;
	Vector3d cellMinFacePos_w(int i, int j, int k) const;
	//Get MAX Faces Positions
	Vector3d cellMaxFacePos_u(int i, int j, int k) const;
	Vector3d cellMaxFacePos_v(int i, int j, int k) const;
	Vector3d cellMaxFacePos_w(int i, int j, int k) const;

	Vector3d cellCenterVel(int i, int j, int k) const;

	//Constants
	//Grid minimum corner
	Vector3d const origin;

	//For each cell the edges are all the same length
	double const cellSize;

	size_t const cellNum_i;
	size_t const cellNum_j;
	size_t const cellNum_k;

	size_t const cellFaceNum_i;
	size_t const cellFaceNum_j;
	size_t const cellFaceNum_k;

	//Velocities stored in the center of all cell faces
	std::vector<double> cellFaceVel_u;
	std::vector<double> cellFaceVel_v;
	std::vector<double> cellFaceVel_w;

	std::vector<double> cellFaceWeightSum_u;
	std::vector<double> cellFaceWeightSum_v;
	std::vector<double> cellFaceWeightSum_w;

	//Pressure stored at the center of each cell
	std::vector<double> cellCenterPressure;

	//Divergence stored at the center of each cell
	std::vector<double> cellCenterDivergence;

	//Store the type for all cells
	std::vector<CellState> cellState;
};

}

#endif /* MAC_GRID_H */