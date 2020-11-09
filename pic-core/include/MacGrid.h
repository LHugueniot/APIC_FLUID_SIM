#pragma once
#ifndef MAC_GRID_H
#define MAC_GRID_H

#include "MathUtils.h"

namespace pic{


enum FaceDim{
	minU = 0,
	maxU = 1,
	minV = 2,
	maxV = 3,
	minW = 4,
	maxW = 5
};

enum CellState{
	SOLID = 0,
	LIQUID = 1,
	AIR = 2
};

struct MacGrid{

public:
	MacGrid(Vector3d const & _origin, uint _cellNum_i, uint _cellNum_j,
		uint _cellNum_k, double _cellSize) :
	origin(_origin),
	CBStart(origin[0] + _cellSize * 1.f,
			origin[1] + _cellSize * 1.f,
			origin[2] + _cellSize * 1.f),
	CBEnd(origin[0] + _cellSize * (double)(_cellNum_i - 1),
		  origin[1] + _cellSize * (double)(_cellNum_j - 1),
		  origin[2] + _cellSize * (double)(_cellNum_k - 1)),
	cellSize(_cellSize),
	invCellSize(1.f/_cellSize),
	cellNum_i(_cellNum_i),
	cellNum_j(_cellNum_j),
	cellNum_k(_cellNum_k),
	cellFaceNum_i(_cellNum_i - 1),
	cellFaceNum_j(_cellNum_j - 1),
	cellFaceNum_k(_cellNum_k - 1),
	cellFaceVel_u(cellFaceNum_i * cellNum_j * cellNum_k, 0),
	cellFaceVel_v(cellNum_i * cellFaceNum_j * cellNum_k, 0),
	cellFaceVel_w(cellNum_i * cellNum_j * cellFaceNum_k, 0),
	cellFaceWeightSum_u(cellFaceNum_i * cellNum_j * cellNum_k, 0),
	cellFaceWeightSum_v(cellNum_i * cellFaceNum_j * cellNum_k, 0),
	cellFaceWeightSum_w(cellNum_i * cellNum_j * cellFaceNum_k, 0),
	cellCenterState(cellNum_i * cellNum_j * cellNum_k, AIR),
	//laplacianSparseNBR(cellNum_i * cellNum_j * cellNum_k, cellNum_i * cellNum_j * cellNum_k),
	cellCenterPressure(cellNum_i * cellNum_j * cellNum_k),
	cellCenterDensity(cellNum_i * cellNum_j * cellNum_k)
	//cellCenterDivergence(cellNum_i * cellNum_j * cellNum_k, 1)
	{

		std::cout<<"Collision box Start = "<<CBStart<<std::endl;
		std::cout<<"Collision box End = "<<CBEnd<<std::endl;

		for (uint i = 0 ; i < cellNum_i ; i += (cellNum_i - 1))
			for (uint j = 0 ; j < cellNum_j ; j++)
				for (uint k = 0 ; k < cellNum_k ; k++)
					cellCenterState.at(cellCenterIdx(i, j, k)) = SOLID;
	
		for (uint i = 0 ; i < cellNum_i ; i ++)
			for (uint j = 0 ; j < cellNum_j ; j += (cellNum_j - 1))
				for (uint k = 0 ; k < cellNum_k ; k++)
					cellCenterState.at(cellCenterIdx(i, j, k)) = SOLID;
	
		for (uint i = 0 ; i < cellNum_i ; i ++)
			for (uint j = 0 ; j < cellNum_j ; j++)
				for (uint k = 0 ; k < cellNum_k ; k += (cellNum_k - 1))
					cellCenterState.at(cellCenterIdx(i, j, k)) = SOLID;
		std::cout<<"Created grid: "<<cellNum_i<<"x"<<cellNum_j<<"x"<<cellNum_k<<std::endl;
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

	//Check if grid coords are not boundary
	bool cellIsOnBounds(int i, int j, int k) const;
	bool cellIsOutOfBounds(int i, int j, int k) const;
	bool cellIsInBounds(int i, int j, int k) const;

	bool isCellInGrid_i(int i) const;
	bool isCellInGrid_j(int j) const;
	bool isCellInGrid_k(int k) const;

	bool isMaxBoundaryCell_i(int i) const;
	bool isMinBoundaryCell_i(int i) const;

	bool isMaxBoundaryCell_j(int j) const;
	bool isMinBoundaryCell_j(int j) const;

	bool isMaxBoundaryCell_k(int k) const;
	bool isMinBoundaryCell_k(int k) const;

	template<FaceDim D>
	int cellFaceIdx(int i, int j, int k) const;

	template<FaceDim D>
	double & cellFaceVel(int i, int j, int k);

	template<FaceDim D>
	bool isValidFace(int i, int j, int k) const;

	void printOutMaxFaceVels();

	void printOutMaxValidFaces();
	void printOutCellStates();

	//Get Face indices as tuples
	tuple3i cellMinFaceIdcs(int i, int j, int k) const;
	tuple3i cellMaxFaceIdcs(int i, int j, int k) const;
	tuple6i cellFaceIdcs(int i, int j, int k) const;
	tuple6i cellFaceIdcs(Vector3d const & worldSpacePos) const;

	template<FaceDim D>
	Vector3d cellFacePos(int i, int j, int k) const;

	//Get Min Faces Positions
	Vector3d cellMinFacePos_u(int i, int j, int k) const;
	Vector3d cellMinFacePos_v(int i, int j, int k) const;
	Vector3d cellMinFacePos_w(int i, int j, int k) const;
	//Get MAX Faces Positions
	Vector3d cellMaxFacePos_u(int i, int j, int k) const;
	Vector3d cellMaxFacePos_v(int i, int j, int k) const;
	Vector3d cellMaxFacePos_w(int i, int j, int k) const;

	//Get cell velocity ref from position
	Vector3d cellCenterVel(int i, int j, int k) const;
	Vector3d cellCenterVel(Vector3d const & worldSpacePos) const;

	CellState & getCellState(int i, int j, int k);
	//CellState getCellState(int i, int j, int k)const{
	//	return cellCenterState[cellCenterIdx(i, j, k)];
	//}

	bool isSolidCell(int i, int j, int k) const{
		return cellCenterState[cellCenterIdx(i, j, k)] == SOLID;
	}
	//Constants
	//Grid minimum corner
	Vector3d const origin;

	//Collision box start
	Vector3d const CBStart;

	//Collision box end
	Vector3d const CBEnd;
	//Grid top corner
	//Vector3d const end;

	//For each cell the edges are all the same length
	double const cellSize;
	double const invCellSize;

	uint const cellNum_i;
	uint const cellNum_j;
	uint const cellNum_k;

	uint const cellFaceNum_i;
	uint const cellFaceNum_j;
	uint const cellFaceNum_k;

	//Velocities stored in the center of all cell faces
	std::vector<double> cellFaceVel_u;
	std::vector<double> cellFaceVel_v;
	std::vector<double> cellFaceVel_w;

	std::vector<double> cellFaceWeightSum_u;
	std::vector<double> cellFaceWeightSum_v;
	std::vector<double> cellFaceWeightSum_w;


//Pressure Force solving stuff

	//Store the type for all cells
	std::vector<CellState> cellCenterState;

	//Sparse Laplacian matrix containing neighbor coeffs 
	//SparseMatrix<double> laplacianSparseNBR;

	//Pressure stored at the center of each cell
	//SparseVector<double> cellCenterPressure;
	std::vector<double> cellCenterPressure;
	std::vector<uint> cellCenterDensity;

	//Divergence stored at the center of each cell
	//SparseMatrix<double> cellCenterDivergence;

	//Pressure solver
	//CholmodSupernodalLLT<SparseMatrix<double>> * cg;
	//CholmodSupernodalLLT<SparseMatrix<double>> cg;
	//SimplicialLLT<SparseMatrix<double>, Lower> cg;
	//ConjugateGradient<SparseMatrix<double>, Eigen::Lower, Eigen::IncompleteCholesky<double>> * cg;
	
	ConjugateGradient<SparseMatrix<double>, Eigen::Lower, Eigen::IdentityPreconditioner> cg;
	//ConjugateGradient<SparseMatrix<double>, Eigen::Lower, Eigen::DiagonalPreconditioner<double>> cg;
	//ConjugateGradient<SparseMatrix<double>, Eigen::Lower, Eigen::IncompleteCholesky<double>> cg;

};

struct FlipMacGrid : MacGrid{

	FlipMacGrid(Vector3d const & _origin, 
				uint _cellNum_i, 
				uint _cellNum_j,
				uint _cellNum_k, 
				double _cellSize) : 
	MacGrid(_origin, 
			_cellNum_i, 
			_cellNum_j, 
			_cellNum_k, 
			_cellSize),
	cellFaceVelOld_u(cellFaceNum_i * cellNum_j * cellNum_k, 0),
	cellFaceVelOld_v(cellNum_i * cellFaceNum_j * cellNum_k, 0),
	cellFaceVelOld_w(cellNum_i * cellNum_j * cellFaceNum_k, 0)
	{
	}

	std::vector<double> cellFaceVelOld_u;
	std::vector<double> cellFaceVelOld_v;
	std::vector<double> cellFaceVelOld_w;

};

}

#endif /* MAC_GRID_H */