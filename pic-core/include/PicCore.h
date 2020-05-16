#pragma once
#ifndef PIC_CORE_H
#define PIC_CORE_H

#include "MacGrid.h"
#include "MacParticles.h"

namespace pic
{

//template<MacGrid::FaceDim F, typename T>
//void cellFaceAttrTransfer(Vector3d const & pos, T const & attr, 
//	tuple8i cellFaceNBRIdcs, vector<T> & gridAttrVector){
//
//}

void generateRandParticlesInAABB(MacParticles & particles, int nParticles,
	double bCorner_x, double bCorner_y, double bCorner_z,	//top corner coords
	double tCorner_x, double tCorner_y, double tCorner_z);	//bottom corner coords

//void cellFaceVelTransfer_u(Vector3d const & pos, Vector3d const & vel, tuple8i const & cellUFaceNBRIdcs, MacGrid & grid);
//void cellFaceVelTransfer_v(Vector3d const & pos, Vector3d const & vel, tuple8i const & cellVFaceNBRIdcs, MacGrid & grid);
//void cellFaceVelTransfer_w(Vector3d const & pos, Vector3d const & vel, tuple8i const & cellWFaceNBRIdcs, MacGrid & grid);

tuple3i getClosestNBRCellIdcs(MacGrid const & grid, 
	Vector3d const & worldSpacePos, int i, int j, int k);

tuple6i getOrderedNBRIdcs(MacGrid const & grid, 
	Vector3d const & worldSpacePos, int i, int j, int k);

tuple8i getStageredCellFaceNBRIdcs_u(MacGrid const & grid,
	int min_i, int min_j, int min_k,
	int max_i, int max_j, int max_k);
tuple8i getStageredCellFaceNBRIdcs_v(MacGrid const & grid,
	int min_i, int min_j, int min_k,
	int max_i, int max_j, int max_k);
tuple8i getStageredCellFaceNBRIdcs_w(MacGrid const & grid,
	int min_i, int min_j, int min_k,
	int max_i, int max_j, int max_k);

void transferAttributes(MacParticles const & particles, MacGrid & grid);
void transferAttributes(MacGrid const & grid, MacParticles & particles);

double calculateSubStep(MacGrid const & grid, double timeStep);

void applyExternalForces(MacGrid & grid, double subStep);

void initializeLaplacianNBRMat(MacGrid & grid);
void initializeCellCenterDivergence(MacGrid & grid);
void applyPressureForces(MacGrid & grid, double subStep);

void advectParticles()
{

}

void pressureSolve()
{



}



//template<typename CollisionProc>
//void timeStep(ParticleAttributes & particleSim,  CollisionProc cp, double dt)
//{
//	int const nParticles = particleSim.velocities_x.size();
//
//	for (int i = 0 ; i < nParticles ; ++i)
//	{
//
//		Vector3d particlePos = cp({particleSim.positions_x[i],
//				particleSim.positions_y[i], particleSim.positions_z[i]},
//				{particleSim.velocities_x[i] * dt, 
//				(particleSim.velocities_y[i] + GRAV_Y) * dt,
//				particleSim.velocities_z[i] * dt});
//
//		particleSim.positions_x[i] = particlePos[0];
//		particleSim.positions_y[i] = particlePos[1];
//		particleSim.positions_z[i] = particlePos[2];
//	}
//}

}
#endif /* PIC_CORE_H */