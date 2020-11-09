#pragma once
#ifndef PIC_CORE_H
#define PIC_CORE_H

#include "MacGrid.h"
#include "Particles.h"
#include "Collisions.h"

namespace pic
{


//=====================================HELPER FUNCTIONS========================================

std::vector<double> AABBRandomParticles(Vector3d const & fieldC000, 
	Vector3d const & fieldC111, uint particleSystemSize);

std::vector<double> AABCubeUniformParticles(Vector3d const & fieldC000, 
	Vector3d const & fieldC111, double interParticleDistance);

tuple3i getClosestNBRCellIdcs(MacGrid const & grid, 
	Vector3d const & worldSpacePos, int i, int j, int k);

tuple6i getCellNBRIdcs(int i, int j, int k);

tuple6i getLocalNBRIdcs(MacGrid const & grid, 
	Vector3d const & worldSpacePos, int i, int j, int k);

//Get cell face Indices
template<FaceDim D>
tuple8i getStageredCellFaceNBRIdcs(MacGrid const & grid, 
	int min_i, int min_j, int min_k, 
	int max_i, int max_j, int max_k){

	return {grid.cellFaceIdx<D>(min_i, min_j, min_k), grid.cellFaceIdx<D>(max_i, min_j, min_k),
			grid.cellFaceIdx<D>(min_i, min_j, max_k), grid.cellFaceIdx<D>(max_i, min_j, max_k),
			grid.cellFaceIdx<D>(min_i, max_j, min_k), grid.cellFaceIdx<D>(max_i, max_j, min_k),
			grid.cellFaceIdx<D>(min_i, max_j, max_k), grid.cellFaceIdx<D>(max_i, max_j, max_k)};
}

tuple8i getStageredCellFaceNBRIdcs_u(MacGrid const & grid,
	int i, int min_j, int min_k, int max_j, int max_k);

tuple8i getStageredCellFaceNBRIdcs_v(MacGrid const & grid,
	int j, int min_i, int min_k, int max_i, int max_k);

tuple8i getStageredCellFaceNBRIdcs_w(MacGrid const & grid,
	int k, int min_i, int min_j, int max_i, int max_j);

void setDefaultCellStates(MacGrid & grid);

void setCellStates(MacGrid & grid, 
				   std::vector<tuple3i> & cellCoords, 
				   std::vector<CellState> & cellCenterStates);


//=============================================================================================
//=====================================FLUID SIM FUNCTIONS=====================================
//=============================================================================================

double calculateSubStep(MacGrid const & grid, double timeStep);

void enforceBoundaryVelocities(MacGrid & grid);
void enforceBoundaryPressure(MacGrid & grid);
void extrapolateBoundaryVelocities(MacGrid & grid);

void applyExternalForces(MacGrid & grid, double timeStep);

void applyPressureForces(MacGrid & grid, double timeStep);
void updateParticles(Particles & particles, double subStep);

//=============================================================================================
//=====================================PIC FUNCTIONS===========================================
//=============================================================================================

void transferAttributes(Particles const & particles, MacGrid & grid);
void transferAttributes(MacGrid const & grid, Particles & particles);
void advanceStep(Particles & particles, MacGrid & grid, double timeStep);

//=============================================================================================
//=====================================FLIP FUNCTIONS==========================================
//=============================================================================================

void transferAttributes(Particles const & particles, FlipMacGrid & grid);
void transferAttributes(FlipMacGrid const & grid, Particles & particles);
void advanceStep(Particles & particles, FlipMacGrid & grid, double timeStep);

//=============================================================================================
//=====================================APIC FUNCTIONS==========================================
//=============================================================================================

void transferAttributes(AffineParticles const & particles, MacGrid & grid);
void transferAttributes(MacGrid const & grid, AffineParticles & particles);
void advanceStep(AffineParticles & particles, MacGrid & grid, double timeStep);


//template<MacGrid::FaceDim F, typename T>
//void cellFaceAttrTransfer(Vector3d const & pos, T const & attr, 
//	tuple8i cellFaceNBRIdcs, vector<T> & gridAttrVector){
//
//}


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