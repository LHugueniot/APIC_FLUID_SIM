
A survey on different computational particle in cell methods.

Main loop psudocode:

void FluidSolver::update(float step) {
    projectVelocitiesToGrid();
    gravitySolve(step);
    enforceBoundary();
    pressureSolve(step);
    enforceBoundary();
    extrapolateVelocity();
    enforceBoundary();
    transferVelocitiesToParticles();
    updateParticlePositions(step);
    resolveCollisions();
    updateCells();

    frame++;
}

while(substep < timestep){

	//Find the maximum timestep such that particles dont advance further than one cell so 
	substep = calculateSubstep()

	//Transfer attributes from particles to grid
	transferParticleToGrid()

	//Change grid values based on
	advectFluid(substep){
		
		//This is solved by APIC?
		applyConvectionForces(substep)

		//Apply Gravity
		applyExternalForces(substep)

		//Maybe apply shearing forces
		viscocitySolve(substep)

		//Solve pressure equation
		applyPressureForces(substep){
			//Linear system we need to solve:
			// Lmat * Pvec = Dvec
			// Where Lmat is the laplacian matrix of neighbours
			// Pmat is the pressure vector
			// Dvec is the divergence vector

			//One half of the matrix equation
			calculateDivergenceMat()

			//Other half of the matrix equation
			calculateLaplacianMatrix()

			//Solve for pressure
			solvePressure()
		}
	}

	//Transfer attributes from grid to particles
	transferGridToParticle()

	particleKinematics()
}

Code snippets:

Two ways to solve

https://github.com/nitronoid/height_map/blob/master/src/main.cu

void apply_sor(
    cusp::csr_matrix<int, real, cusp::device_memory>::const_view di_A,
    cusp::array1d<real, cusp::device_memory>::const_view di_b,
    cusp::array1d<real, cusp::device_memory>::view do_x,
    const real i_w,
    const real i_tol,
    const int i_max_iter,
    const bool verbose)
{

cusp::relaxation::sor<real, cusp::device_memory> M(di_A, i_w);
M(di_A, di_b, do_x);
}


https://eigen.tuxfamily.org/dox/classEigen_1_1CholmodSupernodalLLT.html
https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html#TutorialSparseSolverConcept
//Init
CholmodSupernodalLLT<SparseMatrix<double>> cg;
cg.compute(A);

//Every tick
auto x = cg.solve(B);



Other main loop:
The algorithm that drives the fluid runs through these steps:
1. Calculate the simulation’s time step, t
2. Update the grid based on the presence of marker particles
3. Advance the velocity field, u
	3.1 Apply convection components using a backwards particle trace
	3.2 Apply external force components
	3.3 Apply the viscosity components
	3.4 Calculate the pressure to satisfy the conservation of mass equation
	3.5 Apply the pressure
	3.6 Extrapolate fluid velocities into buffer zone.
	3.7 Set solid cell velocities.
4. Move the particles through the velocity field, u for ttime.
4.1 Check if t extends beyond the next displayed frame time. Advance the particles to the next displayed frame. Display the frame and repeat step 4.1
4.2 Move the particles through the velocity field, u for the remainder of t.

Another main loop:
We initially set all particle positions & velocities
For each time step:
For each MAC grid cell we calculate the weighted average of neighbouring particle velocities
For FLIP: We store the grid velocities
Perform all the non-advection steps of a standard water simulator on the grid
For FLIP: Subtract the new grid velocities from the previously stored velocities, then add the interpolated difference to each particle’s velocity.
For PIC: Interpolate the new grid velocities to the previous particles’ velocities.
Mix Pic and Flip velocities
Move particles through the grid velocity field with an ODE solver, makingsure to push them outside of solid wall boundaries
Write the positions of the particles to an output.