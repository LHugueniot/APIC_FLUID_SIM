#include "common.h"
#include "MacParticles.h"

using namespace pic;

TEST_CASE("MacParticles constructor", "[MacParticles][constructor]") 
{
	std::vector<Vector3d> particlePos {Vector3d(.75f, .75f, .75f)};
	std::vector<Vector3d> particleVel {Vector3d(.75f, .75f, .75f)};
	MacParticles my_particles(1.f, particlePos, particleVel);

	REQUIRE(my_particles.size == 1);
	REQUIRE(Approx(my_particles.mass) == 1.f);
	REQUIRE(my_particles.positions == particlePos);
	REQUIRE(my_particles.velocities == particleVel);
}

TEST_CASE("MacParticles Accessor methods", "[MacParticles][getVec3Ref][getVec3Pos][setVec3Pos]") 
{
	std::vector<double> particlePos {.75f, .75f, .75f};
	std::vector<double> particleVel {.75f, .75f, .75f};
	MacParticles my_particles(1.f, particlePos, particleVel);


	auto start = tnow();
    auto & Vec3Ref = my_particles.getVec3Ref(0);
    Vec3Ref = Vector3d(5, 4, 6);
    auto end = elapsed(start);
	std::cout<< TYPE_NAME(Vec3Ref[0]) << std::endl;
	std::cout<<"my_particles.getVec3Ref[0] took: "<<end<<std::endl;
	REQUIRE(my_particles.positionsTest == std::vector<double>{5, 4, 6});

	//auto start = tnow();
    //auto Vec3Ref = my_particles.getVec3Ref[0];
    //auto end = elapsed(start);
	//std::cout<< TYPE_NAME(Vec3Ref) << std::endl;
	//std::cout<<"my_particles.getVec3Ref[0] took: "<<end<<std::endl;
//
	//my_particles
	//REQUIRE(my_particles.vel == particleVel);
}