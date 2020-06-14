#include "common.h"
#include "Particles.h"

using namespace pic;

TEST_CASE("Particles constructor", "[Particles][constructor]") 
{
	std::vector<double> particlePos {.75f, .75f, .75f};
	std::vector<double> particleVel {.75f, .75f, .75f};
	Particles my_particles(1.f, particlePos, particleVel);

	REQUIRE(my_particles.num == 1);
	REQUIRE(Approx(my_particles.mass) == 1.f);
	REQUIRE(my_particles.positions == particlePos);
	REQUIRE(my_particles.velocities == particleVel);
}

TEST_CASE("Particles Accessor methods", "[Particles][getVec3Ref][getVec3Pos][setVec3Pos]") 
{
	std::vector<double> particlePos {0, 0, 0};
	std::vector<double> particleVel {.75f, .75f, .75f};
	Particles my_particles(1.f, particlePos, particleVel);

	Vector3d toAdd(5, 4, 6);
	auto start = tnow();
    Vector3dRef pos1 = my_particles.getPos(0);
    pos1 += toAdd;
    auto end = elapsed(start);
	std::cout<< TYPE_NAME(pos1[0]) << std::endl;
	std::cout<<"Pos assignement by ref took: "<<end<<std::endl;
	REQUIRE(my_particles.positions == std::vector<double>{5, 4, 6});

	REQUIRE(pos1[0] == 5);
	REQUIRE(pos1[1] == 4);
	REQUIRE(pos1[2] == 6);
/**
	start = tnow();
    Vector3d pos2 = my_particles.getVec3Pos(0);
    my_particles.setVec3Pos(0, pos2 + toAdd);
    end = elapsed(start);
	std::cout<<"Pos assignement by get/set methods took: "<<end<<std::endl;
	REQUIRE(my_particles.positionsTest == std::vector<double>{10, 8, 12});

	start = tnow();
	Vector3d pos3 {
		my_particles.positionsTest[0],
		my_particles.positionsTest[1],
		my_particles.positionsTest[2]
	};
	pos3 += toAdd;
    my_particles.positionsTest[0] = pos3[0];
    my_particles.positionsTest[1] = pos3[1];
    my_particles.positionsTest[2] = pos3[2];
    end = elapsed(start);
	std::cout<<"Pos assignement by direct vec access took: "<<end<<std::endl;
	REQUIRE(my_particles.positionsTest == std::vector<double>{15, 12, 18});
**/

}