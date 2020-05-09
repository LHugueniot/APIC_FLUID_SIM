#include "common.h"
#include "macParticles.h"

using namespace pic;

TEST_CASE("MacParticles constructor", "[MacParticles][constructor]") 
{
	std::vector<Vector3d> particlePos {Vector3d(.75f, .75f, .75f)};
	std::vector<Vector3d> particleVel {Vector3d(.75f, .75f, .75f)};
	MacParticles my_particles(1.f, particlePos, particleVel);

	REQUIRE(my_particles.size == 1);
	REQUIRE(Approx(my_particles.mass) == 1.f);
	REQUIRE(my_particles.pos == particlePos);
	REQUIRE(my_particles.vel == particleVel);
}