#include "common.h"
#include "Collisions.h"

using namespace pic;

TEST_CASE("Box Collision", "[boxCollision]"){
	Vector3d c000 = {0, 0, 0};
	Vector3d c111 = {2, 1, 1};
	Vector3d pos = {.5f, .5f, .5f};
	Vector3d vel = {1, 1, 0};

	auto origPos = pos;
	auto origVel = vel;

	REQUIRE(boxCollisionBounce(c000, c111, pos, vel) == true);

	//REQUIRE(pos == Vector3d(1, 1, .5f));
	REQUIRE(pos == Vector3d(1, 1, .5f));
	REQUIRE(vel == Vector3d(1, -1, 0));
	DEBUG_VAR(pos);
}
