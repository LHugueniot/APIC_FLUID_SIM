#include "common.h"

#include "MathUtils.h"

using namespace pic;

TEST_CASE("is nan", "[isnan]")
{
	REQUIRE(std::isnan(-0) == false);
}


TEST_CASE("int to uint", "[uint_to_int]"){
	unsigned int a = 1;
	int b = 0;
	unsigned int c = -1;
	REQUIRE(std::isnan(-0) == false);
}

TEST_CASE("gradBilinearInterpolation", "[math][gradientTrilinearInterpolation]"){
	tuple8d corners = {1.f, 1.f, -.5f, .125f,.125f, .125f, .125f, .125f};
	Vector3d position(0.5, 0.5, 0.5);

	auto result = gradientTrilinearInterpolation(corners, position);
	DEBUG_VAR(result);
}


TEST_CASE("getDiff")
{
	REQUIRE(Approx(getDiff(3, 0, 10)) == 0.3f);
	REQUIRE(Approx(getDiff(10, 5, 25)) == 0.25f);
	REQUIRE3d(getDiff(Vector3d{.5f,.5f,.5f}, Vector3d{0,0,0}, Vector3d{1,1,1}), Vector3d(.5f,.5f,.5f));
	REQUIRE3d(getDiff(Vector3d{.5f,.5f,.5f}, Vector3d{0,0,0}, Vector3d{2,2,2}), Vector3d(.25f,.25f,.25f));
	REQUIRE3d(getDiff(Vector3d{.5f,.5f,.5f}, Vector3d{0,0,0}, Vector3d{2,2,2}), Vector3d(.25f,.25f,.25f));
	REQUIRE3d(getDiff(Vector3d{.75f,.25f,.55f}, Vector3d{0,0,0}, Vector3d{1,1,1}), Vector3d(.75f,.25f,.55f));
}

TEST_CASE("linearInterpolation")
{
	REQUIRE(Approx(linearInterpolation(0.f, 0.f, 0.f)) == 0.f);
    REQUIRE(Approx(linearInterpolation(0.f, 1.f, .5f)) == .5f);
    REQUIRE(Approx(linearInterpolation(0.f, 2.f, .5f)) == 1.f);
    REQUIRE(Approx(linearInterpolation(0.f, 2.f, .3f)) == 0.6000000238f);
}

TEST_CASE("bilinearInterpolation")
{
    REQUIRE(bilinearInterpolation({1, 1, 1, 1}, .5f, .5f) == 1.0f);

    REQUIRE(bilinearInterpolation({.5f, .5f, .5f, .5f}, .5f, .5f) == .5f);
}

TEST_CASE("trilinearInterpolation") 
{
    REQUIRE(trilinearInterpolation({1, 1, 1, 1, 1, 1, 1, 1},
    								.5f, .5f, .5f) == 1.0f);

    REQUIRE(trilinearInterpolation({.5f, .5f, .5f, .5f,
    								.5f, .5f, .5f, .5f},
    								.5f, .5f, .5f) == .5f);
}

TEST_CASE("getWeights")
{
	//Returns as: w000, w100, w001, w101, w010, w110, w011, w111

	REQUIRE(getWeights(.5f, .5f, .5f) == tuple8d{.125f, .125f, .125f, .125f,.125f, .125f, .125f, .125f});
	REQUIRE(getWeights(Vector3d(.5f, .5f, .5f)) == tuple8d{.125f, .125f, .125f, .125f,.125f, .125f, .125f, .125f});
	REQUIRE(getWeights(Vector3d(.5f, .5f, .5f), Vector3d(0, 0, 0), Vector3d(1, 1, 1)) 
	== tuple8d{.125f, .125f, .125f, .125f,.125f, .125f, .125f, .125f});
	
	REQUIRE(getWeights(0.f, 0.f, 0.f) == tuple8d{1, 0, 0, 0, 0, 0, 0, 0});
	REQUIRE(getWeights(Vector3d(0.f, 0.f, 0.f)) == tuple8d{1, 0, 0, 0, 0, 0, 0, 0});
	REQUIRE(getWeights(Vector3d(0.f, 0.f, 0.f), Vector3d(0, 0, 0), Vector3d(1, 1, 1)) 
	== tuple8d{1, 0, 0, 0, 0, 0, 0, 0});
	
	REQUIRE(getWeights(1.f, 0.f, 0.f) == tuple8d{0, 1, 0, 0, 0, 0, 0, 0});
	REQUIRE(getWeights(Vector3d(1.f, 0.f, 0.f)) == tuple8d{0, 1, 0, 0, 0, 0, 0, 0});
	REQUIRE(getWeights(Vector3d(1.f, 0.f, 0.f), Vector3d(0, 0, 0), Vector3d(1, 1, 1)) 
	== tuple8d{0, 1, 0, 0, 0, 0, 0, 0});
	
	REQUIRE(getWeights(0.f, 1.f, 0.f) == tuple8d{0, 0, 0, 0, 1, 0, 0, 0});
	REQUIRE(getWeights(Vector3d(0.f, 1.f, 0.f)) == tuple8d{0, 0, 0, 0, 1, 0, 0, 0});
	REQUIRE(getWeights(Vector3d(0.f, 1.f, 0.f), Vector3d(0, 0, 0), Vector3d(1, 1, 1)) 
	== tuple8d{0, 0, 0, 0, 1, 0, 0, 0});
	
	REQUIRE(getWeights(0.f, 0.f, 1.f) == tuple8d{0, 0, 1, 0, 0, 0, 0, 0});
	REQUIRE(getWeights(Vector3d(0.f, 0.f, 1.f)) == tuple8d{0, 0, 1, 0, 0, 0, 0, 0});
	REQUIRE(getWeights(Vector3d(0.f, 0.f, 1.f), Vector3d(0, 0, 0), Vector3d(1, 1, 1)) 
	== tuple8d{0, 0, 1, 0, 0, 0, 0, 0});
	
	REQUIRE(getWeights(1.f, 1.f, 0.f) == tuple8d{0, 0, 0, 0, 0, 1, 0, 0});
	REQUIRE(getWeights(Vector3d(1.f, 1.f, 0.f)) == tuple8d{0, 0, 0, 0, 0, 1, 0, 0});
	REQUIRE(getWeights(Vector3d(1.f, 1.f, 0.f), Vector3d(0, 0, 0), Vector3d(1, 1, 1)) 
	== tuple8d{0, 0, 0, 0, 0, 1, 0, 0});
	
	REQUIRE(getWeights(1.f, 0.f, 1.f) == tuple8d{0, 0, 0, 1, 0, 0, 0, 0});
	REQUIRE(getWeights(Vector3d(1.f, 0.f, 1.f)) == tuple8d{0, 0, 0, 1, 0, 0, 0, 0});
	REQUIRE(getWeights(Vector3d(1.f, 0.f, 1.f), Vector3d(0, 0, 0), Vector3d(1, 1, 1)) 
	== tuple8d{0, 0, 0, 1, 0, 0, 0, 0});
	
	REQUIRE(getWeights(0.f, 1.f, 1.f) == tuple8d{0, 0, 0, 0, 0, 0, 1, 0});
	REQUIRE(getWeights(Vector3d(0.f, 1.f, 1.f)) == tuple8d{0, 0, 0, 0, 0, 0, 1, 0});
	REQUIRE(getWeights(Vector3d(0.f, 1.f, 1.f), Vector3d(0, 0, 0), Vector3d(1, 1, 1)) 
	== tuple8d{0, 0, 0, 0, 0, 0, 1, 0});
	
	REQUIRE(getWeights(1.f, 1.f, 1.f) == tuple8d{0, 0, 0, 0, 0, 0, 0, 1});
	REQUIRE(getWeights(Vector3d(1.f, 1.f, 1.f)) == tuple8d{0, 0, 0, 0, 0, 0, 0, 1});
	REQUIRE(getWeights(Vector3d(1.f, 1.f, 1.f), Vector3d(0, 0, 0), Vector3d(1, 1, 1)) 
	== tuple8d{0, 0, 0, 0, 0, 0, 0, 1});

	//c000, c100, c001, c101, c010, c110, c011, c111
	REQUIRE(getWeights(Vector3d(.25f, .25f, .25f), Vector3d(0, 0, 0), Vector3d(1, 1, 1)) 
	== tuple8d{0.421875f, 0.140625f, 0.140625f, 0.046875f, 0.140625f, 0.046875f, 0.046875f, 0.015625f});
}

TEST_CASE("getFromIdcs")
{
	std::vector<double> v{ .5f, .5f, .5f, .5f, 1.f, 1.f, 1.f, 1.f};
	tuple8i i {0, 1, 2, 3, 4, 5,6 ,7};
	auto fromIndices = getFromIdcs(v, i);

	std::cout<<TYPE_NAME(fromIndices[0])<<std::endl;
}

TEST_CASE("project")
{
	//REQUIRE(project(0, 10, 3) == 0.3f);
	//REQUIRE(project(5, 25, 10) == 0.25f);
}

TEST_CASE("intersectRaySphere")
{
	//REQUIRE(intersectRaySphere(0, 10, 3) == 0.3f);
	//REQUIRE(project(5, 25, 10) == 0.25f);
}

/**

**/