#include "common.h"
#include "macGrid.h"

using namespace pic;

TEST_CASE("Deduced return types")
{

	Vector3d a(0,0,0);
	std::cout<< TYPE_NAME(a[0]) << std::endl;
	std::cout<< TYPE_NAME((int)std::floor(a[0])) << std::endl;
	
}

//tuple6i operator+(tuple3i const & a, tuple3i const & b){
//	return {a[0], a[1], a[2], b[0], b[1], b[2]};
//}

tuple6i add(tuple3i const & a, tuple3i const & b){
	return {a[0], a[1], a[2], b[0], b[1], b[2]};
}

TEST_CASE("C++ general")
{
	tuple3i a{1, 2, 3};
	tuple3i b{4, 5, 6};
	//auto c = a + b;

	auto start = tnow();
    auto c1 = a + b;
    auto end = elapsed(start);
	std::cout<< TYPE_NAME(c1) << std::endl;
	std::cout<<"tuple3i + tuple3i took: "<<end<<std::endl;

	start = tnow();
	auto c2 = add(a,b);
	end = elapsed(start);
	std::cout<< TYPE_NAME(c2) << std::endl;
	std::cout<<"add(tuple3i, tuple3i) took:"<<end<<std::endl;

	REQUIRE(c1 == c2);

	std::cout<<c1[0]<<" "<<c2[0]<<std::endl;
	std::cout<<c1[1]<<" "<<c2[1]<<std::endl;
	std::cout<<c1[2]<<" "<<c2[2]<<std::endl;
	std::cout<<c1[3]<<" "<<c2[3]<<std::endl;
	std::cout<<c1[4]<<" "<<c2[4]<<std::endl;
	std::cout<<c1[5]<<" "<<c2[5]<<std::endl;
}

TEST_CASE("MacGrid constructor", "[MacGrid][constructor]") 
{
	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 1.f);
    REQUIRE(my_grid.origin == Vector3d(0,0,0));
    REQUIRE(my_grid.cellSize == 1.f);
    REQUIRE(my_grid.cellNum_i == 10);
	REQUIRE(my_grid.cellNum_j == 10);
	REQUIRE(my_grid.cellNum_k == 10);
	REQUIRE(my_grid.cellFaceNum_i == 10 + 1);
	REQUIRE(my_grid.cellFaceNum_j == 10 + 1);
	REQUIRE(my_grid.cellFaceNum_k == 10 + 1);
}

TEST_CASE("MacGrid Position conversions", "[MacGrid][gridCoord][cellSpacePos][gridSpacePos][cellCenterPos]")
{
	MacGrid my_grid(Vector3d(1,1,1), 5, 5, 5, 1.f);
	Vector3d particlePos{2.5f, 2.5f, 2.5f};
	REQUIRE(my_grid.gridCoord(particlePos) == tuple3i{1,1,1});
	REQUIRE3d(my_grid.cellSpacePos(particlePos, 1, 1, 1), 
		Vector3d(.5f, .5f, .5f));
	REQUIRE3d(my_grid.gridSpacePos(particlePos),
		Vector3d(1.5f, 1.5f, 1.5f));
	REQUIRE3d(my_grid.cellCenterPos(1, 1, 1),
		Vector3d(2.5f, 2.5f, 2.5f));

	MacGrid my_grid2(Vector3d(0,0,0), 5, 5, 5, 3.f);
	particlePos = {2.5f, 2.5f, 2.5f};
	REQUIRE(my_grid2.gridCoord(particlePos) == tuple3i{0,0,0});
	REQUIRE3d( my_grid2.cellSpacePos(particlePos, 0, 0, 0),
		Vector3d(5.f/6.f, 5.f/6.f, 5.f/6.f));
	REQUIRE3d(my_grid2.gridSpacePos(particlePos),
		Vector3d(5.f/6.f, 5.f/6.f, 5.f/6.f));
	REQUIRE3d(my_grid2.cellCenterPos(1, 1, 1), 
		Vector3d(4.5f, 4.5f, 4.5f));

	MacGrid my_grid3(Vector3d(0,0,0), 5, 5, 5, 1.f);
	particlePos = {3.2f, 3.2f, 3.2f};
	REQUIRE(my_grid3.gridCoord(particlePos) == tuple3i{3,3,3});
	REQUIRE3d(my_grid3.cellSpacePos(particlePos, 0, 0, 0),
		Vector3d(3.2f, 3.2f, 3.2f););
	REQUIRE3d(my_grid3.gridSpacePos(particlePos),
		Vector3d(3.2f, 3.2f, 3.2f));
	REQUIRE3d(my_grid3.cellCenterPos(1, 1, 1),
		Vector3d(1.5f, 1.5f, 1.5f));

	particlePos = {0.f, 3.2f, 3.2f};
	REQUIRE(my_grid3.gridCoord(particlePos) == tuple3i{0,3,3});
}

//TEST_CASE("MacGrid CellCenter indexing", "[MacGrid][cellCenterIdx]") 
//{
//	MacGrid my_grid(Vector3d(0,0,0), 5, 5, 5, 1.f);
//	std::vector<int> cellCenterAccessed(my_grid.cellCenterVel.size(), 0);
//	for ( size_t u = 0 ; u < my_grid.cellNum_i ; u++)
//		for ( size_t v = 0 ; v < my_grid.cellNum_j ; v++)
//			for ( size_t w = 0 ; w < my_grid.cellNum_k ; w++){
//				auto idx = my_grid.cellCenterIdx(u, v, w);
//				REQUIRE(my_grid.cellCenterIdx(Vector3d{u + .5f, v + .5f, w + .5f}) == idx);
//				cellCenterAccessed[idx] += 1;
//			}
//	for(auto accessNum : cellCenterAccessed) REQUIRE(accessNum == 1);
//
//}

//TEST_CASE("MacGrid getCellCenterVel", "[MacGrid][getCellCenterVel]")
//{
//	MacGrid my_grid(Vector3d(0,0,0), 1, 1, 1, 1.f);
//	my_grid.cellCenterVel[my_grid.cellCenterIdx(Vector3d{.5f, .5f, .5f})] = Vector3d{10000, 3, 214908};
//	Vector3d worldSpacePos{.75f, .75f, .75f}; 
//	REQUIRE3d(my_grid.getCellCenterVel(worldSpacePos) , Vector3d(10000, 3, 214908));
//}

TEST_CASE("MacGrid CellFace indexing", "[MacGrid][cellMinFaceIdx_u][cellMinFaceIdx_v][cellMinFaceIdx_w]\
	[cellMaxFaceIdx_u][cellMaxFaceIdx_v][cellMaxFaceIdx_w][cellMinFaceIdcs][cellMaxFaceIdcs][cellFaceIdcs]") 
{
	MacGrid my_grid(Vector3d(0,0,0), 5, 5, 5, 1.f);
	std::vector<int> cellFaceAccessed_u(my_grid.cellFaceVel_u.size(), 0);
	std::vector<int> cellFaceAccessed_v(my_grid.cellFaceVel_v.size(), 0);
	std::vector<int> cellFaceAccessed_w(my_grid.cellFaceVel_w.size(), 0);

	for ( size_t u = 0 ; u < my_grid.cellNum_i + 1; u++)
		for ( size_t v = 0 ; v < my_grid.cellNum_j + 1; v++)
			for ( size_t w = 0 ; w < my_grid.cellNum_k + 1; w++){
				auto minidx_u = my_grid.cellMinFaceIdx_u(u, v, w);
				auto minidx_v = my_grid.cellMinFaceIdx_v(u, v, w);
				auto minidx_w = my_grid.cellMinFaceIdx_w(u, v, w);

				REQUIRE(minidx_u < my_grid.cellFaceVel_u.size());
				REQUIRE(minidx_v < my_grid.cellFaceVel_v.size());
				REQUIRE(minidx_w < my_grid.cellFaceVel_w.size());
				cellFaceAccessed_u[minidx_u] += 1;
				cellFaceAccessed_v[minidx_v] += 1;
				cellFaceAccessed_w[minidx_w] += 1;

				auto maxidx_u = my_grid.cellMaxFaceIdx_u(u - 1, v, w);
				auto maxidx_v = my_grid.cellMaxFaceIdx_v(u, v - 1, w);
				auto maxidx_w = my_grid.cellMaxFaceIdx_w(u, v, w - 1);
				REQUIRE(minidx_u == maxidx_u);
				REQUIRE(minidx_v == maxidx_v);
				REQUIRE(minidx_w == maxidx_w);
			}
	for(auto accessNum : cellFaceAccessed_u) REQUIRE(accessNum == 1);
	for(auto accessNum : cellFaceAccessed_v) REQUIRE(accessNum == 1);
	for(auto accessNum : cellFaceAccessed_w) REQUIRE(accessNum == 1);
}

TEST_CASE("MacGrid computeCellFacePos vs getCellMinFacePos", "[MacGrid][minFacePos]")
{
	MacGrid my_grid(Vector3d(0,0,0), 5, 5, 5, 1.f);

	tp start;
	ns end;
	//Compute is faster for some reason (only one function call? and getCellMinFacePos 
	//still needs to do a computation as well as an array fetch which is parallelised upon compulation through optimisation)
	//old get function: my_grid.getCellMinFacePos_v(0, 0, 0);

	start = tnow();
	auto computedMinFacePos_u = my_grid.cellMinFacePos_u(0, 0, 0);
    end = elapsed(start);
	std::cout<<"cellMinFacePos_u took: "<<end<<std::endl;
	REQUIRE3d(computedMinFacePos_u, Vector3d(.0f, .5f, .5f));

	//my_grid.cellFacePos_u[my_grid.cellMinFaceIdx_u(0, 0, 0)] = computedMinFacePos_u;

	//start = tnow();
	//auto gotMinFacePos_u = my_grid.getCellMinFacePos_u(0, 0, 0);
    //end = elapsed(start);
	//std::cout<<"getCellMinFacePos_u took: "<<end<<std::endl;

	auto computedMinFacePos_v = my_grid.cellMinFacePos_v(0, 0, 0);
	REQUIRE3d(computedMinFacePos_v, Vector3d(.5f, .0f, .5f));
	
	auto computedMinFacePos_w = my_grid.cellMinFacePos_w(0, 0, 0);
	REQUIRE3d(computedMinFacePos_w, Vector3d(.5f, .5f, .0f));
	
	auto computedMaxFacePos_v = my_grid.cellMaxFacePos_v(0, 0, 0);
	REQUIRE3d(computedMaxFacePos_v, Vector3d(.5f, 1.f, .5f));
	
	auto computedMaxFacePos_w = my_grid.cellMaxFacePos_w(0, 0, 0);
	REQUIRE3d(computedMaxFacePos_w, Vector3d(.5f, .5f, 1.f));
}