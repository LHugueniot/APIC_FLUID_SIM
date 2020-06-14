#include "common.h"
#include "MacGrid.h"

using namespace pic;

TEST_CASE("Deduced return types"){

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

TEST_CASE("C++ general"){

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

TEST_CASE("MacGrid constructor", "[MacGrid][constructor]"){

	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 1.f);
    REQUIRE(my_grid.origin == Vector3d(0,0,0));
    REQUIRE(my_grid.cellSize == 1.f);
    REQUIRE(my_grid.cellNum_i == 10);
	REQUIRE(my_grid.cellNum_j == 10);
	REQUIRE(my_grid.cellNum_k == 10);
	REQUIRE(my_grid.cellFaceNum_i == 9);
	REQUIRE(my_grid.cellFaceNum_j == 9);
	REQUIRE(my_grid.cellFaceNum_k == 9);
}

TEST_CASE("MacGrid Position conversions", "[MacGrid][gridCoord][cellSpacePos][gridSpacePos][cellCenterPos]"){

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

TEST_CASE("MacGrid getCellCenterVel", "[MacGrid][getCellCenterVel]"){

	MacGrid my_grid(Vector3d(0,0,0), 10, 10, 10, 1.f);
	Vector3d pWorldSpacePos{1.5f, 1.5f, 1.5f}; 

	Vector3d pVel{10000, 3, 214908};
	DEBUG();
	auto [min_u, min_v, min_w, max_u, max_v, max_w] = my_grid.cellFaceIdcs(1,1,1);
	DEBUG();
	my_grid.cellFaceVel_u[min_u] = 10000.f;
	my_grid.cellFaceVel_u[max_u] = 10000.f;
	DEBUG();
	my_grid.cellFaceVel_v[min_v] = 2.f;
	my_grid.cellFaceVel_v[max_v] = 4.f;
DEBUG();

	my_grid.cellFaceVel_w[min_w] = 214900.f;
	my_grid.cellFaceVel_w[max_w] = 214916.f;

	//std::cout<<"cellFaceIdcs\n"<<min_u<<"\n"<<max_u<<"\n"<<min_v<<"\n"<<max_v<<"\n"<<min_w<<"\n"<<max_w<<"\n";

	REQUIRE3d(my_grid.cellCenterVel(pWorldSpacePos) , Vector3d(10000, 3, 214908));
}


TEST_CASE("MacGrid isBoundaryCell", "[MacGrid][isBoundaryCell][cellIsOutOfBounds][cellIsOnBounds]"){
	MacGrid my_grid(Vector3d(0,0,0), 5, 5, 5, 1.f);

	tp start;
	ns end;

	start = tnow();
	for ( size_t i = 1 ; i < my_grid.cellFaceNum_i; i++)
		for ( size_t j = 1 ; j < my_grid.cellFaceNum_j; j++)
			for ( size_t k = 1 ; k < my_grid.cellFaceNum_k; k++){
				auto cellIsInBounds = my_grid.cellIsInBounds(i,j,k);
				REQUIRE(cellIsInBounds == true);
			}
    end = elapsed(start);
	std::cout<<"isBoundaryCell took: "<<end<<std::endl;

	REQUIRE(my_grid.cellIsOutOfBounds(0, 0, 0) == false);
	REQUIRE(my_grid.cellIsOutOfBounds(1, 1, 1) == false);
	REQUIRE(my_grid.cellIsOutOfBounds(-1, 1, 1) == true);
	REQUIRE(my_grid.cellIsOutOfBounds(5, 5, 5) == true);
	REQUIRE(my_grid.cellIsOutOfBounds(6, 6, 6) == true);


	REQUIRE(my_grid.cellIsInBounds(0, 0, 0) == true);
	REQUIRE(my_grid.cellIsInBounds(1, 1, 1) == true);
	REQUIRE(my_grid.cellIsInBounds(3, 4, 2) == true);
	REQUIRE(my_grid.cellIsInBounds(-1, 1, 1) == false);
	REQUIRE(my_grid.cellIsInBounds(5, 5, 5) == false);
	REQUIRE(my_grid.cellIsInBounds(6, 6, 6) == false);
}

TEST_CASE("MacGrid is valid cell face", "[MacGrid][isValidFace]"){

	MacGrid my_grid(Vector3d(0,0,0), 5, 5, 5, 1.f);

	REQUIRE(my_grid.isValidFace<maxU>(0, 0, 0) == true);
	REQUIRE(my_grid.isValidFace<maxU>(3, 0, 0) == true);
	REQUIRE(my_grid.isValidFace<maxU>(3, 4, 0) == true);
	REQUIRE(my_grid.isValidFace<maxU>(3, 4, 4) == true);
	REQUIRE(my_grid.isValidFace<maxU>(3, 0, 4) == true);
	REQUIRE(my_grid.isValidFace<maxU>(1, 1, 1) == true);
	REQUIRE(my_grid.isValidFace<maxU>(0, 1, 1) == true);
	REQUIRE(my_grid.isValidFace<maxU>(4, 1, 1) == false);
	REQUIRE(my_grid.isValidFace<maxU>(4, 1, 1) == false);
	REQUIRE(my_grid.isValidFace<maxU>(4, 4, 4) == false);
	REQUIRE(my_grid.isValidFace<maxU>(3, 4, 4) == true);

	REQUIRE(my_grid.isValidFace<maxV>(0, 0, 0) == true);
	REQUIRE(my_grid.isValidFace<maxV>(1, 1, 1) == true);
	REQUIRE(my_grid.isValidFace<maxV>(1, 0, 1) == true);
	REQUIRE(my_grid.isValidFace<maxV>(1, 4, 1) == false);
	REQUIRE(my_grid.isValidFace<maxV>(4, 4, 4) == false);
	REQUIRE(my_grid.isValidFace<maxV>(4, 3, 4) == true);

	REQUIRE(my_grid.isValidFace<maxW>(0, 0, 0) == true);
	REQUIRE(my_grid.isValidFace<maxW>(1, 1, 1) == true);
	REQUIRE(my_grid.isValidFace<maxW>(1, 1, 0) == true);
	REQUIRE(my_grid.isValidFace<maxW>(1, 1, 4) == false);
	REQUIRE(my_grid.isValidFace<maxW>(4, 4, 4) == false);
	REQUIRE(my_grid.isValidFace<maxW>(4, 4, 3) == true);

	REQUIRE(my_grid.isValidFace<minU>(0, 0, 0) == false);
	REQUIRE(my_grid.isValidFace<minU>(0, 1, 1) == false);
	REQUIRE(my_grid.isValidFace<minU>(1, 1, 1) == true);
	REQUIRE(my_grid.isValidFace<minU>(4, 1, 1) == true);
	REQUIRE(my_grid.isValidFace<minU>(4, 4, 4) == true);
	REQUIRE(my_grid.isValidFace<minU>(5, 4, 4) == false);
	REQUIRE(my_grid.isValidFace<minU>(4, 5, 4) == false);

	REQUIRE(my_grid.isValidFace<minV>(0, 0, 0) == false);
	REQUIRE(my_grid.isValidFace<minV>(1, 0, 1) == false);
	REQUIRE(my_grid.isValidFace<minV>(1, 1, 1) == true);
	REQUIRE(my_grid.isValidFace<minV>(1, 4, 1) == true);
	REQUIRE(my_grid.isValidFace<minV>(4, 4, 4) == true);
	REQUIRE(my_grid.isValidFace<minV>(4, 5, 4) == false);
	REQUIRE(my_grid.isValidFace<minV>(4, 4, 5) == false);

	REQUIRE(my_grid.isValidFace<minW>(0, 0, 0) == false);
	REQUIRE(my_grid.isValidFace<minW>(1, 1, 0) == false);
	REQUIRE(my_grid.isValidFace<minW>(1, 1, 1) == true);
	REQUIRE(my_grid.isValidFace<minW>(1, 1, 4) == true);
	REQUIRE(my_grid.isValidFace<minW>(4, 4, 4) == true);
	REQUIRE(my_grid.isValidFace<minW>(4, 4, 5) == false);
	REQUIRE(my_grid.isValidFace<minW>(5, 4, 4) == false);
}


TEST_CASE("MacGrid CellFace indexing", "[MacGrid][cellFaceIdx<minU>][cellFaceIdx<minV>][cellFaceIdx<minW>]\
	[cellFaceIdx<maxU>][cellFaceIdx<maxV>][cellFaceIdx<maxW>][cellMinFaceIdcs][cellMaxFaceIdcs][cellFaceIdcs]"){


	for ( uint cellNum_i = 5 ; cellNum_i < 10; cellNum_i++)
		for ( uint cellNum_j = 5 ; cellNum_j < 10; cellNum_j++)
			for ( uint cellNum_k = 5 ; cellNum_k < 10; cellNum_k++){

				DEBUG_VAR(cellNum_i);
				DEBUG_VAR(cellNum_j);
				DEBUG_VAR(cellNum_k);

				MacGrid my_grid(Vector3d(0,0,0), cellNum_i, cellNum_j, cellNum_k, 1.f);
				std::vector<int> cellFaceAccessed_u(my_grid.cellFaceVel_u.size(), 0);
				std::vector<int> cellFaceAccessed_v(my_grid.cellFaceVel_v.size(), 0);
				std::vector<int> cellFaceAccessed_w(my_grid.cellFaceVel_w.size(), 0);
			
				DEBUG_VAR(cellFaceAccessed_u.size());
				DEBUG_VAR(cellFaceAccessed_v.size());
				DEBUG_VAR(cellFaceAccessed_w.size());
				//my_grid.printOutMaxFaceVels();
				for ( uint i = 0 ; i < my_grid.cellNum_i; i++)
					for ( uint j = 0 ; j < my_grid.cellNum_j; j++)
						for ( uint k = 0 ; k < my_grid.cellNum_k; k++){
							//DEBUG();
							if(my_grid.isValidFace<maxU>(i, j, k)){
								//my_grid.cellFaceVel<maxU>(i, j, k) = 1;
								auto maxidx_u = my_grid.cellFaceIdx<maxU>(i, j, k);
								cellFaceAccessed_u[maxidx_u] += 1;
								assert(cellFaceAccessed_u.at(maxidx_u) >= 1);
							}
							if(my_grid.isValidFace<maxV>(i, j, k)){
								//my_grid.cellFaceVel<maxV>(i, j, k) = 1;
								auto maxidx_v = my_grid.cellFaceIdx<maxV>(i, j, k);
								cellFaceAccessed_v[maxidx_v] += 1;
								assert(cellFaceAccessed_v.at(maxidx_v) >= 1);
							}
							if(my_grid.isValidFace<maxW>(i, j, k)){
								//my_grid.cellFaceVel<maxW>(i, j, k) = 1;
								auto maxidx_w = my_grid.cellFaceIdx<maxW>(i, j, k);
								cellFaceAccessed_w[maxidx_w] += 1;
								assert(cellFaceAccessed_w.at(maxidx_w) >= 1);
							}
						}
			
				DEBUG_VAR(cellNum_i);
				DEBUG_VAR(cellNum_j);
				DEBUG_VAR(cellNum_k);
				for(int i = 0 ; i < cellFaceAccessed_u.size() ; i++){
					if(cellFaceAccessed_u[i] != 1){
						DEBUG_VAR(cellFaceAccessed_u);
						DEBUG_VAR(i);
					}
					REQUIRE(cellFaceAccessed_u[i] == 1);
				}
				for(int i = 0 ; i < cellFaceAccessed_v.size() ; i++){
					if(cellFaceAccessed_v[i] != 1){
						DEBUG_VAR(cellFaceAccessed_v);
						DEBUG_VAR(i);
					}
					REQUIRE(cellFaceAccessed_v[i] == 1);
				}
				for(int i = 0 ; i < cellFaceAccessed_w.size() ; i++){
					if(cellFaceAccessed_w[i] != 1){
						DEBUG_VAR(cellFaceAccessed_w);
						DEBUG_VAR(i);
					}
					REQUIRE(cellFaceAccessed_w[i] == 1);

				}

			}

}

TEST_CASE("MacGrid CellFace template indexing benchmarking", "[MacGrid][cellFaceIdx]"){

	MacGrid my_grid(Vector3d(0,0,0), 5, 5, 5, 1.f);

	auto start = tnow();
	auto t_minidx_u = my_grid.cellFaceIdx<minU>(1, 1, 1);
	auto t_minidx_v = my_grid.cellFaceIdx<minV>(1, 1, 1);
	auto t_minidx_w = my_grid.cellFaceIdx<minW>(1, 1, 1);
	auto t_maxidx_u = my_grid.cellFaceIdx<maxU>(1, 1, 1);
	auto t_maxidx_v = my_grid.cellFaceIdx<maxV>(1, 1, 1);
	auto t_maxidx_w = my_grid.cellFaceIdx<maxW>(1, 1, 1);
    auto end = elapsed(start);
	std::cout<<"templated cellFaceIdx took:"<<end<<std::endl;
		
}

TEST_CASE("MacGrid CellFaceVelocity direct accessors.", "[MacGrid][cellFaceVel]"){

	MacGrid my_grid(Vector3d(0,0,0), 5, 5, 5, 1.f);

	my_grid.cellFaceVel_u [my_grid.cellFaceIdx<minU>(1, 1, 1)] = 11;
	my_grid.cellFaceVel_v [my_grid.cellFaceIdx<minV>(1, 1, 1)] = 12;
	my_grid.cellFaceVel_w [my_grid.cellFaceIdx<minW>(1, 1, 1)] = 13;
	my_grid.cellFaceVel_u [my_grid.cellFaceIdx<maxU>(1, 1, 1)] = 14;
	my_grid.cellFaceVel_v [my_grid.cellFaceIdx<maxV>(1, 1, 1)] = 15;
	my_grid.cellFaceVel_w [my_grid.cellFaceIdx<maxW>(1, 1, 1)] = 16;

	REQUIRE(my_grid.cellFaceVel_u[my_grid.cellFaceIdx<minU>(1, 1, 1)] == my_grid.cellFaceVel<minU>(1, 1, 1));
	REQUIRE(my_grid.cellFaceVel_v[my_grid.cellFaceIdx<minV>(1, 1, 1)] == my_grid.cellFaceVel<minV>(1, 1, 1));
	REQUIRE(my_grid.cellFaceVel_w[my_grid.cellFaceIdx<minW>(1, 1, 1)] == my_grid.cellFaceVel<minW>(1, 1, 1));
	REQUIRE(my_grid.cellFaceVel_u[my_grid.cellFaceIdx<maxU>(1, 1, 1)] == my_grid.cellFaceVel<maxU>(1, 1, 1));
	REQUIRE(my_grid.cellFaceVel_v[my_grid.cellFaceIdx<maxV>(1, 1, 1)] == my_grid.cellFaceVel<maxV>(1, 1, 1));
	REQUIRE(my_grid.cellFaceVel_w[my_grid.cellFaceIdx<maxW>(1, 1, 1)] == my_grid.cellFaceVel<maxW>(1, 1, 1));
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

	//my_grid.cellFacePos_u[my_grid.cellFaceIdx<minU>(0, 0, 0)] = computedMinFacePos_u;

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



TEST_CASE("MacGrid bound checking", "[MacGrid][boundChecking]") 
{
	MacGrid my_grid(Vector3d(0,0,0), 5, 6, 10, 1.f);

	REQUIRE(my_grid.cellIsOnBounds(0,0,0));
	REQUIRE(my_grid.cellIsOnBounds(0,5,0));
	REQUIRE(my_grid.cellIsOnBounds(4,5,9));

	REQUIRE(my_grid.cellIsOutOfBounds(-1, -1, -1));
	REQUIRE(my_grid.cellIsOutOfBounds(1, 1, 10));
	REQUIRE(!my_grid.cellIsOutOfBounds(2, 3, 4));

	REQUIRE(!my_grid.cellIsInBounds(-1, 0, 0));
	REQUIRE(!my_grid.cellIsInBounds(-1, -1, -1));
	REQUIRE(my_grid.cellIsInBounds(1, 1, 5));
	REQUIRE(my_grid.cellIsInBounds(2, 3, 4));
	REQUIRE(my_grid.cellIsInBounds(4, 5, 9));

	REQUIRE(my_grid.isCellInGrid_i(0));
	REQUIRE(!my_grid.isCellInGrid_i(5));

	REQUIRE(my_grid.isCellInGrid_j(0));
	REQUIRE(my_grid.isCellInGrid_j(2));
	REQUIRE(!my_grid.isCellInGrid_j(6));
	REQUIRE(!my_grid.isCellInGrid_j(10));

	REQUIRE(my_grid.isCellInGrid_k(0));
	REQUIRE(my_grid.isCellInGrid_k(1));
	REQUIRE(!my_grid.isCellInGrid_k(-1));
	REQUIRE(!my_grid.isCellInGrid_k(10));

	REQUIRE(my_grid.isMaxBoundaryCell_i(4));
	REQUIRE(!my_grid.isMaxBoundaryCell_i(5));
	REQUIRE(!my_grid.isMaxBoundaryCell_i(2));

	REQUIRE(my_grid.isMinBoundaryCell_i(0));
	REQUIRE(!my_grid.isMinBoundaryCell_i(-1));
	REQUIRE(!my_grid.isMinBoundaryCell_i(2));

	REQUIRE(my_grid.isMaxBoundaryCell_j(5));
	REQUIRE(!my_grid.isMaxBoundaryCell_j(2));
	REQUIRE(!my_grid.isMaxBoundaryCell_j(6));

	REQUIRE(my_grid.isMinBoundaryCell_j(0));
	REQUIRE(!my_grid.isMinBoundaryCell_j(-2));
	REQUIRE(!my_grid.isMinBoundaryCell_j(2));

	REQUIRE(my_grid.isMaxBoundaryCell_k(9));
	REQUIRE(!my_grid.isMaxBoundaryCell_k(10));
	REQUIRE(!my_grid.isMaxBoundaryCell_k(8));

	REQUIRE(my_grid.isMinBoundaryCell_k(0));
	REQUIRE(!my_grid.isMinBoundaryCell_k(5));
	REQUIRE(!my_grid.isMinBoundaryCell_k(-2));

}

TEST_CASE("MacGrid check if cell is solid", "[MacGrid][isSolidCell]"){

	MacGrid my_grid(Vector3d(0,0,0), 5, 6, 10, 1.f);

	my_grid.cellCenterState[my_grid.cellCenterIdx(0,0,0)] = SOLID;

	REQUIRE(!SOLID == true);
	DEBUG_VAR(my_grid.cellCenterState[my_grid.cellCenterIdx(0,0,0)]);
	REQUIRE(!my_grid.cellCenterState[my_grid.cellCenterIdx(0,0,0)] == true);
	REQUIRE(my_grid.isSolidCell(0,0,0));
}


TEST_CASE("MacGrid get cell state ref method.", "[MacGrid][getCellState]"){

	MacGrid my_grid(Vector3d(0,0,0), 5, 6, 10, 1.f);
	my_grid.cellCenterState[my_grid.cellCenterIdx(0,0,0)] = SOLID;
	REQUIRE(my_grid.getCellState(0,0,0) == SOLID);

}

TEST_CASE("MacGrid grid coordinate rounding error test", "[MacGrid][gridCoord][error]"){

	MacGrid my_grid(Vector3d(0,0,0), 5, 5, 5, 1.f);

	REQUIRE(my_grid.gridCoord(Vector3d(0.00999999, 0.26, 1)) == tuple3i{0,0,1});


	REQUIRE(my_grid.gridCoord(Vector3d(1,1,1)) == tuple3i{1,1,1});
	REQUIRE(my_grid.gridCoord(Vector3d(1.f,1.f,1.f)) == tuple3i{1,1,1});


	int i = 1;
	std::cout<<(double)(i + 1.f)<<std::endl;
}