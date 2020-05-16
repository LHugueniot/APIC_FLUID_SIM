cc_library(
	srcs = glob(["src/main.cpp"]),
	name = "test_main",
	hdrs = glob(["include/catch2.hpp", "include/common.h"]),
	strip_include_prefix = "include",
	alwayslink = True,
	deps = ["//:pic_core"]
	)

cc_test(
	name = "MathUtils",
	srcs = ["src/testMathUtils.cpp"],
	deps = [":test_main"])

cc_test( 
	name = "MacGrid",
	srcs = ["src/testMacGrid.cpp"],
	deps = [":test_main"])

cc_test(
	name = "PicCore",
	srcs = ["src/testPicCore.cpp"],
	deps = [":test_main"])

cc_test(
	name = "MacParticles",
	srcs = ["src/testMacParticles.cpp"],
	deps = [":test_main"])

cc_test(
	name = "PicStreamIO",
	srcs = ["src/testParticleCacheIO.cpp"],
	deps = [":test_main"])

test_suite(
	name = "all",
	tests = ["MathUtils",
			 "testMacGrid",
			 "testPicCore"])