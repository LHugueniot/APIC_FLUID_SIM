cc_library(srcs = glob(["src/main.cpp"]),
	hdrs =["third-party/catch2.hpp"],
	strip_include_prefix="third-party",
	name="test-main",
	alwayslink = True
	)

cc_test(srcs = ["src/linear_interpolation.cpp",
				"src/testMathUtils.cpp",
				"src/testPicCore.cpp"],
	name="linearInterpolation",
	deps=[":test-main", "//:pic_core"])
