glm_path = "glm-master"

#print ( glob(["glm-master/glm/**/*.hpp", "glm-master/glm/*.hpp"]))
#print("{}/**/*.hpp".format(glm_path))

#cc_library(name="glm-lib",
#	hdrs = glob(["glm-master/glm/**/*.hpp",
#				"glm-master/glm/**/*.h",
#				"glm-master/glm/**/*.inl",
#				"glm-master/glm/*.inl",
#				"glm-master/glm/*.hpp"]),
#	strip_include_prefix = glm_path)
#deps=["@gl", "@glew", "@sdl2", ":glm-lib"])

cc_library(srcs = glob(["src/*.cpp"]),
	hdrs = glob(["include/*.h"]),
	strip_include_prefix = "include",
	name="pv-lib",
	deps=["@gl", "@glew", "@sdl2"])

cc_binary(name="pv",
	deps=["//:pv-lib"])