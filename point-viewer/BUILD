cc_library(srcs = glob(["src/*.cpp"]),
	hdrs = glob(["include/*.h"]),
	strip_include_prefix = "include",
	name="pv_lib",
	deps=["@gl", "@glew", "@sdl2", "@pic_core"])

#print (glob(["src/*.cpp"]))

cc_binary(name="pv",
	deps=["//:pv_lib"])

#cc_library(
#    name = "pic-core",
#    srcs = ["some_lib.cc"],
#    hdrs = ["include/some_lib.h"],
#    copts = ["-Ilegacy/some_lib/include"],
#)
