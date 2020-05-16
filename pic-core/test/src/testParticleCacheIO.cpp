#include "common.h"

#include <ParticleCacheIO.h>

using namespace pic;

TEST_CASE("Pic Input Output Regexes", "[POINT_REGEX][FRAME_REGEX][FPS_REGEX]")
{   
    std::string line("p 1.000000 1.000000 1.000000");
    std::smatch regexMatch;
    std::regex_search(line, regexMatch, POINT_REGEX);
    REQUIRE(std::stof(regexMatch[1].str()) == 1.f);
    REQUIRE(std::stof(regexMatch[2].str()) == 1.f);
    REQUIRE(std::stof(regexMatch[3].str()) == 1.f);

    line = "frame 1";
    std::regex_search(line, regexMatch, FRAME_REGEX);
    REQUIRE(std::stoi(regexMatch[1].str()) == 1);

}

TEST_CASE("Pic Input Output test", "[PicOuput][ParticleCacheI]")
{
    std::vector<Vector3d> particlesFrame1{
        Vector3d(1.543f,1,1),
        Vector3d(1,114,1.f),
        Vector3d(1,12,1),
        Vector3d(12431,1,1)
    };
    std::vector<Vector3d> particlesFrame2{
        Vector3d(1,1,1324),
        Vector3d(121,212,121),
        Vector3d(671,1,3241),
        Vector3d(57651,1,1456)
    };

    {
        ParticleCacheO a("/mnt/c/Users/lucciano/Desktop/dev/MAJOR2020/pic-core/test/frame_files/cache.pic", 24);
        a.writePositions(1, particlesFrame1);
        a.writePositions(2, particlesFrame2);
    }
    std::vector<Vector3d> readParticles1;
    std::vector<Vector3d> readParticles2;
    {
        ParticleCacheI in("/mnt/c/Users/lucciano/Desktop/dev/MAJOR2020/pic-core/test/frame_files/cache.pic");
        in.readPositions(1, readParticles1);
        in.readPositions(2, readParticles2);
    }
    REQUIRE(particlesFrame1 == readParticles1);
    REQUIRE(particlesFrame2 == readParticles2);
}

