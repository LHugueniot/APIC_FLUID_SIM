#include "ParticleCacheIO.h"

namespace pic{

//==============================================================================================

bool ParticleCacheO::Initialize(){
    std::cout<<outputFilePath<<std::endl;
    pointFile = std::ofstream(outputFilePath);
    if (!pointFile)
        return false;
    pointFile<<"#Particle in cell sim cache."<<std::endl;
    pointFile<<"fps "<<std::to_string(fps)<<std::endl;
    return true;
}

bool ParticleCacheO::writePositions( uint frame, std::vector<Vector3d> const & pointPositions) {

    pointFile<<"frame "<<std::to_string(frame)<<std::endl;
    pointFile<<"{"<<std::endl;
    for(uint i = 0 ; i< pointPositions.size() ; i++){
        auto & pPos = pointPositions[i];

        auto x = pPos[0];
        auto y = pPos[1];
        auto z = pPos[2];

        pointFile<<"p "<<std::to_string(x)<<" "<<
        std::to_string(y)<<" "<<
        std::to_string(z)<<" "<<
        std::endl;
    }
    pointFile<<"}"<<std::endl;
    return true;
    
}

//==============================================================================================

bool ParticleCacheI::Initialize(){
    pointFile = std::ifstream(_inputFilePath);
    if (!pointFile)
        return false;
    return true;
}


bool ParticleCacheI::readPositions(uint frame, std::vector<Vector3d> & pointPositions) {
    std::string line;
    std::smatch regexMatch;
    bool frameScopeFound = false;
    while(std::getline(pointFile, line))
        if (std::regex_search(line, regexMatch, FRAME_REGEX) 
                && std::stoi(regexMatch[1].str()) == frame) {
            frameScopeFound = true;
            break;
        }
    std::cout<<"frameScopeFound = "<<frameScopeFound<<std::endl;
    if(!frameScopeFound)return false;

    while (std::getline(pointFile, line)) {
        if (line.find("}") == 0)
            break;
        if (std::regex_search(line, regexMatch, POINT_REGEX))
            pointPositions.push_back(
                Vector3d(
                    std::stof(regexMatch[1].str()), 
                    std::stof(regexMatch[2].str()),
                    std::stof(regexMatch[3].str())
                    )
                );
    }

    return true;
}

/**
//bool PicAlembicO::writePoints( uint frame, uint particleNum,
//    std::vector<Vector3d> const & pointPositions, 
//    std::vector<Vector3d> const & pointVelocities,
//    std::vector<tuple3d> const & pointColours){
//
//    std::vector<V3f> positions;
//    std::vector<V3f> velocities;
//    std::vector<C3f> colours;
//    std::vector<Alembic::Util::uint64_t> ids;
//
//    pSchema.set( OPointsSchema::Sample(V3fArraySample(positions), UInt64ArraySample(ids)) );
//    velProp.set(velocities);
//    rgbProp.set(colours);
//
//    for(uint i = 0 ; i< particleNum ; i++){
//        auto & p = pointPositions[i];
//        auto & v = pointVelocities[i];
//        auto & c = pointColours[i];
//
//        std::cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<"\n";
//        positions.push_back(V3f(p[0], p[1], p[2]));
//        velocities.push_back(V3f(v[0], v[1], v[2]));
//        colours.push_back(C3f(c[0], c[1], c[2]));
//        ids.push_back(i);
//    }
//
//    std::cout<<"Wrote "<<particleNum<<" particles to frame: "<<frame<<std::endl;
//}

bool PicAlembicI::readPoints(index_t frame,
    std::vector<Vector3d> & pointPositions, 
    std::vector<Vector3d> & pointVelocities){

    IPointsSchema::Sample psamp;
    pSchema.get( psamp, frame );

    auto particleNum = psamp.getPositions()->size();
    std::cout<<"particleNum:\n"<<particleNum<<std::endl;

    pointPositions.resize(particleNum);

    auto & inPoses = (*(psamp.getPositions()));
    for(uint i = 0 ; i< particleNum ; i++){

        pointPositions[i] = {inPoses.get()[i][0], inPoses.get()[i][1], inPoses.get()[i][2] };
        std::cout<<inPoses[i]<<std::endl;
    }

    pointVelocities.resize(particleNum);
    auto & inVels = *psamp.getVelocities();
    for(uint i = 0 ; i< particleNum ; i++)
        pointVelocities[i] = {inVels[i][0], inVels[i][1], inVels[i][2] };
}

bool PicAlembicI::readPoints(index_t frame,
    std::vector<Vector3d> & pointPositions, 
    std::vector<Vector3d> & pointVelocities,
    std::vector<tuple3d> & pointColours){

    IPointsSchema::Sample psamp;
    pSchema.get( psamp, frame );

    size_t particleNum = psamp.getPositions()->size();

    pointPositions.resize(particleNum);
    auto & inPoses = *psamp.getPositions();
    for(uint i = 0 ; i< particleNum ; i++)
        pointPositions[i] = {inPoses[i][0], inPoses[i][1], inPoses[i][2] };

    pointVelocities.resize(particleNum);
    auto & inVels = *psamp.getVelocities();
    for(uint i = 0 ; i< particleNum ; i++)
        pointVelocities[i] = {inVels[i][0], inVels[i][1], inVels[i][2] };

    //pointVelocities.resize(particleNum);
    //auto & inVels = psamp.getVelocities();
    //for(int i = 0 ; i< particleNum ; i++)
    //    pointVelocities[i] = {inVels[i][0], inVels[i][1], inVels[i][2] };
}
**/
}

