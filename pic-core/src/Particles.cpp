#include "Particles.h"

namespace pic{

template<>
Vector3dRef AffineParticles::getAffine<X>(uint i){
	return Vector3dRef(&affineVec_x[i * 3]);
}
template<>
Vector3dRef AffineParticles::getAffine<Y>(uint i){
	return Vector3dRef(&affineVec_y[i * 3]);
}
template<>
Vector3dRef AffineParticles::getAffine<Z>(uint i){
	return Vector3dRef(&affineVec_z[i * 3]);
}

template<>
Vector3d AffineParticles::getAffine<X>(uint i) const{
		return Vector3d{affineVec_x[i * 3], affineVec_x[i * 3 + 1], affineVec_x[i * 3 + 2]};
}
template<>
Vector3d AffineParticles::getAffine<Y>(uint i) const{
		return Vector3d{affineVec_y[i * 3], affineVec_y[i * 3 + 1], affineVec_y[i * 3 + 2]};
}
template<>
Vector3d AffineParticles::getAffine<Z>(uint i) const{
		return Vector3d{affineVec_z[i * 3], affineVec_z[i * 3 + 1], affineVec_z[i * 3 + 2]};
}

}