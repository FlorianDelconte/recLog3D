#include "CylindricalCoordinateSystem.h"

void CylindricalCoordinateSystem::init(){
    computeOxyNormals();
    computeThetaAxis();
}

CylindricalPoint CylindricalCoordinateSystem::xyz2Cylindrical(const Z3i::RealPoint &xyz){
    unsigned int segmentId = getSegment(xyz);
    Z3i::RealPoint aDir = getDirectionVector(segmentId);
    Z3i::RealPoint vectRadial = getRadialVector(xyz, aDir, centerline[segmentId]);
    double ra = vectRadial.norm();
    double angle = acos(vectRadial.dot(thetaAxis[segmentId])/
            thetaAxis[segmentId].norm()/vectRadial.norm());
    Z3i::RealPoint u = thetaAxis[segmentId].crossProduct(aDir);
    if (u.dot(vectRadial) < 0)
    {
        angle = 2 * M_PI - angle;
    }
    CylindricalPoint cyl(ra, angle, xyz[2]);


    cyl.segmentId = segmentId;

    return cyl;
}

Z3i::RealPoint
CylindricalCoordinateSystem::getRadialVector(const Z3i::RealPoint &aPoint, const Z3i::RealPoint &aDirection, const Z3i::RealPoint &p0){
    double dist = aDirection.dot(aPoint - p0);
    Z3i::RealPoint proj = p0 + dist*aDirection;
    return aPoint - proj;
}


unsigned int
CylindricalCoordinateSystem::getSegment(const Z3i::RealPoint &aPoint){
    unsigned int nbSegment = centerline.size() - 1;
    double lastSign = 1;
    for (int i = 1; i < nbSegment; i++){
        Z3i::RealPoint aVect = aPoint - centerline.at(i);
        double sign = aVect.dot(oxyNormals[i]);
        if(sign*lastSign <= 0){
            return i - 1;
        }
        lastSign = sign;
    }
    return nbSegment - 1;
}


Z3i::RealPoint
CylindricalCoordinateSystem::getDirectionVector(const unsigned int &segmentId){
    Z3i::RealPoint vectDir = centerline[segmentId + 1] - centerline[segmentId];
    return vectDir/vectDir.norm();
}
void
CylindricalCoordinateSystem::computeOxyNormals(){
    oxyNormals.resize(centerline.size() - 1);
    //std::vector<Z3i::RealPoint> ns(centerline.size() - 1);
    for (unsigned int i = 0; i < centerline.size() - 1; i++){
        Z3i::RealPoint vectDir = centerline.at(i + 1) - centerline.at(i);
        if(i == 0){
            oxyNormals[i] = vectDir.getNormalized();
        }else{
            Z3i::RealPoint previousVectDir = centerline.at(i) - centerline.at(i - 1);
            oxyNormals[i] = (previousVectDir + vectDir).getNormalized();
        }
    }
}

void
CylindricalCoordinateSystem::computeThetaAxis(){
    thetaAxis.resize(centerline.size() - 1);
    unsigned int seedSegmentId = getSegment(markPoint);;
    Z3i::RealPoint initDir = getDirectionVector(seedSegmentId);
    Z3i::RealPoint seedVector = getRadialVector(markPoint, initDir, centerline[seedSegmentId]);
    seedVector /= seedVector.norm();
    Z3i::RealPoint lastVectMark = seedVector;
    //normal vector of local plane O_ix_iy_i
    //dummy
    for (unsigned int i = seedSegmentId; i < centerline.size() - 1; i++){
        //Z3i::RealPoint Oi = centerline.at(i);
        //translation of lastVectorMark to O_i
        Z3i::RealPoint Ozi = getDirectionVector(i);
        //projection of trVtMark on Ozi
        double len = lastVectMark.dot(Ozi);
        if( i == seedSegmentId ){
        }
        //project trVtMark on Oixizi
        Z3i::RealPoint projvt = lastVectMark - len*Ozi;

        thetaAxis[i] = projvt.getNormalized();
        lastVectMark = thetaAxis[i];
    }
    if(seedSegmentId > 0){
        lastVectMark = seedVector;
        for (int i = seedSegmentId - 1; i >= 0 ; i--){
            //Z3i::RealPoint Ozi = Ozs[i];
            Z3i::RealPoint Ozi = getDirectionVector(i);
            //projection of trVtMark on Ozi
            double len = lastVectMark.dot(Ozi);
            //project trVtMark on Oixizi
            Z3i::RealPoint projvt = lastVectMark - len*Ozi;

            thetaAxis[i] = projvt.getNormalized();
            lastVectMark = thetaAxis[i];
        }
    }
}
