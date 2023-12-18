#ifndef CYLINDRICAL_COORDINATE_SYSTEM_H
#define CYLINDRICAL_COORDINATE_SYSTEM_H

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include "CylindricalPoint.h"


using namespace DGtal;


class CylindricalCoordinateSystem{

// ----------------------- Standard methods ------------------------------
public:

    CylindricalCoordinateSystem(const std::vector<Z3i::RealPoint> &cen, const Z3i::RealPoint &seedPoint) : centerline(cen), markPoint(seedPoint){
        init();
    }

    CylindricalPoint xyz2Cylindrical(const Z3i::RealPoint &xyz);

    Z3i::RealPoint getDirectionVector(const unsigned int &segmentId);

    unsigned int getSegment(const Z3i::RealPoint &aPoint);

protected:
    void init();

    void computeThetaAxis();

    void computeOxyNormals();

    Z3i::RealPoint getRadialVector(const Z3i::RealPoint &aPoint, const Z3i::RealPoint &aDirection, const Z3i::RealPoint &p0);

    //////////////////////////////////////////
    //z coordinate
    std::vector<Z3i::RealPoint> centerline;
    //used to determine the theta
    Z3i::RealPoint markPoint;

    std::vector<Z3i::RealPoint> thetaAxis;
    std::vector<Z3i::RealPoint> oxyNormals;
};

#endif //end CENTERLINE_H
