#ifndef CYLINDRICAL_POINT_H
#define CYLINDRICAL_POINT_H




//struct cylindrical point
class CylindricalPoint{
    public:
        CylindricalPoint(double r, double a, double h): radius(r), angle(a), height(h){}
        CylindricalPoint(){}

        CylindricalPoint operator+(const CylindricalPoint& otherPoint){
            CylindricalPoint sumPoint;
            sumPoint.radius = this->radius + otherPoint.radius;
            sumPoint.angle = this->angle + otherPoint.angle;
            sumPoint.height = this->height + otherPoint.height;
            return sumPoint;
        }

        CylindricalPoint operator/(const double& scala){
            //@TODO: handle /0
            CylindricalPoint rs;
            if( scala != 0.0 ){
                rs.radius = this->radius / scala;
                rs.angle = this->angle / scala;
                rs.height = this->height / scala; 
            }
            return rs;
        }


        double radius;
        double angle;
        double height;
        unsigned int segmentId;
};

struct CylindricalPointThetaOrder {
    bool operator() (CylindricalPoint p1, CylindricalPoint p2){
        return p1.angle < p2.angle;
    }
};

#endif
