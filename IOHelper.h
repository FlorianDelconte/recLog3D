#ifndef IO_HELPER_H
#define IO_HELPER_H
#include <iostream>
#include <utility>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/writers/MeshWriter.h"
#include "DGtal/shapes/Mesh.h"

#include "CylindricalPoint.h"

using namespace DGtal;

typedef DGtal::PointVector<3U, double, std::array<double, 6UL>> LongRealPoint;
class IOHelper{

public:
    IOHelper(){
    }

    static void export2Text(const std::vector<DGtal::Z3i::RealPoint> &v, const std::string &filename);

    static void export2Text(const std::vector<DGtal::Z3i::RealPoint> &v,std::vector<Z3i::RealPoint> c, const std::string &filename);

    static void export2Text(const std::vector<unsigned int> &idPcl, const std::vector<unsigned int> &idSector, const std::string &filename);

    static void import(const std::string &filePath,std::vector<unsigned int> &lCId,std::vector<unsigned int> &lCSId);
};

#endif
