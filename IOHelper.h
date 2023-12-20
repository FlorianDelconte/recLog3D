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

    template<typename T, typename = typename
            std::enable_if<std::is_arithmetic<T>::value>::type>
    static void export2Text(const std::vector<T> &xs, const std::vector<T> &ys, const std::string &filename){
        std::ofstream outStream;
        outStream.open(filename.c_str(), std::ofstream::out);
        if( xs.size() == ys.size() ){
            for(unsigned int i = 0; i < xs.size();i++){
                outStream << xs.at(i) << " "<< ys.at(i) <<std::endl;
            }
        }
        outStream.close();
    }

    static void export2Text(const std::vector<DGtal::Z3i::RealPoint> &v, const std::string &filename);

    static void export2Text(const std::vector<DGtal::Z3i::RealPoint> &v,std::vector<Z3i::RealPoint> c, const std::string &filename);

    static void export2Text(const std::vector<unsigned int> &idPcl, const std::vector<unsigned int> &idSector, const std::string &filename);

    static void import(const std::string &filePath,std::vector<unsigned int> &lCId,std::vector<unsigned int> &lCSId);
};

#endif
