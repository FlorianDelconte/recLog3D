#include <iostream>
#include <set>
#include <map>
#include <utility>
#include <fstream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include "IOHelper.h"

//using namespace DGtal;

//OK
void IOHelper::export2Text(const std::vector<unsigned int> &idPcl, const std::vector<unsigned int> &idSector, const std::string &filename){
  std::ofstream outStream;
  outStream.open(filename.c_str(), std::ofstream::out);
  for(unsigned int i = 0; i < idPcl.size();i++){
    unsigned int idpoint=idPcl[i];
    unsigned int idsector=idSector[i];
    outStream << idpoint << " "<< idsector<<std::endl;
  }
  outStream.close();
}

//OK
void IOHelper::export2Text(const std::vector<DGtal::Z3i::RealPoint> &v,std::vector<Z3i::RealPoint>c, const std::string &filename){
    std::ofstream outStream;
    outStream.open(filename.c_str(), std::ofstream::out);
    for(int i = 0; i < v.size();i++){

        Z3i::RealPoint p = v[i];
        Z3i::RealPoint col= c[i];
        outStream << p[0] << " "<< p[1] << " "<< p[2]<< " "
                  << col[0]<< " "<<col[1]<< " "<<col[2]<<std::endl;
    }
    outStream.close();
}

//OK
void IOHelper::export2Text(const std::vector<DGtal::Z3i::RealPoint> &v, const std::string &filename){
    std::ofstream outStream;
    outStream.open(filename.c_str(), std::ofstream::out);
    for(int i = 0; i < v.size();i++){
        Z3i::RealPoint p = v[i];
        outStream << p[0] << " "<< p[1] << " "<< p[2]<< std::endl;
    }
    outStream.close();
}

void IOHelper::import(const std::string &filePath,std::vector<unsigned int> &lCId,std::vector<unsigned int> &lCSId){
    std::ifstream f;

    f.open(filePath);

    std::string line, val;

    while (std::getline (f, line)) {
        std::vector<unsigned int> v;
        std::stringstream s (line);
        while (getline (s, val, ' '))
            v.push_back(std::stoi((val)));
        lCId.push_back(v[0]);
        lCSId.push_back(v[1]);
    }
}
