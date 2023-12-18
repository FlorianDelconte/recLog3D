#include <ctime>
#include <iostream>
#include <fstream>
#include <utility>
#include <stdlib.h>
#include <string>

#include <pcl/common/common.h>
#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/kdtree/kdtree.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/ModelCoefficients.h>
#include <pcl/point_types.h>
#include <pcl/features/normal_3d.h>
#include <pcl/segmentation/extract_clusters.h>

#include "DGtal/base/Common.h"
#include "DGtal/base/BasicTypes.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/writers/VolWriter.h"

#include "IOHelper.h"
#include "CylindricalCoordinateSystem.h"

#include "ext/CLI11.hpp"

using namespace DGtal;

/**
* extract the part of centerline thats is in the bb give in paramater
* c : the centerline, a vector of Z3i point
* bb : the bounding box, a pair of point (downRight,upLeft)
* Allows speed up for the cylindrical coordinate conversion
*/
std::vector<Z3i::RealPoint>
constrainCenterline(std::vector<Z3i::RealPoint> c,std::pair<Z3i::RealPoint ,Z3i::RealPoint>  bb){
    std::vector<Z3i::RealPoint>sub_c;
    Z3i::RealPoint lb=bb.first;
    Z3i::RealPoint ub=bb.second;
    trace.info()<<c.size()<<std::endl;
    for(int i =0;i<c.size();i++){
        if(c[i][0]>lb[0] && c[i][1]>lb[1] && c[i][2]>lb[2] &&
           c[i][0]<ub[0] && c[i][1]<ub[1] && c[i][2]<ub[2] ){
            sub_c.push_back(c[i]);

        }
    }
    return sub_c;
}

/**
* to compute the bounding box from vector of Z3i point
* l : the vector of points
*/
std::pair<Z3i::RealPoint ,Z3i::RealPoint>
computeBBXYZ(std::vector<Z3i::RealPoint> l){
    // init extrem values
    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();
    double minZ = std::numeric_limits<double>::max();
    double maxZ = std::numeric_limits<double>::lowest();
    for (const Z3i::RealPoint& point : l) {
        double x = point[0];
        double y = point[1];
        double z = point[2];
        //maj extrema values
        if (x < minX) minX = x;
        if (x > maxX) maxX = x;
        if (y < minY) minY = y;
        if (y > maxY) maxY = y;
        if (z < minZ) minZ = z;
        if (z > maxZ) maxZ = z;
    }
    //the bb is a pair (downRight, upLeft)
    std::pair<Z3i::RealPoint ,Z3i::RealPoint> p;
    p.first=Z3i::RealPoint(minX, minY, minZ);
    p.second=Z3i::RealPoint(maxX, maxY, maxZ);
    return p;
}

/**
* for rescale a vector of Z3i Point by scale paramater
* pts : the vector of Z3i point to rescal
* scale : the scale factor
*/
void
scaleCloud(std::vector<Z3i::RealPoint> &pts, double scale){
  for(unsigned int i = 0; i <  pts.size(); i++){
    Z3i::RealPoint &pi =  pts[i];
    pi *= scale;
  }
}
/**
* Distrbute each point in a vector (that's represent the sector).
* pcl : a vector of Z3i point
* id : a vector of unisgned int. The cleaned id point ( after noise processing)
* idS : a vector of unsigned int. The scan id point, number from 0 to 3 for each id. (That's represent the scan number)
* nbP : the number of sector (the lenght of sector is compute by height of the lof / by nbP)
* nbS : the number of scan. (3 for 'pin sylvestre')
*/
std::vector<std::vector<std::vector<unsigned int>>>
toSector(const std::vector<Z3i::RealPoint> &pcl, const std::vector<unsigned int> &id, const std::vector<unsigned int> &idS, int nbSec, int nbSca){
    // a vector of size : nbSec * nbSca * nbPoint
    std::vector<std::vector<std::vector<unsigned int> > > sectors(nbSec,std::vector<std::vector<unsigned int> >(nbSca));
    // compute the height of a sector
    std::pair<Z3i::RealPoint ,Z3i::RealPoint> bb=computeBBXYZ(pcl);
    Z3i::RealPoint lb=bb.first;
    Z3i::RealPoint ub=bb.second;
    double minZ = lb[2];
    double maxZ = ub[2];
    double rangeZ = maxZ - minZ;
    double sliceHeight = rangeZ / nbSec;
    // including input pcl for each parts
    for(int i=0;i<id.size();i++){
        auto fc = pcl[id[i]];
        unsigned int indexSec = floor((fc[2]-minZ)/sliceHeight);
        sectors[indexSec][idS[i]].push_back(id[i]);
    }
    return sectors;
}

int
main(int argc,char **argv)
{

    //path to PCL log
    std::string logFile="/volWork/these/DATA/THDATA/30_08_data_avricourt/FINAL/results/byInd/WSPS4/RGB_R_WSPS4_Cloud.xyz";
    //path to PCL centerline
    std::string centerlineFile="/volWork/these/DATA/THDATA/30_08_data_avricourt/FINAL/results/centerline/WSPS1.asc_centerline.xyz";
    //output prefixe
    std::string outputFile="WSPS4";
    //path to clean id log
    std::string logCleanId="/volWork/these/source/reconstruction/data/WSPS4_clean.id";
    // cell's size of discretisation on radius axis
    double rCellsS=0.05;//in meter
    // cell's size of discretisation on angle axis
    double aCellsS=(5*M_PI)/180;//in rad
    // cell's size of discretisation on heigth axis
    double zCellsS=0.05;//in meter
    //number of sector
    //(CARE : 2 SECTOR DELETED IN NOISE PROCESS -> first and last sector) not a good idea
    unsigned int nbSector {18};
    //if log and centerline are not on the same scale, neeed to multiply centerline by :
    double SCALE=0.001;
    //number of scans
    int nbScans = 3;
    /*******************************/
    /* parse command line using CLI*/
    /*******************************/
    CLI::App app;
    app.description("Allowed options are: ");
    app.add_option("-i, --input", logFile , "log xyz file");
    app.add_option("-c, --centerline", centerlineFile, "centerline xyz file");
    app.add_option("-r, --radiusCellsSize", rCellsS, "cell size of discretisation on radius axis");
    app.add_option("-a, --angleCellsSize", aCellsS, "cell size of discretisation on angle axis");
    app.add_option("-z, --zCellsSize", zCellsS, "cell size of discretisation on height axis");
    app.add_option("-n, --numberSec", nbSector, "number of sector for log analyze");
    app.add_option("-o, --output", outputFile, "output prefixe");
    app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);
    /*******************/
    /*LOAD ID clean log*/
    /*******************/
    std::vector<unsigned int> cleanLogID=std::vector<unsigned int>();
    std::vector<unsigned int> cleanLogScanID=std::vector<unsigned int>();
    IOHelper::import(logCleanId,cleanLogID,cleanLogScanID);
    /**************/
    /*LOAD PCL log*/
    /**************/
    trace.info()<<"read XYZ..."<<std::endl;
    std::vector<Z3i::RealPoint> logPcl = PointListReader<Z3i::RealPoint>::getPointsFromFile(logFile);
    /*******************************/
    /*LAUNCH CENTERLINE AND RESCALE*/
    /*******************************/
    std::vector<Z3i::RealPoint> centerlinePcl = PointListReader<Z3i::RealPoint>::getPointsFromFile(centerlineFile);
    scaleCloud(centerlinePcl,SCALE);
    /**************************************/
    /*CREATE ONE VECTOR BY SECTOR BY SCANS*/
    /**************************************/
    std::vector< std::vector < std::vector< unsigned int > > > IdBSBS = toSector(logPcl,cleanLogID,cleanLogScanID,nbSector,nbScans);
    /*******************************/
    /*RECONSTRUCT PROCESS BY SECTOR*/
    /*******************************/
    for(int sec=0 ; sec < nbSector ; sec ++ ){
      //The sector pcl point
      std::vector<Z3i::RealPoint> localXYZ = std::vector<Z3i::RealPoint> ();
      std::vector<unsigned int> localScan = std::vector<unsigned int> ();
      for (int sca = 0;sca<nbScans;sca++){
        for (int id = 0;id<IdBSBS[sec][sca].size();id++){
          localXYZ.push_back(logPcl[IdBSBS[sec][sca][id]]);
          localScan.push_back(sca);
        }
      }
      //Boundary the centerline
      std::vector<Z3i::RealPoint> subCenterlinePcl=constrainCenterline(centerlinePcl,computeBBXYZ(localXYZ));
      //Convert to cylindrical
      std::vector<CylindricalPoint> sectorCYL;
      //use sub centerline allow speedUp on cylindrical conversion (cause of lenght of centerline for each point)
      CylindricalCoordinateSystem ccs(subCenterlinePcl, Z3i::RealPoint(0.0,0.0,0.0));
      for(unsigned int i = 0; i < localXYZ.size(); i++){
          CylindricalPoint cylP = ccs.xyz2Cylindrical(localXYZ[i]);
          sectorCYL.push_back(cylP);
      }

    }
    /*****/
    /*Out*/
    /*****/
    /*
    std::vector<Z3i::RealPoint> logPclClean;
    for(int i =0;i<cleanLogID.size();i++){
      logPclClean.push_back(logPcl[cleanLogID[i]]);
    }
    IOHelper::export2Text(logPclClean,outputFile+"_clean.xyz");
    */
}
