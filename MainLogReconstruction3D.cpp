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
    unsigned int nbP {20};
    //if log and centerline are not on the same scale, neeed to multiply centerline by :
    double SCALE=0.001;
    //number of scans
    int NBS = 4;
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
    app.add_option("-n, --numberSec", nbP, "number of sector for log analyze");
    app.add_option("-o, --output", outputFile, "output prefixe");
    app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);
    /*******************/
    /*LOAD ID clean log*/
    /*******************/
    std::vector<unsigned int> globalLogID=std::vector<unsigned int>();
    std::vector<unsigned int> globalLogScanID=std::vector<unsigned int>();
    IOHelper::import(logCleanId,globalLogID,globalLogScanID);
    /**************/
    /*LOAD PCL log*/
    /**************/
    trace.info()<<"read XYZ..."<<std::endl;
    std::vector<Z3i::RealPoint> logPcl = PointListReader<Z3i::RealPoint>::getPointsFromFile(logFile);

}
