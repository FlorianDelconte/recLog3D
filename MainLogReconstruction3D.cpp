#include <ctime>
#include <iostream>
#include <fstream>
#include <utility>
#include <stdlib.h>
#include <string>

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
    for(int i =0;i<c.size();i++){
        if(c[i][0]>lb[0] && c[i][1]>lb[1] && c[i][2]>lb[2] &&
           c[i][0]<ub[0] && c[i][1]<ub[1] && c[i][2]<ub[2] ){
            sub_c.push_back(c[i]);

        }
    }
    return sub_c;
}

/**
* to compute the bounding box from vector of Z3i point in CYL coordinate
* l : the vector of points
* TODO : template |CYL | XYZ |
*/
std::pair<CylindricalPoint ,CylindricalPoint>
computeBBCYL(std::vector<CylindricalPoint> l){
    // Initialisation des valeurs extrêmes.

    double minA = l[0].angle;
    double maxA = l[0].angle;
    double minR = l[0].radius;
    double maxR = l[0].radius;
    double minZ = l[0].height;
    double maxZ = l[0].height;
    // tous les points pour calculer les valeurs extrêmes.
    for (const CylindricalPoint& point : l) {
        double a = point.angle;
        double r = point.radius;
        double z = point.height;
        // Mise à jour des valeurs minimales et maximales.
        if (a < minA) minA = a;
        if (a > maxA) maxA = a;
        if (r < minR) minR = r;
        if (r > maxR) maxR = r;
        if (z < minZ) minZ = z;
        if (z > maxZ) maxZ = z;
    }
    std::pair<CylindricalPoint ,CylindricalPoint> p;
    p.first=CylindricalPoint(minR, minA, minZ);
    p.second=CylindricalPoint(maxR, maxA, maxZ);
    return p;
}

/**
* to compute the bounding box from vector of Z3i point in XYZ coordinate
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
/**
* get the scan number thats is the most present in the vector
* i_ps : the vector of pair : id point | id scan
*/
int
getBestScanId(std::vector<std::pair<int,int>> i_ps){
    int bestScan=-1;
    //4 scans max
    std::vector<int> counterScan=std::vector<int>(4,0);
    //count
    for(int i =0;i<i_ps.size();i++){
        counterScan[i_ps[i].second]+=1;
    }
    //max
    for(int i =0;i<counterScan.size();i++){
        if(counterScan[i]>bestScan){
            bestScan=i;
        }
    }
    return bestScan;
}

/**
*Reconstruct log by cylindrical sector loop. Search for each part of the sector, the best density (number of point in a part of cylindrical sector), delete all other in the part.
/* discretisationMap : the discret cells
/* rS : radius dimension (number of cells)
/* aS : angle dimension (number of cells)
/* zS : Height dimension (number of cells)
/* return a vector of index (points accepted) and the corresponding scan id
*/
std::vector<std::pair<int,int>>
reconstruct_by_BestScanInSecCyl(const std::vector<std::vector<std::vector<std::vector<std::pair<int,int>>>>> discretisationMap,
                        int rS,int aS,int zS){
    //log <indexPCL|indexScan> to reconstruct
    std::vector<std::pair<int,int>> logRec;
    for (int z = 0; z < zS; ++z){
        for (int a = 0; a < aS; ++a){
            for (int r = 0; r < rS; ++r){
                //vérifie que ce n'est pas une cellule vide
                if(discretisationMap[r][a][z].size()>0){
                    std::vector<std::pair<int,int>> id_PS=discretisationMap[r][a][z];
                    //calcul le scan le plus présent dans le secteur cylindrique
                    int scan_id=getBestScanId(id_PS);
                    //ajoute les points du bestScan a la reconstruction
                    for(int i = 0;i<id_PS.size();i++){
                        //si le numero du point est dans le bestScan
                        if(id_PS[i].second == scan_id){
                            //le point pcl
                            logRec.push_back(id_PS[i]);
                        }
                    }
                }
            }
        }
    }
    return logRec;
}

int
main(int argc,char **argv)
{

    //path to PCL log
    std::string logFile="../data/RGB_R_WSPS4_Cloud.xyz";
    //path to PCL centerline
    std::string centerlineFile="../data/WSPS4_centerline.xyz";
    //output prefixe
    std::string outputFile="WSPS4";
    //path to clean id log
    std::string logCleanFile="../data/WSPS4_clean.id";
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
    app.add_option("-i, --input", logFile , "log, xyz file");
    app.add_option("-c, --centerline", centerlineFile, "centerline, xyz file");
    app.add_option("-p, --preCleanedFile", logCleanFile, "cleaned log, id file");
    app.add_option("-r, --radiusCellsSize", rCellsS, "cell size of discretisation on radius axis");
    app.add_option("-a, --angleCellsSize", aCellsS, "cell size of discretisation on angle axis");
    app.add_option("-z, --zCellsSize", zCellsS, "cell size of discretisation on height axis");
    app.add_option("-n, --numberSec", nbSector, "number of sector for log reconstruction");
    app.add_option("-o, --output", outputFile, "output prefixe");
    app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);

    double load_duration=0.;
    double rec_duration=0.;
          /******/
          /*LOAD*/
          /******/
    trace.beginBlock("load ...");
    /***************/
    /* id clean log*/
    /***************/
    std::vector<unsigned int> cleanLogID=std::vector<unsigned int>();
    std::vector<unsigned int> cleanLogScanID=std::vector<unsigned int>();
    IOHelper::import(logCleanFile,cleanLogID,cleanLogScanID);
    /**********/
    /* pcl log*/
    /**********/
    std::vector<Z3i::RealPoint> logPcl = PointListReader<Z3i::RealPoint>::getPointsFromFile(logFile);
    /*************************/
    /* cenrerline and rescale*/
    /*************************/
    std::vector<Z3i::RealPoint> centerlinePcl = PointListReader<Z3i::RealPoint>::getPointsFromFile(centerlineFile);
    scaleCloud(centerlinePcl,SCALE);
    /********************************/
    /* one vector by sector by scans*/
    /********************************/
    std::vector< std::vector < std::vector< unsigned int > > > IdBSBS = toSector(logPcl,cleanLogID,cleanLogScanID,nbSector,nbScans);
    load_duration+=trace.endBlock();
    trace.info()<<"number of points in brut log : "<<logPcl.size()<<std::endl;
    trace.info()<<"number of points in clean log : "<<cleanLogID.size()<<std::endl;
    trace.info()<<"number of points in centerline : "<<centerlinePcl.size()<<std::endl;
    trace.info()<<"total time load : "<<load_duration<<"ms"<<std::endl;
          /*******************************/
          /*RECONSTRUCT PROCESS BY SECTOR*/
          /*******************************/
    trace.beginBlock("process ...");
    //colors (1 for each scans)
    Z3i::RealPoint palette[nbScans] {Z3i::RealPoint(255,0,0),Z3i::RealPoint(0,255,0),
                                Z3i::RealPoint(0,0,255),Z3i::RealPoint(255,255,0)};
    //global log reconstruct
    std::vector<unsigned int> globalRecId=std::vector<unsigned int>();
    std::vector<unsigned int> globalRecScan=std::vector<unsigned int>();
    //for each sector...
    for(int sec=0 ; sec < nbSector ; sec ++ ){
      trace.info()<<"------------"<<std::endl;
      trace.info()<<"extract local sector n°" <<sec<<" ...";
      //local  pcl point
      std::vector<Z3i::RealPoint> localXYZ = std::vector<Z3i::RealPoint> ();
      //local  id scan
      std::vector<unsigned int> localScan = std::vector<unsigned int> ();
      //local  id point
      std::vector<unsigned int> localId = std::vector<unsigned int> ();
      //pas très opti... redondant
      for (int sca = 0;sca<nbScans;sca++){
        for (int id = 0;id<IdBSBS[sec][sca].size();id++){
          localXYZ.push_back(logPcl[IdBSBS[sec][sca][id]]);
          localScan.push_back(sca);
          localId.push_back(IdBSBS[sec][sca][id]);
        }
      }
      trace.info()<<"OK"<<std::endl;
      //Boundary the centerline
      trace.info()<<"bound centerline...";
      std::vector<Z3i::RealPoint> subCenterlinePcl=constrainCenterline(centerlinePcl,computeBBXYZ(localXYZ));
      trace.info()<<"OK"<<std::endl;
      //Convert to cylindrical
      trace.info()<<"convert to cylindrical..."<<std::endl;;
      std::vector<CylindricalPoint> localCYL;
      //use sub centerline allow speedUp on cylindrical conversion (cause of lenght of centerline for each point)
      CylindricalCoordinateSystem ccs(subCenterlinePcl, Z3i::RealPoint(0.0,0.0,0.0));
      for(unsigned int i = 0; i < localXYZ.size(); i++){
          trace.progressBar(i, localXYZ.size());
          CylindricalPoint cylP = ccs.xyz2Cylindrical(localXYZ[i]);
          localCYL.push_back(cylP);
      }trace.info()<<std::endl;
      //Init discretisation map : each cells of discretisation is a pair of int : <index of point, scan number>
      trace.info()<<"init local discretisation...";
      std::pair<CylindricalPoint ,CylindricalPoint> bbcyl =computeBBCYL(localCYL);
      int rSize=floor((bbcyl.second.radius - bbcyl.first.radius)/rCellsS)+1;
      int aSize=floor((bbcyl.second.angle - bbcyl.first.angle)/aCellsS)+1;
      int zSize=floor((bbcyl.second.height - bbcyl.first.height)/zCellsS)+1;
      std::vector<std::vector<std::vector<std::vector<std::pair<int,int>>>>> discretisationMap(rSize);
      for (int r = 0; r < rSize; ++r){
          discretisationMap[r]=std::vector<std::vector<std::vector<std::pair<int,int>>>>(aSize);
          for (int a = 0; a < aSize; ++a){
              discretisationMap[r][a]=std::vector<std::vector<std::pair<int,int>>>(zSize);
              for (int z = 0; z < zSize; ++z){
                  discretisationMap[r][a][z]=std::vector<std::pair<int,int>>();
              }
          }
      }
      trace.info()<<"OK"<<std::endl;
      //insert id of point and scans in discretisation
      trace.info()<<"insert points in discretisation...";
      for(int i = 0; i < localCYL.size(); i++){
        CylindricalPoint mpCurrent=localCYL[i];
        int indr=floor((mpCurrent.radius- bbcyl.first.radius)/rCellsS);
        int inda=floor((mpCurrent.angle- bbcyl.first.angle)/aCellsS);
        int indz=floor((mpCurrent.height- bbcyl.first.height)/zCellsS);
        std::pair<int,int> id_PS;
        id_PS.first=localId[i];
        id_PS.second=localScan[i];
        discretisationMap[indr][inda][indz].push_back(id_PS);
      }
      trace.info()<<"OK"<<std::endl;
      //Use density by radius to reconstruct log
      trace.info()<<"reconstruction...";
      std::vector<std::pair<int,int>> rec=reconstruct_by_BestScanInSecCyl(discretisationMap,rSize,aSize,zSize);
      trace.info()<<"OK"<<std::endl;
      //add to global
      trace.info()<<"convert to global...";
      for (int i =0;i< rec.size();i++ ){
          globalRecId.push_back(rec[i].first);
          globalRecScan.push_back(rec[i].second);
      }
      trace.info()<<"OK"<<std::endl;
    }
    rec_duration+=trace.endBlock();
    trace.info()<<"number of points in reconstruted log : "<<globalRecId.size()<<std::endl;
    trace.info()<<"total time reconstruction : "<<rec_duration<<"ms"<<std::endl;
    /*****/
    /*Out*/
    /*****/
    //global log id
    IOHelper::export2Text(globalRecId,globalRecScan,outputFile+"_rec.id");
    //global log reconstruct with color per scan
    std::vector<Z3i::RealPoint> globalRecXYZ=std::vector<Z3i::RealPoint>();
    std::vector<Z3i::RealPoint> globalRecCOLOR=std::vector<Z3i::RealPoint>();

    for (int i =0;i< globalRecId.size();i++ ){
        globalRecXYZ.push_back(logPcl[globalRecId[i]]);
        globalRecCOLOR.push_back(palette[globalRecScan[i]]);
    }
    IOHelper::export2Text(globalRecXYZ,globalRecCOLOR,outputFile+"_rec_colorByScan.xyz");

}
