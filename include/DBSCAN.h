#ifndef DBSCAN_H
#define DBSCAN_H

#include <vector>
#include <cmath>
using namespace std;

#define UNCLASSIFIED -1
// #define CORE_POINT 1
// #define BORDER_POINT 2
#define NOISE -2
#define SUCCESS 0
#define FAILURE -3
// #define MERGE_DR
const double MERGE_CLUSTER_DR = 0.6;
const int N_phicorr = 8;
const int N_rcorr = 6;

const double r_corr[N_rcorr] = { 1.05, 1.1, 1.15, 1.2, 1.25, 1.3 };
const double phi_corr[N_phicorr] = { 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.80, 0.85 };

struct cscCluster
{
  float x, y, z, t, tTotal, tWire, tWirePruned, eta, phi;//t is t_strip, tWire, tTotal is sum
  int nCscSegments;
  float jetVeto, calojetVeto, muonVeto;
  int maxChamber, maxChamberSegment, nChamber;
  int maxStation, maxStationSegment, nStation, nStation5, nStation10, nStation10perc;
  float avgStation, avgStation5, avgStation10, avgStation10perc;
  int nCscSegmentChamberPlus11, nCscSegmentChamberPlus12, nCscSegmentChamberPlus13, nCscSegmentChamberPlus21, nCscSegmentChamberPlus22, nCscSegmentChamberPlus31, nCscSegmentChamberPlus32, nCscSegmentChamberPlus41, nCscSegmentChamberPlus42;
  int nCscSegmentChamberMinus11, nCscSegmentChamberMinus12, nCscSegmentChamberMinus13, nCscSegmentChamberMinus21, nCscSegmentChamberMinus22, nCscSegmentChamberMinus31, nCscSegmentChamberMinus32, nCscSegmentChamberMinus41, nCscSegmentChamberMinus42;

  float Me11Ratio, Me12Ratio;
  float MajorAxis, MinorAxis, EtaSpread, PhiSpread, EtaPhiSpread, XYSpread, DeltaRSpread;
  float XSpread, YSpread, ZSpread, TSpread, TWireSpread, TTotalSpread, TTotalSpreadPruned,RSpread;

  // float XSpread_phi0p5, YSpread_phi0p5, XYSpread_phi0p5, PhiSpread_phi0p5, EtaPhiSpread_phi0p5;
  // float XSpread_phi0p55, YSpread_phi0p55, XYSpread_phi0p55, PhiSpread_phi0p55, EtaPhiSpread_phi0p55;
  // float XSpread_phi0p6, YSpread_phi0p6, XYSpread_phi0p6, PhiSpread_phi0p6, EtaPhiSpread_phi0p6;
  // float XSpread_phi0p65, YSpread_phi0p65, XYSpread_phi0p65, PhiSpread_phi0p65, EtaPhiSpread_phi0p65;
  // float XSpread_phi0p7, YSpread_phi0p7, XYSpread_phi0p7, PhiSpread_phi0p7, EtaPhiSpread_phi0p7;
  // float XSpread_phi0p75, YSpread_phi0p75, XYSpread_phi0p75, PhiSpread_phi0p75, EtaPhiSpread_phi0p75;
  //
  // float XSpread_phi0p7_r1p3, YSpread_phi0p7_r1p3, XYSpread_phi0p7_r1p3, EtaSpread_phi0p7_r1p3, PhiSpread_phi0p7_r1p3, EtaPhiSpread_phi0p7_r1p3, RSpread_phi0p7_r1p3;
  // float XSpread_phi0p7_r1p2, YSpread_phi0p7_r1p2, XYSpread_phi0p7_r1p2, EtaSpread_phi0p7_r1p2, PhiSpread_phi0p7_r1p2, EtaPhiSpread_phi0p7_r1p2, RSpread_phi0p7_r1p2;
  // float XSpread_phi0p7_r1p25, YSpread_phi0p7_r1p25, XYSpread_phi0p7_r1p25, EtaSpread_phi0p7_r1p25, PhiSpread_phi0p7_r1p25, EtaPhiSpread_phi0p7_r1p25, RSpread_phi0p7_r1p25;
  // float XSpread_phi0p7_r1p15, YSpread_phi0p7_r1p15, XYSpread_phi0p7_r1p15, EtaSpread_phi0p7_r1p15, PhiSpread_phi0p7_r1p15, EtaPhiSpread_phi0p7_r1p15, RSpread_phi0p7_r1p15;
  // float XSpread_phi0p7_r1p1, YSpread_phi0p7_r1p1, XYSpread_phi0p7_r1p1, EtaSpread_phi0p7_r1p1, PhiSpread_phi0p7_r1p1, EtaPhiSpread_phi0p7_r1p1, RSpread_phi0p7_r1p1;
  //
  // float XSpread_r1p2, YSpread_r1p2, XYSpread_r1p2, EtaSpread_r1p2, EtaPhiSpread_r1p2, RSpread_r1p2, PhiSpread_r1p2;

  float XSpread_corr[N_phicorr][N_rcorr];
  float YSpread_corr[N_phicorr][N_rcorr];
  float XYSpread_corr[N_phicorr][N_rcorr];
  float EtaSpread_corr[N_phicorr][N_rcorr];
  float PhiSpread_corr[N_phicorr][N_rcorr];
  float EtaPhiSpread_corr[N_phicorr][N_rcorr];
  float RSpread_corr[N_phicorr][N_rcorr];

  float vertex_r, vertex_z, vertex_dis, vertex_chi2;
  int vertex_n, vertex_n1, vertex_n5, vertex_n10, vertex_n15, vertex_n20;
  vector<int>segment_id;
};

// struct largest_nCsc_cluster_
// {
//   inline bool operator() (const cscCluster& c1, const cscCluster& c2){return c1.nCscSegments > c2.nCscSegments;}
// } largest_nCsc_cluster;

typedef struct Point_
{
    float x, y, z, t, twire;  // X, Y, Z position
    float eta,phi;
    float dirX, dirY, dirZ;
    int station, chamber;
    int clusterID;  // clustered ID
}Point;

class DBSCAN {
public:
    DBSCAN(unsigned int minPts, float eps, vector<Point> points){
        m_minPoints = minPts;
        m_epsilon = eps;
        m_points = points;
        m_pointSize = points.size();
    }
    ~DBSCAN(){}
    int nClusters;
    vector<int>clusterSize;
    vector<int>cscLabels;
    vector<float>clusterEta;
    vector<float>clusterPhi;
    vector<float>clusterX;
    vector<float>clusterY;
    vector<float>clusterZ;
    vector<float>clusterTime;
    vector<float>clusterTimeWire;
    vector<float>clusterTimeWirePruned;
    vector<float>clusterTimeTotal;
    vector<float>clusterMajorAxis;
    vector<float>clusterMinorAxis;
    vector<float>clusterXSpread;
    vector<float>clusterYSpread;
    vector<float>clusterXYSpread;
    vector<float>clusterRSpread;

    vector<float>clusterZSpread;
    vector<float>clusterTimeSpread;
    vector<float>clusterTimeTotalSpread;
    vector<float>clusterTimeWireSpread;
    vector<float>clusterTimeTotalSpreadPruned;
    vector<float>clusterEtaPhiSpread;
    vector<float>clusterEtaSpread;
    vector<float>clusterPhiSpread;
    vector<float>clusterDeltaRSpread;




    // vector<float>clusterXSpread_phi0p5;
    // vector<float>clusterYSpread_phi0p5;
    // vector<float>clusterXYSpread_phi0p5;
    // vector<float>clusterPhiSpread_phi0p5;
    // vector<float>clusterEtaPhiSpread_phi0p5;
    //
    // vector<float>clusterXSpread_phi0p55;
    // vector<float>clusterYSpread_phi0p55;
    // vector<float>clusterXYSpread_phi0p55;
    // vector<float>clusterPhiSpread_phi0p55;
    // vector<float>clusterEtaPhiSpread_phi0p55;
    //
    // vector<float>clusterXSpread_phi0p6;
    // vector<float>clusterYSpread_phi0p6;
    // vector<float>clusterXYSpread_phi0p6;
    // vector<float>clusterPhiSpread_phi0p6;
    // vector<float>clusterEtaPhiSpread_phi0p6;
    //
    // vector<float>clusterXSpread_phi0p65;
    // vector<float>clusterYSpread_phi0p65;
    // vector<float>clusterXYSpread_phi0p65;
    // vector<float>clusterPhiSpread_phi0p65;
    // vector<float>clusterEtaPhiSpread_phi0p65;
    //
    // vector<float>clusterXSpread_phi0p7;
    // vector<float>clusterYSpread_phi0p7;
    // vector<float>clusterXYSpread_phi0p7;
    // vector<float>clusterPhiSpread_phi0p7;
    // vector<float>clusterEtaPhiSpread_phi0p7;
    //
    // vector<float>clusterXSpread_phi0p75;
    // vector<float>clusterYSpread_phi0p75;
    // vector<float>clusterXYSpread_phi0p75;
    // vector<float>clusterPhiSpread_phi0p75;
    // vector<float>clusterEtaPhiSpread_phi0p75;
    //
    // vector<float>clusterXSpread_phi0p7_r1p3;
    // vector<float>clusterYSpread_phi0p7_r1p3;
    // vector<float>clusterXYSpread_phi0p7_r1p3;
    // vector<float>clusterEtaSpread_phi0p7_r1p3;
    // vector<float>clusterPhiSpread_phi0p7_r1p3;
    // vector<float>clusterEtaPhiSpread_phi0p7_r1p3;
    // vector<float>clusterRSpread_phi0p7_r1p3;
    //
    // vector<float>clusterXSpread_phi0p7_r1p2;
    // vector<float>clusterYSpread_phi0p7_r1p2;
    // vector<float>clusterXYSpread_phi0p7_r1p2;
    // vector<float>clusterEtaSpread_phi0p7_r1p2;
    // vector<float>clusterPhiSpread_phi0p7_r1p2;
    // vector<float>clusterEtaPhiSpread_phi0p7_r1p2;
    // vector<float>clusterRSpread_phi0p7_r1p2;
    //
    // vector<float>clusterXSpread_phi0p7_r1p1;
    // vector<float>clusterYSpread_phi0p7_r1p1;
    // vector<float>clusterXYSpread_phi0p7_r1p1;
    // vector<float>clusterEtaSpread_phi0p7_r1p1;
    // vector<float>clusterPhiSpread_phi0p7_r1p1;
    // vector<float>clusterEtaPhiSpread_phi0p7_r1p1;
    // vector<float>clusterRSpread_phi0p7_r1p1;
    //
    // vector<float>clusterXSpread_phi0p7_r1p15;
    // vector<float>clusterYSpread_phi0p7_r1p15;
    // vector<float>clusterXYSpread_phi0p7_r1p15;
    // vector<float>clusterEtaSpread_phi0p7_r1p15;
    // vector<float>clusterPhiSpread_phi0p7_r1p15;
    // vector<float>clusterEtaPhiSpread_phi0p7_r1p15;
    // vector<float>clusterRSpread_phi0p7_r1p15;
    //
    // vector<float>clusterXSpread_phi0p7_r1p25;
    // vector<float>clusterYSpread_phi0p7_r1p25;
    // vector<float>clusterXYSpread_phi0p7_r1p25;
    // vector<float>clusterEtaSpread_phi0p7_r1p25;
    // vector<float>clusterPhiSpread_phi0p7_r1p25;
    // vector<float>clusterEtaPhiSpread_phi0p7_r1p25;
    // vector<float>clusterRSpread_phi0p7_r1p25;
    //
    // vector<float>clusterXSpread_r1p2;
    // vector<float>clusterYSpread_r1p2;
    // vector<float>clusterXYSpread_r1p2;
    // vector<float>clusterEtaSpread_r1p2;
    // vector<float>clusterEtaPhiSpread_r1p2;
    // vector<float>clusterRSpread_r1p2;
    // vector<float>clusterPhiSpread_r1p2;

    vector<vector<vector<float>>>clusterXSpread_corr;
    vector<vector<vector<float>>>clusterYSpread_corr;
    vector<vector<vector<float>>>clusterXYSpread_corr;
    vector<vector<vector<float>>>clusterEtaSpread_corr;
    vector<vector<vector<float>>>clusterPhiSpread_corr;
    vector<vector<vector<float>>>clusterEtaPhiSpread_corr;
    vector<vector<vector<float>>>clusterRSpread_corr;


    vector<float>clusterVertexR;
    vector<float>clusterVertexZ;
    vector<float>clusterVertexDis;
    vector<float>clusterVertexChi2;
    vector<int>clusterVertexN;

    vector<int>clusterVertexN1cm;
    vector<int>clusterVertexN5cm;
    vector<int>clusterVertexN10cm;
    vector<int>clusterVertexN15cm;
    vector<int>clusterVertexN20cm;


    vector<cscCluster> CscCluster;

    int run();
    double deltaPhi(double phi1, double phi2);
    double deltaR(double eta1, double phi1, double eta2, double phi2);

    int result();
    int clusterMoments();
    void clear_clusters();
    void sort_clusters();
    void merge_clusters();//change cluster id of the points if the clusters are close, and clear all the vectors at the end.
    int vertexing();
    vector<int> calculateCluster(Point point);
    int expandCluster(Point point, int clusterID);
    inline double calculateDistance(Point pointCore, Point pointTarget);

    int getTotalPointSize() {return m_pointSize;}
    int getMinimumClusterSize() {return m_minPoints;}
    int getEpsilonSize() {return m_epsilon;}
private:
    vector<Point> m_points;
    unsigned int m_pointSize;
    unsigned int m_minPoints;
    float m_epsilon;
};

#endif // DBSCAN_H
