#ifndef DBSCAN_H
#define DBSCAN_H

#include <vector>
#include <cmath>

#define UNCLASSIFIED -1
// #define CORE_POINT 1
// #define BORDER_POINT 2
#define NOISE -2
#define SUCCESS 0
#define FAILURE -3

using namespace std;

typedef struct Point_
{
    float x, y, z;  // X, Y, Z position
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


    int run();
    double deltaPhi(double phi1, double phi2);
    int result(int &nClusters, int cscLabels[], int clusterSize[], float clusterX[], float clusterY[], float clusterZ[], float clusterEta[], float clusterPhi[], float clusterRadius[]);
    int clusterMoments(int &nClusters, float clusterMajorAxis[], float clusterMinorAxis[], float clusterXSpread[], float clusterYSpread[], float clusterZSpread[], float clusterEtaSpread[], float clusterPhiSpread[],float clusterEtaPhiSpread[], float clusterX[], float clusterY[], float clusterZ[], float clusterEta[], float clusterPhi[], int clusterSize[]);
    int vertexing(vector<float> cscX,vector<float> cscY, vector<float> cscZ, vector<float> cscDirX, vector<float> cscDirY, vector<float> cscDirZ, float &clusterVertexR, float &clusterVertexZ, float &clusterVertexDis,float &clusterVertexChi2, int &clusterVertexN,int &clusterVertexN1cm, int &clusterVertexN5cm, int &clusterVertexN10cm, int &clusterVertexN15cm, int &clusterVertexN20cm);
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
