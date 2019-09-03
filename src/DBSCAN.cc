#include "DBSCAN.h"
#include "TMath.h"
#include <iostream>

int DBSCAN::run()
{
    int clusterID = 1;
    vector<Point>::iterator iter;
    for(iter = m_points.begin(); iter != m_points.end(); ++iter)
    {
        if ( iter->clusterID == UNCLASSIFIED )
        {
            if ( expandCluster(*iter, clusterID) != FAILURE )
            {
                clusterID += 1;
            }
        }
    }

    return (clusterID-1);
}
double DBSCAN::deltaPhi(double phi1, double phi2)
{
  double dphi = phi1-phi2;
  while (dphi > TMath::Pi())
  {
    dphi -= TMath::TwoPi();
  }
  while (dphi <= -TMath::Pi())
  {
    dphi += TMath::TwoPi();
  }
  return dphi;
};
int DBSCAN::result(int &nClusters, int cscLabels[], int clusterSize[], float clusterX[], float clusterY[], float clusterZ[], float clusterEta[], float clusterPhi[]){

  for (unsigned int i = 0;i < m_pointSize;i++)
  {
    cscLabels[i] = m_points[i].clusterID;
    // std::cout<<cscLabels.size()<<
    // std::cout<<"csclabel"<<cscLabels[i]<<std::endl;
  }
  for(int i = 0; i < nClusters; i++)
  {
    float avg_x(0.0), avg_y(0.0), avg_z(0.0);
    vector<Point>::iterator iter;
    for(iter = m_points.begin(); iter != m_points.end(); ++iter)
    {

      if ( iter->clusterID == i+1 )
      {
          clusterSize[i] ++;
          avg_x += iter->x;
          avg_y += iter->y;
          avg_z += iter->z;
      }
    }
    avg_x = avg_x/clusterSize[i];
    avg_y = avg_y/clusterSize[i];
    avg_z = avg_z/clusterSize[i];
    clusterX[i] = avg_x;
    clusterY[i] = avg_y;
    clusterZ[i] = avg_z;

    clusterPhi[i] = atan(avg_y/avg_x);
    if  (avg_x < 0.0){
      clusterPhi[i] = TMath::Pi() + clusterPhi[i];
    }
    clusterPhi[i] = deltaPhi(clusterPhi[i],0.0);
    double theta = atan(sqrt(pow(avg_x,2)+pow(avg_y,2))/abs(avg_z));
    clusterEta[i] = -1.0*TMath::Sign(1.0, avg_z)*log(tan(theta/2));
    // std::cout<<clusterPhi[i]<<", "<<clusterEta[i]<<", "<<clusterSize[i]<<std::endl;
  }
  return 0;
}


int DBSCAN::expandCluster(Point point, int clusterID)
{
    vector<int> clusterSeeds = calculateCluster(point);//neighbors, including itself

    if ( clusterSeeds.size() < m_minPoints )
    {
        point.clusterID = NOISE;
        return FAILURE;
    }
    else//core points
    {
        int index = 0, indexCorePoint = 0;
        vector<int>::iterator iterSeeds;
        //loop through neighbors of core points
        for( iterSeeds = clusterSeeds.begin(); iterSeeds != clusterSeeds.end(); ++iterSeeds)
        {
            m_points.at(*iterSeeds).clusterID = clusterID;//setting all the neighbors to the same clusterID
            //get the index of the core point itself
            if (m_points.at(*iterSeeds).x == point.x && m_points.at(*iterSeeds).y == point.y && m_points.at(*iterSeeds).z == point.z )
            {
                indexCorePoint = index;
            }
            ++index;
        }

        clusterSeeds.erase(clusterSeeds.begin()+indexCorePoint);
        //now clusterSeeds only contains neighbors, not itself; loop through the neighbors
        for( vector<int>::size_type i = 0, n = clusterSeeds.size(); i < n; ++i )
        {
            vector<int> clusterNeighors = calculateCluster(m_points.at(clusterSeeds[i]));
            //if this neighbor point is a core point
            if ( clusterNeighors.size() >= m_minPoints )
            {
                vector<int>::iterator iterNeighors;
                for ( iterNeighors = clusterNeighors.begin(); iterNeighors != clusterNeighors.end(); ++iterNeighors )
                {
                    if ( m_points.at(*iterNeighors).clusterID == UNCLASSIFIED || m_points.at(*iterNeighors).clusterID == NOISE )
                    {
                        if ( m_points.at(*iterNeighors).clusterID == UNCLASSIFIED )
                        {
                            clusterSeeds.push_back(*iterNeighors);
                            n = clusterSeeds.size();
                        }
                        m_points.at(*iterNeighors).clusterID = clusterID;
                    }
                }
            }
        }

        return SUCCESS;
    }
}

vector<int> DBSCAN::calculateCluster(Point point)
{
    int index = 0;
    vector<Point>::iterator iter;
    vector<int> clusterIndex;
    Point minimum;
    float min_dis = 1000000.0;
    for( iter = m_points.begin(); iter != m_points.end(); ++iter)
    {
        if (min_dis > calculateDistance(point, *iter)) minimum  = *iter;
        if ( calculateDistance(point, *iter) <= m_epsilon )
        {
            clusterIndex.push_back(index);
        }
        index++;
    }
    // if (index > 10) std::cout<<calculateDistance(point, minimum)<<std::endl;
    return clusterIndex;
}

inline double DBSCAN::calculateDistance( Point pointCore, Point pointTarget )
{
    return sqrt(pow(pointCore.x - pointTarget.x,2)+pow(pointCore.y - pointTarget.y,2)+pow(pointCore.z - pointTarget.z,2));
}
