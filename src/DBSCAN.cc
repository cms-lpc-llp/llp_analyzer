#include "DBSCAN.h"
#include "TMath.h"
#include <iostream>
#include "TVector3.h"
#include "TGraph.h"
#include "TF1.h"

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
int DBSCAN::result(int &nClusters, int cscLabels[], int clusterSize[], float clusterX[], float clusterY[], float clusterZ[], float clusterEta[], float clusterPhi[], float clusterRadius[]){

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

    // go through all the points again to get the farthest point from cluster
    float max_dis = 0.0;
    // vector<Point>::iterator iter;
    Point center;
    center.x = avg_x;
    center.y = avg_y;
    center.z = avg_z;
    center.clusterID =  i+1;
    for(iter = m_points.begin(); iter != m_points.end(); ++iter)
    {

      if ( iter->clusterID == i+1 )
      {
        float temp_dis = calculateDistance(*iter,center);
        if ( temp_dis > max_dis )max_dis = temp_dis;
      }
    }
    clusterRadius[i] = max_dis;

  }
  return 0;
}
//input x,y,z of segments in cluster, avgX, avgY, avgZ
int DBSCAN::vertexing(vector<float> cscX,vector<float> cscY, vector<float> cscZ, vector<float> cscDirX, vector<float> cscDirY, vector<float> cscDirZ, float &clusterVertexR, float &clusterVertexZ, int &clusterVertexN, float &clusterVertexDis)
{

  TVector3 vecDir;
  TVector3 vecCsc;
  TF1 *fit = new TF1();

  TGraph *gr = new TGraph(cscX.size());

  for(unsigned int i = 0;i < cscX.size();i++)
  {
    vecCsc.SetXYZ(cscX[i],cscY[i],cscZ[i]);
    vecDir.SetXYZ(cscDirX[i],cscDirY[i],cscDirZ[i]);
    double slope = (vecCsc.Perp() - (vecCsc+vecDir).Perp()) / (vecCsc.Z() - (vecCsc+vecDir).Z());
    double beta = -1.0*vecCsc.Z() * slope + vecCsc.Perp();
    gr->SetPoint(i,slope,beta);
  }
  // Here's where the fit happens
  // Process is repeated until a good vertex is found or fewer than 3 segments left
  // The segments form a line in RZ plane
  // Ideally, the vertex (z,r) is a point on every line in the cluster
  // So using all the measured segment slopes (a_i) and a vertex guess (z,r) gives predicted intercept (b_i)
  // r = a_i*z + b_i --> b_i = r - a_i*z
  // Fitting (a_i, b_i) points to "b_i = p[0] + p[1]*a_i" by minimizing chi-squared between measured andß
  // predicted intercepts gives best (r,z) estimate
  // Then calculate distance of closest approach of all segment lines and vertex
  // If max distance > 30cm throw out farthest point and refit until all segments within 30cm or < 3 segments left

  bool goodVtx = false;
  double distance = 0.0;
  double maxDistance = 0.0;
  double farSegment = -1;

  // cout << "starting with " << gr->GetN() << " segments" << endl;
  while (!goodVtx && gr->GetN()>=3){
  	maxDistance = 0.0;
  	gr->Fit("1++x","Q");
  	fit = gr->GetFunction("1++x");
  	for (int i=0; i<gr->GetN(); i++){
  	  // Distance from point (x_0, y_0) to line ax + by + c = 0 is |ax_0 + by_0 + c| / sqrt(a^2 + b^2)
  	  // Here: a = slope = "X"; b = -1,;c = intercept = "Y"; (x_0, y_0) = (z, r) = (-p[1], p[0]);
  	  distance = abs(gr->GetX()[i]*fit->GetParameter(1)*-1.0 + -1.0*fit->GetParameter(0) + gr->GetY()[i]);
  	  distance = distance / sqrt(pow(gr->GetX()[i],2)+pow(-1.0,2));
  	  if (distance > maxDistance){
  	    maxDistance = distance;
  	    farSegment = i;
  	  }
  	}
  	if (maxDistance < 30.0){
  	  goodVtx = true;
  	}
  	else {
  	  gr->RemovePoint(farSegment);
  	}
  }
  if (goodVtx && gr->GetN()>=3){
    clusterVertexR = fit->GetParameter(0);
    clusterVertexZ = -1.0*fit->GetParameter(1);
    clusterVertexN = gr->GetN();
    clusterVertexDis = maxDistance;
  }
  else{
    clusterVertexR = 0.0;
    clusterVertexZ = 0.0;
    clusterVertexN = 0;
    clusterVertexDis =-999;
  }
  // cout << "vertex? " << goodVtx << ", # segments = " << gr->GetN() << ", maxDistance = " << maxDistance << endl;
  // cout << "vertex R, Z, N:  " << clusterVertexR <<", " << clusterVertexZ << ",  " << clusterVertexN << endl;
  // cout << "maxDistance = " << maxDistance << endl;

  return 0;

}
int DBSCAN::clusterMoments(int &nClusters, float clusterMajorAxis[], float clusterMinorAxis[], float clusterXSpread[], float clusterYSpread[], float clusterZSpread[], float clusterEtaSpread[], float clusterPhiSpread[],float clusterEtaPhiSpread[], float clusterX[], float clusterY[], float clusterZ[], float clusterEta[], float clusterPhi[], int clusterSize[])
{


  for(int i = 0; i < nClusters; i++)
  {
    float m11(0.0), m12(0.0), m22(0.0);
    vector<Point>::iterator iter;
    for(iter = m_points.begin(); iter != m_points.end(); ++iter)
    {


      if ( iter->clusterID == i+1 )
      {
          float phi = atan(iter->y/iter->x);
          if  (iter->x < 0.0){
            phi = TMath::Pi() + phi;
          }
          phi = deltaPhi(phi,0.0);
          float eta = atan(sqrt(pow(iter->x,2)+pow(iter->y,2))/abs(iter->z));
          eta = -1.0*TMath::Sign(1.0, iter->z)*log(tan(eta/2));
          m11 += (eta-clusterEta[i])*(eta-clusterEta[i]);
          m12 += (eta-clusterEta[i])* deltaPhi(phi,clusterPhi[i]);
          m22 += deltaPhi(phi,clusterPhi[i])*deltaPhi(phi,clusterPhi[i]);
          clusterXSpread[i] += (iter->x - clusterX[i]) * (iter->x - clusterX[i]);
          clusterYSpread[i] += (iter->y - clusterY[i]) * (iter->y - clusterY[i]);
          clusterZSpread[i] += (iter->z - clusterZ[i]) * (iter->z - clusterZ[i]);

      }
    }
    float a = (m11+m22)/2;
    float b = 0.5*sqrt((m11+m22)*(m11+m22)-4*(m11*m22-m12*m12));
    clusterXSpread[i] = sqrt(clusterXSpread[i]/(float)clusterSize[i]);
    clusterYSpread[i] = sqrt(clusterYSpread[i]/(float)clusterSize[i]);
    clusterZSpread[i] = sqrt(clusterZSpread[i]/(float)clusterSize[i]);
    clusterEtaSpread[i] = sqrt(m11/clusterSize[i]);
    clusterEtaPhiSpread[i] = sqrt(abs(m12)/clusterSize[i]);
    clusterPhiSpread[i] = sqrt(m22/clusterSize[i]);
    clusterMajorAxis[i] = sqrt((a+b)/clusterSize[i]);
    clusterMinorAxis[i] = sqrt((a-b)/clusterSize[i]);
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
