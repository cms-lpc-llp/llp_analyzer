#include "DBSCAN.h"
#include "TMath.h"
#include <iostream>
#include "TVector3.h"
#include "TGraph.h"
#include "TF1.h"
struct largest_nCsc_cluster_
{
  inline bool operator() (const cscCluster& c1, const cscCluster& c2){return c1.nCscSegments > c2.nCscSegments;}
} largest_nCsc_cluster;

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
    nClusters = clusterID-1;
    return (clusterID-1);
}
void DBSCAN::clear_clusters(){
  // nClusters = 0;
  clusterSize.clear();
  cscLabels.clear();
  clusterEta.clear();
  clusterPhi.clear();
  clusterX.clear();
  clusterY.clear();
  clusterZ.clear();
  clusterTime.clear();
  clusterTimeWire.clear();
  clusterTimeWirePruned.clear();

  clusterTimeTotal.clear();
  clusterMajorAxis.clear();
  clusterMinorAxis.clear();
  clusterXSpread.clear();
  clusterYSpread.clear();
  clusterXYSpread.clear();
  clusterRSpread.clear();
  clusterZSpread.clear();
  clusterTimeSpread.clear();
  clusterTimeTotalSpread.clear();
  clusterTimeTotalSpreadPruned.clear();

  clusterTimeWireSpread.clear();
  clusterEtaPhiSpread.clear();
  clusterEtaSpread.clear();
  clusterPhiSpread.clear();
  clusterDeltaRSpread.clear();

  clusterXSpread_corr.clear();
  clusterYSpread_corr.clear();
  clusterXYSpread_corr.clear();
  clusterEtaSpread_corr.clear();
  clusterPhiSpread_corr.clear();
  clusterEtaPhiSpread_corr.clear();
  clusterRSpread_corr.clear();
  clusterVertexR.clear();
  clusterVertexZ.clear();
  clusterVertexDis.clear();
  clusterVertexChi2.clear();
  clusterVertexN.clear();
  clusterVertexN1cm.clear();
  clusterVertexN5cm.clear();
  clusterVertexN10cm.clear();
  clusterVertexN15cm.clear();
  clusterVertexN20cm.clear();
  CscCluster.clear();
}
int DBSCAN::result(){

  // for (unsigned int i = 0;i < m_pointSize;i++)
  // {
  //   cscLabels.push_back(m_points[i].clusterID);
  //
  // }
  for(int i = 0; i < nClusters; i++)
  {
    float avg_x(0.0), avg_y(0.0), avg_z(0.0), avg_tWire(0.0), avg_tWirePruned(0.0), avg_t(0.0), avg_tTotal(0.0),tTotalSpreadPruned(0.0);
    float avg_eta(0.0), avg_phi(0.0);
    int size(0);
    vector<float> wireTimes;
    vector<float> stripTimes;

    vector<Point>::iterator iter;
    for(iter = m_points.begin(); iter != m_points.end(); ++iter)
    {

      if ( iter->clusterID == i+1 )
      {

          avg_x += iter->x;
          avg_y += iter->y;
          avg_z += iter->z;
          avg_t += iter->t;
          avg_tWire += iter->twire;
          wireTimes.push_back(iter->twire);
          stripTimes.push_back(iter->t);
          size ++;
      }
    }
    avg_x = avg_x/size;
    avg_y = avg_y/size;
    avg_z = avg_z/size;
    avg_tTotal = (avg_t + avg_tWire)/(2 * size);
    avg_t = avg_t/size;
    avg_tWire = avg_tWire/size;

    // prune wire time
    //The wire times have a long tail that has to be pruned.  The strip times (tpeak) are fine
    bool modified = true;
    while (modified) {
      modified = false;
      double maxDiff = -1;
      std::vector<float>::iterator maxHit;
      for (std::vector<float>::iterator itWT = wireTimes.begin(); itWT != wireTimes.end(); ++itWT) {
        float diff = fabs(*itWT - avg_tTotal);
        if (diff > maxDiff) {
          maxDiff = diff;
          maxHit = itWT;
        }
      }
      if (maxDiff > 26) {
        int N = size + wireTimes.size();
        avg_tTotal = (avg_tTotal * N - (*maxHit)) / (N - 1);
        wireTimes.erase(maxHit);
        modified = true;
      }
    }
    for (std::vector<float>::iterator itWT = wireTimes.begin(); itWT != wireTimes.end(); ++itWT) {
      tTotalSpreadPruned += (*itWT - avg_tTotal)*(*itWT - avg_tTotal);
      avg_tWirePruned += *itWT;
    }
    for (std::vector<float>::iterator itWT = stripTimes.begin(); itWT != stripTimes.end(); ++itWT) {
      tTotalSpreadPruned += (*itWT - avg_tTotal)*(*itWT - avg_tTotal);
    }
    tTotalSpreadPruned = sqrt(tTotalSpreadPruned/(stripTimes.size()+wireTimes.size()));
    avg_tWirePruned = avg_tWirePruned/wireTimes.size();
    // calculate cluster eta and phi
    avg_phi = atan(avg_y/avg_x);
    if  (avg_x < 0.0){
      avg_phi = TMath::Pi() + avg_phi;
    }
    avg_phi = deltaPhi(avg_phi,0.0);
    avg_eta = atan(sqrt(pow(avg_x,2)+pow(avg_y,2))/abs(avg_z));
    avg_eta = -1.0*TMath::Sign(1.0, avg_z)*log(tan(avg_eta/2));

    clusterEta.push_back(avg_eta);
    clusterPhi.push_back(avg_phi);
    clusterX.push_back(avg_x);
    clusterY.push_back(avg_y);
    clusterZ.push_back(avg_z);
    clusterTime.push_back(avg_t);
    clusterTimeWire.push_back(avg_tWire);
    clusterTimeWirePruned.push_back(avg_tWirePruned);
    clusterTimeTotal.push_back(avg_tTotal);
    clusterSize.push_back(size);
    clusterTimeTotalSpreadPruned.push_back(tTotalSpreadPruned);

  }
  return 0;
}
int DBSCAN::clusterMoments()
{


  for(int i = 0; i < nClusters; i++)
  {
    float m11(0.0), m12(0.0), m22(0.0);
    float XSpread(0.0), YSpread(0.0), ZSpread(0.0), TSpread(0.0), TTotalSpread(0.0), TWireSpread(0.0), XYSpread(0.0), RSpread(0.0), DeltaRSpread(0.0);

    // float XSpread_phi0p5(0.0), YSpread_phi0p5(0.0), XYSpread_phi0p5(0.0), PhiSpread_phi0p5(0.0), EtaPhiSpread_phi0p5(0.0);
    // float XSpread_phi0p55(0.0), YSpread_phi0p55(0.0), XYSpread_phi0p55(0.0), PhiSpread_phi0p55(0.0), EtaPhiSpread_phi0p55(0.0);
    // float XSpread_phi0p6(0.0), YSpread_phi0p6(0.0), XYSpread_phi0p6(0.0), PhiSpread_phi0p6(0.0), EtaPhiSpread_phi0p6(0.0);
    // float XSpread_phi0p65(0.0), YSpread_phi0p65(0.0), XYSpread_phi0p65(0.0), PhiSpread_phi0p65(0.0), EtaPhiSpread_phi0p65(0.0);
    // float XSpread_phi0p7(0.0), YSpread_phi0p7(0.0), XYSpread_phi0p7(0.0), PhiSpread_phi0p7(0.0), EtaPhiSpread_phi0p7(0.0);
    // float XSpread_phi0p75(0.0), YSpread_phi0p75(0.0), XYSpread_phi0p75(0.0), PhiSpread_phi0p75(0.0), EtaPhiSpread_phi0p75(0.0);
    //
    // float XSpread_phi0p7_r1p3(0.0), YSpread_phi0p7_r1p3(0.0), XYSpread_phi0p7_r1p3(0.0), RSpread_phi0p7_r1p3(0.0),
    // PhiSpread_phi0p7_r1p3(0.0),EtaSpread_phi0p7_r1p3(0.0), EtaPhiSpread_phi0p7_r1p3(0.0);
    // float XSpread_phi0p7_r1p2(0.0), YSpread_phi0p7_r1p2(0.0), XYSpread_phi0p7_r1p2(0.0), RSpread_phi0p7_r1p2(0.0),
    // PhiSpread_phi0p7_r1p2(0.0),EtaSpread_phi0p7_r1p2(0.0), EtaPhiSpread_phi0p7_r1p2(0.0);
    // float XSpread_phi0p7_r1p1(0.0), YSpread_phi0p7_r1p1(0.0), XYSpread_phi0p7_r1p1(0.0), RSpread_phi0p7_r1p1(0.0),
    // PhiSpread_phi0p7_r1p1(0.0),EtaSpread_phi0p7_r1p1(0.0), EtaPhiSpread_phi0p7_r1p1(0.0);
    // float XSpread_phi0p7_r1p15(0.0), YSpread_phi0p7_r1p15(0.0), XYSpread_phi0p7_r1p15(0.0), RSpread_phi0p7_r1p15(0.0),
    // PhiSpread_phi0p7_r1p15(0.0),EtaSpread_phi0p7_r1p15(0.0), EtaPhiSpread_phi0p7_r1p15(0.0);
    // float XSpread_phi0p7_r1p25(0.0), YSpread_phi0p7_r1p25(0.0), XYSpread_phi0p7_r1p25(0.0), RSpread_phi0p7_r1p25(0.0),
    // PhiSpread_phi0p7_r1p25(0.0),EtaSpread_phi0p7_r1p25(0.0), EtaPhiSpread_phi0p7_r1p25(0.0);
    //
    // float XSpread_r1p2(0.0), YSpread_r1p2(0.0), XYSpread_r1p2(0.0), RSpread_r1p2(0.0),
    // PhiSpread_r1p2(0.0), EtaSpread_r1p2(0.0), EtaPhiSpread_r1p2(0.0);


    vector<vector<float>> XSpread_corr( N_phicorr , vector<float> (N_rcorr, 0));
    vector<vector<float>> YSpread_corr( N_phicorr , vector<float> (N_rcorr, 0));
    vector<vector<float>> XYSpread_corr( N_phicorr , vector<float> (N_rcorr, 0));
    vector<vector<float>> EtaSpread_corr( N_phicorr , vector<float> (N_rcorr, 0));
    vector<vector<float>> PhiSpread_corr( N_phicorr , vector<float> (N_rcorr, 0));
    vector<vector<float>> EtaPhiSpread_corr( N_phicorr , vector<float> (N_rcorr, 0));
    vector<vector<float>> RSpread_corr( N_phicorr , vector<float> (N_rcorr, 0));


    vector<Point>::iterator iter;
    for(iter = m_points.begin(); iter != m_points.end(); ++iter)
    {


      if ( iter->clusterID == i+1 )
      {
          // float phi = atan(iter->y/iter->x);
          // if  (iter->x < 0.0){
          //   phi = TMath::Pi() + phi;
          // }
          // phi = deltaPhi(phi,0.0);
          // float eta = atan(sqrt(pow(iter->x,2)+pow(iter->y,2))/abs(iter->z));
          // eta = -1.0*TMath::Sign(1.0, iter->z)*log(tan(eta/2));
          m11 += (iter->eta-clusterEta[i])*(iter->eta-clusterEta[i]);
          m12 += (iter->eta-clusterEta[i])* deltaPhi(iter->phi,clusterPhi[i]);
          m22 += deltaPhi(iter->phi,clusterPhi[i])*deltaPhi(iter->phi,clusterPhi[i]);
          DeltaRSpread +=  pow(deltaR(clusterEta[i], clusterPhi[i], iter->eta, iter->phi),2);

          XYSpread += (iter->x - clusterX[i])*(iter->y - clusterY[i]);
          XSpread += (iter->x - clusterX[i]) * (iter->x - clusterX[i]);
          YSpread += (iter->y - clusterY[i]) * (iter->y - clusterY[i]);
          ZSpread += (iter->z - clusterZ[i]) * (iter->z - clusterZ[i]);
          TSpread += (iter->t - clusterTime[i]) * (iter->t - clusterTime[i]);
          TWireSpread += (iter->twire - clusterTimeWire[i]) * (iter->twire - clusterTimeWire[i]);
          TTotalSpread += (iter->t - clusterTimeTotal[i]) * (iter->t - clusterTimeTotal[i])+
          (iter->twire - clusterTimeTotal[i]) * (iter->twire - clusterTimeTotal[i]);


          float radius = sqrt(pow(iter->x, 2) + pow(iter->y, 2));
          RSpread += pow(radius-sqrt(clusterX[i]*clusterX[i]+clusterY[i]*clusterY[i]),2);


          // phi+R correction
          for (int j = 0; j < N_phicorr; j++)
          {
            for (int k = 0; k < N_rcorr; k++)
            {
              float phi = deltaPhi(iter->phi,clusterPhi[i])*phi_corr[j] + clusterPhi[i];
              float r = (radius-sqrt(clusterX[i]*clusterX[i]+clusterY[i]*clusterY[i]))*r_corr[k] + sqrt(clusterX[i]*clusterX[i]+clusterY[i]*clusterY[i]);
              if (r < 0) cout<<"r is negative!!!"<<r<<","<<r_corr[k]<<","<<j<<","<<k<<endl;
              float x = r*cos(phi);
              float y = r*sin(phi);
              double theta = atan(r/abs(iter->z));
              double eta = -1.0*TMath::Sign(1.0, iter->z)*log(tan(theta/2));
              PhiSpread_corr[j][k] += deltaPhi(phi,clusterPhi[i])*deltaPhi(phi,clusterPhi[i]);
              EtaPhiSpread_corr[j][k] += (eta-clusterEta[i])* deltaPhi(phi,clusterPhi[i]);
              EtaSpread_corr[j][k] += (eta-clusterEta[i])* (eta-clusterEta[i]);
              XYSpread_corr[j][k] += (x - clusterX[i])*(y - clusterY[i]);
              XSpread_corr[j][k] += (x - clusterX[i]) * (x - clusterX[i]);
              YSpread_corr[j][k] += (y - clusterY[i]) * (y - clusterY[i]);
              RSpread_corr[j][k] += pow(r - sqrt(clusterX[i]*clusterX[i]+clusterY[i]*clusterY[i]), 2);
            }
          }




      }
    }
    float a = (m11+m22)/2;
    float b = 0.5*sqrt((m11+m22)*(m11+m22)-4*(m11*m22-m12*m12));
    clusterXSpread.push_back(sqrt(XSpread/(float)clusterSize[i]));
    clusterYSpread.push_back(sqrt(YSpread/(float)clusterSize[i]));
    clusterZSpread.push_back(sqrt(ZSpread/(float)clusterSize[i]));
    clusterRSpread.push_back(sqrt(RSpread/(float)clusterSize[i]));
    clusterDeltaRSpread.push_back(sqrt(DeltaRSpread/(float)clusterSize[i]));
    clusterXYSpread.push_back(sqrt(abs(XYSpread)/(float)clusterSize[i]));
    clusterTimeSpread.push_back(sqrt(TSpread/(float)clusterSize[i]));
    clusterTimeTotalSpread.push_back(sqrt(TTotalSpread/(float)clusterSize[i]/2));
    clusterTimeWireSpread.push_back(sqrt(TWireSpread/(float)clusterSize[i]));
    clusterEtaSpread.push_back(sqrt(abs(m11)/clusterSize[i]));
    clusterEtaPhiSpread.push_back(sqrt(abs(m12)/clusterSize[i]));
    clusterPhiSpread.push_back(sqrt(abs(m22)/clusterSize[i]));
    clusterMajorAxis.push_back(sqrt((a+b)/clusterSize[i]));
    clusterMinorAxis.push_back(sqrt((a-b)/clusterSize[i]));
    for (int j = 0; j < N_phicorr; j++)
    {

      for (int k = 0; k < N_rcorr; k++)
      {
        PhiSpread_corr[j][k] = sqrt(abs(PhiSpread_corr[j][k])/(float)clusterSize[i]);
        EtaPhiSpread_corr[j][k] = sqrt(abs(EtaPhiSpread_corr[j][k])/(float)clusterSize[i]);
        EtaSpread_corr[j][k] = sqrt(abs(EtaSpread_corr[j][k])/(float)clusterSize[i]);
        XYSpread_corr[j][k] = sqrt(abs(XYSpread_corr[j][k])/(float)clusterSize[i]);
        XSpread_corr[j][k] = sqrt(abs(XSpread_corr[j][k])/(float)clusterSize[i]);
        YSpread_corr[j][k] = sqrt(abs(YSpread_corr[j][k])/(float)clusterSize[i]);
        RSpread_corr[j][k] = sqrt(abs(RSpread_corr[j][k])/(float)clusterSize[i]);
      }
    }
    clusterPhiSpread_corr.push_back(PhiSpread_corr);
    clusterEtaPhiSpread_corr.push_back(EtaPhiSpread_corr);
    clusterEtaSpread_corr.push_back(EtaSpread_corr);
    clusterXYSpread_corr.push_back(XYSpread_corr);
    clusterXSpread_corr.push_back(XSpread_corr);
    clusterYSpread_corr.push_back(YSpread_corr);
    clusterRSpread_corr.push_back(RSpread_corr);


  }
  return 0;
}
//input x,y,z of segments in cluster, avgX, avgY, avgZ
void DBSCAN::sort_clusters()
{
  for(int i = 0; i < nClusters; i++){
    vector<int>cscStations;
    vector<int>cscStations_copy;
    vector<int>cscChambers;
    vector<int>cscChambers_copy;
    vector<int>cscLayersPlus11;
    vector<int>cscLayersPlus12;
    vector<int>cscLayersPlus13;
    vector<int>cscLayersPlus21;
    vector<int>cscLayersPlus22;
    vector<int>cscLayersPlus31;
    vector<int>cscLayersPlus32;
    vector<int>cscLayersPlus41;
    vector<int>cscLayersPlus42;

    vector<int>cscLayersMinus11;
    vector<int>cscLayersMinus12;
    vector<int>cscLayersMinus13;
    vector<int>cscLayersMinus21;
    vector<int>cscLayersMinus22;
    vector<int>cscLayersMinus31;
    vector<int>cscLayersMinus32;
    vector<int>cscLayersMinus41;
    vector<int>cscLayersMinus42;

    vector<int>segment_index;
    int nSegments_Me11 = 0;
    int nSegments_Me12 = 0;

    cscCluster tmpCluster;
    tmpCluster.nCscSegmentChamberPlus11 = 0;
    tmpCluster.nCscSegmentChamberPlus12 = 0;
    tmpCluster.nCscSegmentChamberPlus13 = 0;
    tmpCluster.nCscSegmentChamberPlus21 = 0;
    tmpCluster.nCscSegmentChamberPlus22 = 0;
    tmpCluster.nCscSegmentChamberPlus31 = 0;
    tmpCluster.nCscSegmentChamberPlus32 = 0;
    tmpCluster.nCscSegmentChamberPlus41 = 0;
    tmpCluster.nCscSegmentChamberPlus42 = 0;
    tmpCluster.nCscSegmentChamberMinus11 = 0;
    tmpCluster.nCscSegmentChamberMinus12 = 0;
    tmpCluster.nCscSegmentChamberMinus13 = 0;
    tmpCluster.nCscSegmentChamberMinus21 = 0;
    tmpCluster.nCscSegmentChamberMinus22 = 0;
    tmpCluster.nCscSegmentChamberMinus31 = 0;
    tmpCluster.nCscSegmentChamberMinus32 = 0;
    tmpCluster.nCscSegmentChamberMinus41 = 0;
    tmpCluster.nCscSegmentChamberMinus42 = 0;

    // tmpCluster.nDtSegmentStation1 = 0;
    // tmpCluster.nDtSegmentStation2 = 0;
    // tmpCluster.nDtSegmentStation3 = 0;
    // tmpCluster.nDtSegmentStation4 = 0;
    for(unsigned int l=0; l < m_pointSize; l++){
      if (m_points[l].clusterID == i+1){
        cscStations.push_back(m_points[l].station);
        cscChambers.push_back(m_points[l].chamber);
        cscStations_copy.push_back(m_points[l].station);
        cscChambers_copy.push_back(m_points[l].chamber);
        segment_index.push_back(l);

        if (m_points[l].chamber == 11) cscLayersPlus11.push_back(m_points[l].layer);
        if (m_points[l].chamber == 12) cscLayersPlus12.push_back(m_points[l].layer);
        if (m_points[l].chamber == 13) cscLayersPlus13.push_back(m_points[l].layer);
        if (m_points[l].chamber == 21) cscLayersPlus21.push_back(m_points[l].layer);
      	if (m_points[l].chamber == 22) cscLayersPlus22.push_back(m_points[l].layer);
        if (m_points[l].chamber == 31) cscLayersPlus31.push_back(m_points[l].layer);
      	if (m_points[l].chamber == 32) cscLayersPlus32.push_back(m_points[l].layer);
        if (m_points[l].chamber == 41) cscLayersPlus41.push_back(m_points[l].layer);
      	if (m_points[l].chamber == 42) cscLayersPlus42.push_back(m_points[l].layer);

        if (m_points[l].chamber == -11) cscLayersMinus11.push_back(m_points[l].layer);
        if (m_points[l].chamber == -12) cscLayersMinus12.push_back(m_points[l].layer);
        if (m_points[l].chamber == -13) cscLayersMinus13.push_back(m_points[l].layer);
        if (m_points[l].chamber == -21) cscLayersMinus21.push_back(m_points[l].layer);
      	if (m_points[l].chamber == -22) cscLayersMinus22.push_back(m_points[l].layer);
        if (m_points[l].chamber == -31) cscLayersMinus31.push_back(m_points[l].layer);
      	if (m_points[l].chamber == -32) cscLayersMinus32.push_back(m_points[l].layer);
        if (m_points[l].chamber == -41) cscLayersMinus41.push_back(m_points[l].layer);
      	if (m_points[l].chamber == -42) cscLayersMinus42.push_back(m_points[l].layer);

        if (abs(m_points[l].chamber) == 11) nSegments_Me11++;
        if (abs(m_points[l].chamber) == 12) nSegments_Me12++;
      	if (m_points[l].chamber == 11) tmpCluster.nCscSegmentChamberPlus11++;
      	if (m_points[l].chamber == 12) tmpCluster.nCscSegmentChamberPlus12++;
      	if (m_points[l].chamber == 13) tmpCluster.nCscSegmentChamberPlus13++;
      	if (m_points[l].chamber == 21) tmpCluster.nCscSegmentChamberPlus21++;
      	if (m_points[l].chamber == 22) tmpCluster.nCscSegmentChamberPlus22++;
      	if (m_points[l].chamber == 31) tmpCluster.nCscSegmentChamberPlus31++;
      	if (m_points[l].chamber == 32) tmpCluster.nCscSegmentChamberPlus32++;
      	if (m_points[l].chamber == 41) tmpCluster.nCscSegmentChamberPlus41++;
      	if (m_points[l].chamber == 42) tmpCluster.nCscSegmentChamberPlus42++;
      	if (m_points[l].chamber == -11) tmpCluster.nCscSegmentChamberMinus11++;
      	if (m_points[l].chamber == -12) tmpCluster.nCscSegmentChamberMinus12++;
      	if (m_points[l].chamber == -13) tmpCluster.nCscSegmentChamberMinus13++;
      	if (m_points[l].chamber == -21) tmpCluster.nCscSegmentChamberMinus21++;
      	if (m_points[l].chamber == -22) tmpCluster.nCscSegmentChamberMinus22++;
      	if (m_points[l].chamber == -31) tmpCluster.nCscSegmentChamberMinus31++;
      	if (m_points[l].chamber == -32) tmpCluster.nCscSegmentChamberMinus32++;
      	if (m_points[l].chamber == -41) tmpCluster.nCscSegmentChamberMinus41++;
      	if (m_points[l].chamber == -42) tmpCluster.nCscSegmentChamberMinus42++;


        // if (m_points[l].station == 1) tmpCluster.nDtSegmentStation1++;
        // if (m_points[l].station == 2) tmpCluster.nDtSegmentStation2++;
        // if (m_points[l].station == 3) tmpCluster.nDtSegmentStation3++;
        // if (m_points[l].station == 4) tmpCluster.nDtSegmentStation4++;

      }
    }
    tmpCluster.x = clusterX[i];
    tmpCluster.y = clusterY[i];
    tmpCluster.z = clusterZ[i];
    tmpCluster.eta = clusterEta[i];
    tmpCluster.phi = clusterPhi[i];
    tmpCluster.t = clusterTime[i];
    tmpCluster.tTotal = clusterTimeTotal[i];
    tmpCluster.tWire = clusterTimeWire[i];
    tmpCluster.tWirePruned = clusterTimeWirePruned[i];

    tmpCluster.MajorAxis = clusterMajorAxis[i];
    tmpCluster.MinorAxis = clusterMinorAxis[i];
    tmpCluster.XSpread = clusterXSpread[i];
    tmpCluster.XYSpread = clusterXYSpread[i];
    tmpCluster.RSpread = clusterRSpread[i];


    tmpCluster.YSpread = clusterYSpread[i];
    tmpCluster.ZSpread = clusterZSpread[i];
    tmpCluster.TSpread = clusterTimeSpread[i];
    tmpCluster.TWireSpread = clusterTimeWireSpread[i];
    tmpCluster.TTotalSpread = clusterTimeTotalSpread[i];
    tmpCluster.TTotalSpreadPruned = clusterTimeTotalSpreadPruned[i];

    tmpCluster.EtaPhiSpread = clusterEtaPhiSpread[i];
    tmpCluster.EtaSpread = clusterEtaSpread[i];
    tmpCluster.DeltaRSpread = clusterDeltaRSpread[i];
    tmpCluster.PhiSpread = clusterPhiSpread[i];
    tmpCluster.nCscSegments = clusterSize[i];
    tmpCluster.Me11Ratio = 1.0*nSegments_Me11/clusterSize[i];
    tmpCluster.Me12Ratio = 1.0*nSegments_Me12/clusterSize[i];

    // tmpCluster.cscChambers = cscChambers;
    for (int j = 0; j < N_phicorr; j++)
    {
      for (int k = 0; k < N_rcorr; k++)
      {
        tmpCluster.XSpread_corr[j][k] = clusterXSpread_corr[i][j][k];
        tmpCluster.YSpread_corr[j][k] = clusterYSpread_corr[i][j][k];
        tmpCluster.XYSpread_corr[j][k] = clusterXYSpread_corr[i][j][k];
        tmpCluster.EtaSpread_corr[j][k] = clusterEtaSpread_corr[i][j][k];
        tmpCluster.PhiSpread_corr[j][k] = clusterPhiSpread_corr[i][j][k];
        tmpCluster.EtaPhiSpread_corr[j][k] = clusterEtaPhiSpread_corr[i][j][k];
        tmpCluster.RSpread_corr[j][k] = clusterRSpread_corr[i][j][k];
      }
    }



    // count number of cscLayers
    std::sort(cscLayersPlus11.begin(), cscLayersPlus11.end());
    std::sort(cscLayersPlus12.begin(), cscLayersPlus12.end());
    std::sort(cscLayersPlus13.begin(), cscLayersPlus13.end());
    std::sort(cscLayersPlus21.begin(), cscLayersPlus21.end());
    std::sort(cscLayersPlus22.begin(), cscLayersPlus22.end());
    std::sort(cscLayersPlus31.begin(), cscLayersPlus31.end());
    std::sort(cscLayersPlus32.begin(), cscLayersPlus32.end());
    std::sort(cscLayersPlus41.begin(), cscLayersPlus41.end());
    std::sort(cscLayersPlus42.begin(), cscLayersPlus42.end());

    std::sort(cscLayersMinus11.begin(), cscLayersMinus11.end());
    std::sort(cscLayersMinus12.begin(), cscLayersMinus12.end());
    std::sort(cscLayersMinus13.begin(), cscLayersMinus13.end());
    std::sort(cscLayersMinus21.begin(), cscLayersMinus21.end());
    std::sort(cscLayersMinus22.begin(), cscLayersMinus22.end());
    std::sort(cscLayersMinus31.begin(), cscLayersMinus31.end());
    std::sort(cscLayersMinus32.begin(), cscLayersMinus32.end());
    std::sort(cscLayersMinus41.begin(), cscLayersMinus41.end());
    std::sort(cscLayersMinus42.begin(), cscLayersMinus42.end());


    tmpCluster.nLayersChamberPlus11 = std::unique(cscLayersPlus11.begin(), cscLayersPlus11.end())-cscLayersPlus11.begin();
    tmpCluster.nLayersChamberPlus12 = std::unique(cscLayersPlus12.begin(), cscLayersPlus12.end())-cscLayersPlus12.begin();
    tmpCluster.nLayersChamberPlus13 = std::unique(cscLayersPlus13.begin(), cscLayersPlus13.end())-cscLayersPlus13.begin();
    tmpCluster.nLayersChamberPlus21 = std::unique(cscLayersPlus21.begin(), cscLayersPlus21.end())-cscLayersPlus21.begin();
    tmpCluster.nLayersChamberPlus22 = std::unique(cscLayersPlus22.begin(), cscLayersPlus22.end())-cscLayersPlus22.begin();
    tmpCluster.nLayersChamberPlus31 = std::unique(cscLayersPlus31.begin(), cscLayersPlus31.end())-cscLayersPlus31.begin();
    tmpCluster.nLayersChamberPlus32 = std::unique(cscLayersPlus32.begin(), cscLayersPlus32.end())-cscLayersPlus32.begin();
    tmpCluster.nLayersChamberPlus41 = std::unique(cscLayersPlus41.begin(), cscLayersPlus41.end())-cscLayersPlus41.begin();
    tmpCluster.nLayersChamberPlus42 = std::unique(cscLayersPlus42.begin(), cscLayersPlus42.end())-cscLayersPlus42.begin();
    tmpCluster.nLayersChamberMinus11 = std::unique(cscLayersMinus11.begin(), cscLayersMinus11.end())-cscLayersMinus11.begin();
    tmpCluster.nLayersChamberMinus12 = std::unique(cscLayersMinus12.begin(), cscLayersMinus12.end())-cscLayersMinus12.begin();
    tmpCluster.nLayersChamberMinus13 = std::unique(cscLayersMinus13.begin(), cscLayersMinus13.end())-cscLayersMinus13.begin();
    tmpCluster.nLayersChamberMinus21 = std::unique(cscLayersMinus21.begin(), cscLayersMinus21.end())-cscLayersMinus21.begin();
    tmpCluster.nLayersChamberMinus22 = std::unique(cscLayersMinus22.begin(), cscLayersMinus22.end())-cscLayersMinus22.begin();
    tmpCluster.nLayersChamberMinus31 = std::unique(cscLayersMinus31.begin(), cscLayersMinus31.end())-cscLayersMinus31.begin();
    tmpCluster.nLayersChamberMinus32 = std::unique(cscLayersMinus32.begin(), cscLayersMinus32.end())-cscLayersMinus32.begin();
    tmpCluster.nLayersChamberMinus41 = std::unique(cscLayersMinus41.begin(), cscLayersMinus41.end())-cscLayersMinus41.begin();
    tmpCluster.nLayersChamberMinus42 = std::unique(cscLayersMinus42.begin(), cscLayersMinus42.end())-cscLayersMinus42.begin();
    tmpCluster.segment_id = segment_index;

    // count the number of chambers and max chamber segments
    std::vector<int>::iterator chamber_it;
    std::sort(cscChambers.begin(), cscChambers.end());
    chamber_it = std::unique(cscChambers.begin(), cscChambers.end());
    cscChambers.resize( std::distance(cscChambers.begin(),chamber_it) );
    int max_chamber = 999; // station with the maximum number of cscsegment in this cluster
    int max_chamber_segment = 0; // station with the maximum number of cscsegment in this cluster
    tmpCluster.nChamber = 0;
    for (unsigned int l = 0; l < cscChambers.size(); l++)
    {
      int counter = 0;
      for(unsigned int j = 0; j < cscChambers_copy.size(); j ++)
      {
        if (cscChambers_copy[j] == cscChambers[l]) counter++;
      }
      if (counter>max_chamber_segment)
      {
        max_chamber_segment = counter;
        max_chamber = cscChambers[l];
      }
      if(counter>5)tmpCluster.nChamber++;
    }
    tmpCluster.maxChamber = max_chamber;
    tmpCluster.maxChamberSegment = max_chamber_segment;
    // count the number of chambers and max chamber segments
    std::vector<int>::iterator station_it;
    std::sort(cscStations.begin(), cscStations.end());
    station_it = std::unique(cscStations.begin(), cscStations.end());
    cscStations.resize( std::distance(cscStations.begin(),station_it) );//list of unique stations
    int max_station = 999; // station with the maximum number of cscsegment in this cluster
    int max_station_segment = 0; // station with the maximum number of cscsegment in this cluster
    tmpCluster.nStation = 0;
    tmpCluster.nStation5 = 0;
    tmpCluster.nStation10 = 0;
    tmpCluster.nStation10perc = 0;
    tmpCluster.avgStation = 0.0;
    tmpCluster.avgStation5 = 0.0;
    tmpCluster.avgStation10 = 0.0;
    tmpCluster.avgStation10perc = 0.0;
    int nSeg10perc = 0;
    int nSeg5 = 0;
    int nSeg10 = 0;
    for (unsigned int l = 0; l < cscStations.size(); l++)
    {
      int counter = 0;
      for(unsigned int j = 0; j < cscStations_copy.size(); j ++)
      {
        if (cscStations_copy[j] == cscStations[l]) counter++;
      }
      if (counter>max_station_segment)
      {
        max_station_segment = counter;
        max_station = cscStations[l];
      }
      tmpCluster.nStation++;
      tmpCluster.avgStation += counter * cscStations[l];

      if(counter>=5.0){
        tmpCluster.avgStation5 += counter * cscStations[l];
        tmpCluster.nStation5++;
        nSeg5 += counter;

      }
      if(counter>=10.0){
        tmpCluster.avgStation10 += counter * cscStations[l];
        tmpCluster.nStation10++;
        nSeg10 += counter;
      }
      if(1.0*counter/clusterSize[i] > 0.1)
      {
        tmpCluster.avgStation10perc += counter * cscStations[l];
        tmpCluster.nStation10perc++;
        nSeg10perc += counter;
      }
    }
    tmpCluster.avgStation10perc = 1.0* tmpCluster.avgStation10perc/nSeg10perc;
    tmpCluster.avgStation5 = 1.0* tmpCluster.avgStation5/nSeg5;
    tmpCluster.avgStation10 = 1.0* tmpCluster.avgStation10/nSeg10;
    tmpCluster.avgStation = 1.0* tmpCluster.avgStation/clusterSize[i];

    tmpCluster.maxStation = max_station;
    tmpCluster.maxStationSegment = max_station_segment;

    CscCluster.push_back(tmpCluster);


  }

  //sort the clusters by size
  sort(CscCluster.begin(), CscCluster.end(), largest_nCsc_cluster);

}


void DBSCAN::merge_clusters()
{
  // clear all the cluster variables
  //change cluster ID of points
  bool modified = true;
  while(modified){
    modified = false;
    float mindR = 15;
    int cluster1 = 999;
    int cluster2 = 999;

    for(unsigned int i = 0; i < clusterEta.size(); i++){
      for(unsigned int j = i+1; j < clusterEta.size(); j++){
        float current_dR = deltaR(clusterEta[i], clusterPhi[i], clusterEta[j], clusterPhi[j]);
        if(current_dR<mindR)
        {
          mindR = current_dR;
          cluster1 = i;
          cluster2 = j;

        }
      }
    }
    if (mindR < MERGE_CLUSTER_DR){
      vector<Point>::iterator iter;
      int count = 0;
      for(iter = m_points.begin(); iter != m_points.end(); ++iter)
      {
        if ( iter->clusterID == cluster2+1 ){
          iter->clusterID = cluster1+1;
          count++;
        }
        if ( iter->clusterID > cluster2+1 )iter->clusterID = iter->clusterID-1;
      }
      clusterEta.erase(clusterEta.begin() + cluster2);
      clusterPhi.erase(clusterPhi.begin() + cluster2);
      nClusters--;
      modified = true;
      // can't use cscCluster, because its sorted, but the other vectors and clusterID are not sorted.
    }
  }
  clear_clusters(); // clear all the vectors, but nClusters is kept to keep track of the number of clusters.
}

int DBSCAN::vertexing()
{
  for(int i = 0; i < nClusters; i++){
    TVector3 vecDir;
    TVector3 vecCsc;
    TF1 *fit = new TF1();

    TGraph *gr = new TGraph(clusterSize[i]);

    for(unsigned int j = 0;j < m_pointSize;j++)
    {
      if ( m_points[j].clusterID == i+1 ){
        vecCsc.SetXYZ(m_points[j].x, m_points[j].y, m_points[j].z);
        vecDir.SetXYZ(m_points[j].dirX, m_points[j].dirY, m_points[j].dirZ);
        double slope = (vecCsc.Perp() - (vecCsc+vecDir).Perp()) / (vecCsc.Z() - (vecCsc+vecDir).Z());
        double beta = -1.0*vecCsc.Z() * slope + vecCsc.Perp();
        gr->SetPoint(j,slope,beta);
      }

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
    int vertexN1cm(0), vertexN5cm(0), vertexN10cm(0), vertexN15cm(0), vertexN20cm(0), vertexN(0);
    float vertexChi2(0.0), vertexDis(0.0), vertexZ(0.0), vertexR(0.0);


    // cout << "starting with " << gr->GetN() << " segments" << endl;
    while (!goodVtx && gr->GetN()>=3){
    	maxDistance = 0.0;
      // gr->Fit("pol 1");
      // fit = gr->GetFunction("pol 1");
    	gr->Fit("1++x","Q");
    	fit = gr->GetFunction("1++x");
    	for (int j=0; j<gr->GetN(); j++){
    	  // Distance from point (x_0, y_0) to line ax + by + c = 0 is |ax_0 + by_0 + c| / sqrt(a^2 + b^2)
    	  // Here: a = slope = "X"; b = -1,;c = intercept = "Y"; (x_0, y_0) = (z, r) = (-p[1], p[0]);
    	  distance = abs(gr->GetX()[j]*fit->GetParameter(1)*-1.0 + -1.0*fit->GetParameter(0) + gr->GetY()[j]);
    	  distance = distance / sqrt(pow(gr->GetX()[j],2)+pow(-1.0,2));
    	  if (distance > maxDistance){
    	    maxDistance = distance;
    	    farSegment = j;
    	  }
    	}
    	if (maxDistance < 30.0){
    	  goodVtx = true;
    	}
    	else {
    	  gr->RemovePoint(farSegment);
    	}
    }
    // count number of segment in distance 1,5,10,20cm
    for (int j=0; j<gr->GetN(); j++){
      // Distance from point (x_0, y_0) to line ax + by + c = 0 is |ax_0 + by_0 + c| / sqrt(a^2 + b^2)
      // Here: a = slope = "X"; b = -1,;c = intercept = "Y"; (x_0, y_0) = (z, r) = (-p[1], p[0]);
      distance = abs(gr->GetX()[j]*fit->GetParameter(1)*-1.0 + -1.0*fit->GetParameter(0) + gr->GetY()[j]);
      distance = distance / sqrt(pow(gr->GetX()[j],2)+pow(-1.0,2));
      if (distance < 1.0) vertexN1cm ++;
      if (distance < 5.0) vertexN5cm ++;
      if (distance < 10.0) vertexN10cm ++;
      if (distance < 15.0) vertexN15cm ++;
      if (distance < 20.0) vertexN20cm ++;

    }
    if (goodVtx && gr->GetN()>=3){
      vertexR = fit->GetParameter(0);
      vertexZ = -1.0*fit->GetParameter(1);
      vertexN = gr->GetN();
      vertexDis = maxDistance;
      vertexChi2 = fit->GetChisquare() ;
    }
    else{
      vertexR = 0.0;
      vertexZ = 0.0;
      vertexN = 0;
      vertexDis =-999;
      vertexChi2 = -999.;

    }
    clusterVertexR.push_back(vertexR);
    clusterVertexZ.push_back(vertexZ);
    clusterVertexN.push_back(vertexN);
    clusterVertexDis.push_back(vertexDis);
    clusterVertexChi2.push_back(vertexChi2);

    clusterVertexN1cm.push_back(vertexN1cm);
    clusterVertexN5cm.push_back(vertexN5cm);
    clusterVertexN10cm.push_back(vertexN10cm);
    clusterVertexN15cm.push_back(vertexN15cm);
    clusterVertexN20cm.push_back(vertexN20cm);
  }
  // cout << "vertex? " << goodVtx << ", # segments = " << gr->GetN() << ", maxDistance = " << maxDistance << endl;
  // cout << "vertex R, Z, N:  " << clusterVertexR <<", " << clusterVertexZ << ",  " << clusterVertexN << endl;
  // cout << "maxDistance = " << maxDistance << endl;

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

            // if (m_points.at(*iterSeeds).x == point.x && m_points.at(*iterSeeds).y == point.y && m_points.at(*iterSeeds).z == point.z )
            if (m_points.at(*iterSeeds).eta == point.eta && m_points.at(*iterSeeds).phi == point.phi )
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
    // return sqrt(pow(pointCore.x - pointTarget.x,2)+pow(pointCore.y - pointTarget.y,2)+pow(pointCore.z - pointTarget.z,2));
    // return sqrt(pow(pointCore.eta - pointTarget.eta,2)+pow(deltaPhi(pointCore.phi, pointTarget.phi),2)+pow((pointCore.t - pointTarget.t)/100.0,2));
    return sqrt(pow(pointCore.eta - pointTarget.eta,2)+pow(deltaPhi(pointCore.phi, pointTarget.phi),2));

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
double DBSCAN::deltaR(double eta1, double phi1, double eta2, double phi2) {
  double dphi = deltaPhi(phi1,phi2);
  double deta = eta1 - eta2;
  return sqrt( dphi*dphi + deta*deta);
}
