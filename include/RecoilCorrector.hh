//================================================================================================
//
// Class and tools for recoil corrections
//
//  * Defines RecoilCorrector class to access and apply MET corrections based on Z recoil
//
//________________________________________________________________________________________________

#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TRandom.h>
#include <TString.h>
#include <TVector2.h>
#include <TFitResult.h>
#include <iostream>


//--------------------------------------------------------------------------------------------------
class RecoilCorrector
{
public:
  RecoilCorrector(TString fname,                                                 // file name of recoil fits
                  TString fname_Wp="", TString fname_Wm="", TString fname_Z="",  // recoil fits from MC for W/Z correction
		  Int_t iSeed=0xDEADBEEF);                                       // seed for random number generator
  ~RecoilCorrector();
  void Correct(Double_t &pfmet, Double_t &pfmetphi,  // reference to variables to store corrected MET and phi(MET)
               Double_t genWPt, Double_t genWPhi,    // GEN W boson pT and phi
	       Double_t lepPt,  Double_t lepPhi,     // lepton pT and phi
	       Double_t nsigma,                      // # of sigmas on fit uncertainty for systematics studies (0 = nominal correction)
	       Int_t charge);                        // lepton charge

protected:
  TF1 *fcnPFu1mean, *fcnPFu1sigma1, *fcnPFu1sigma2, *fcnPFu1sigma0;
  TF1 *fcnPFu2mean, *fcnPFu2sigma1, *fcnPFu2sigma2, *fcnPFu2sigma0;
  
  TFitResult *fitresPFu1mean, *fitresPFu1sigma1, *fitresPFu1sigma2,  *fitresPFu1sigma0;
  TFitResult *fitresPFu2mean, *fitresPFu2sigma1, *fitresPFu2sigma2,  *fitresPFu2sigma0;
  
  Double_t pfu1meanCov[2][2];
  Double_t pfu2meanCov[2][2];
  Double_t pfu1sigma1Cov[3][3];
  Double_t pfu2sigma1Cov[3][3];
  Double_t pfu1sigma2Cov[3][3];
  Double_t pfu2sigma2Cov[3][3];
  Double_t pfu1sigma0Cov[3][3];
  Double_t pfu2sigma0Cov[3][3]; 

  TF1 *fcnPFu1mean_Wp, *fcnPFu1sigma1_Wp, *fcnPFu1sigma2_Wp, *fcnPFu1sigma0_Wp;
  TF1 *fcnPFu2mean_Wp, *fcnPFu2sigma1_Wp, *fcnPFu2sigma2_Wp, *fcnPFu2sigma0_Wp;

  TF1 *fcnPFu1mean_Wm, *fcnPFu1sigma1_Wm, *fcnPFu1sigma2_Wm, *fcnPFu1sigma0_Wm;
  TF1 *fcnPFu2mean_Wm, *fcnPFu2sigma1_Wm, *fcnPFu2sigma2_Wm, *fcnPFu2sigma0_Wm;

  TF1 *fcnPFu1mean_Z, *fcnPFu1sigma1_Z, *fcnPFu1sigma2_Z, *fcnPFu1sigma0_Z;
  TF1 *fcnPFu2mean_Z, *fcnPFu2sigma1_Z, *fcnPFu2sigma2_Z, *fcnPFu2sigma0_Z;
};

//--------------------------------------------------------------------------------------------------
Double_t sigmaFunc(Double_t *x, Double_t *par) {
  // par[0]: quadratic coefficient
  // par[1]: linear coefficient
  // par[2]: constant term
  
  Double_t a  = par[0];
  Double_t b  = par[1];
  Double_t c  = par[2];
    
  return a*x[0]*x[0] + b*x[0] + c;
}

//--------------------------------------------------------------------------------------------------
Double_t dMean(const TF1 *fcn, const Double_t x, const TFitResultPtr fs) {
  Double_t df[2];
  df[0] = 1;
  df[1] = x;
  Double_t err2 = df[0]*df[0]*(fs->GetCovarianceMatrix()[0][0]) 
                  + df[1]*df[1]*(fs->GetCovarianceMatrix()[1][1]) 
		  + 2.0*df[0]*df[1]*(fs->GetCovarianceMatrix()[0][1]);
  assert(err2>=0);
  return sqrt(err2);
}

Double_t dMean(const TF1 *fcn, const Double_t x, const Double_t cov[2][2]) {
  Double_t df[2];
  df[0] = 1;
  df[1] = x;
  Double_t err2 = df[0]*df[0]*(cov[0][0]) 
                  + df[1]*df[1]*(cov[1][1]) 
		  + 2.0*df[0]*df[1]*(cov[0][1]);
  assert(err2>=0);
  return sqrt(err2);
}

//--------------------------------------------------------------------------------------------------
Double_t dSigma(const TF1 *fcn, const Double_t x, const TFitResultPtr fs) {
  Double_t df[3];
  Double_t a  = fcn->GetParameter(0);
  Double_t b  = fcn->GetParameter(1);
  Double_t c  = fcn->GetParameter(2);
  
  df[0] = x*x;
  df[1] = x;
  df[2] = 1;
  
  Double_t err2=0;
  for(Int_t i=0; i<3; i++) {
    err2 += df[i]*df[i]*(fs->GetCovarianceMatrix()[i][i]);
    for(Int_t j=i+1; j<3; j++) {
      err2 += 2.0*df[i]*df[j]*(fs->GetCovarianceMatrix()[i][j]);
    }
  }
  assert(err2>=0);
  return sqrt(err2);
}

Double_t dSigma(const TF1 *fcn, const Double_t x, const Double_t cov[3][3]) {
  Double_t df[3];
  Double_t a  = fcn->GetParameter(0);
  Double_t b  = fcn->GetParameter(1);
  Double_t c  = fcn->GetParameter(2);
  
  df[0] = x*x;
  df[1] = x;
  df[2] = 1;
  
  Double_t err2=0;
  for(Int_t i=0; i<3; i++) {
    err2 += df[i]*df[i]*(cov[i][i]);
    for(Int_t j=i+1; j<3; j++) {
      err2 += 2.0*df[i]*df[j]*(cov[i][j]);
    }
  }
  assert(err2>=0);
  return sqrt(err2);
}

//==================================================================================================

RecoilCorrector::RecoilCorrector(TString fname, TString fname_Wp, TString fname_Wm, TString fname_Z, Int_t seed)
{
  gRandom->SetSeed(seed);

  TFile infile(fname);

  fcnPFu1mean  = (TF1*)(infile.Get("fcnPFu1mean")->Clone("fcnPFu1mean"));
  fcnPFu2mean  = (TF1*)(infile.Get("fcnPFu2mean")->Clone("fcnPFu2mean"));
        
  TF1* tmp;
  Double_t xmin,xmax;

  //
  // Load PF-recoil functions
  //
  tmp = (TF1*)infile.Get("fcnPFu1sigma1");
  tmp->GetRange(xmin,xmax);
  fcnPFu1sigma1 = new TF1("fcnPFu1sigma1",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
  for(int i=0; i<tmp->GetNpar(); i++) 
    fcnPFu1sigma1->SetParameter(i,tmp->GetParameter(i));
  
  tmp = (TF1*)infile.Get("fcnPFu1sigma2");
  tmp->GetRange(xmin,xmax);
  fcnPFu1sigma2 = new TF1("fcnPFu1sigma2",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
  for(int i=0; i<tmp->GetNpar(); i++) 
    fcnPFu1sigma2->SetParameter(i,tmp->GetParameter(i));
  
  tmp = (TF1*)infile.Get("fcnPFu1sigma0");
  tmp->GetRange(xmin,xmax);
  fcnPFu1sigma0 = new TF1("fcnPFu1sigma0",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
  for(int i=0; i<tmp->GetNpar(); i++) 
    fcnPFu1sigma0->SetParameter(i,tmp->GetParameter(i));
  
  tmp = (TF1*)infile.Get("fcnPFu2sigma1");
  tmp->GetRange(xmin,xmax);
  fcnPFu2sigma1 = new TF1("fcnPFu2sigma1",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
  for(int i=0; i<tmp->GetNpar(); i++) 
    fcnPFu2sigma1->SetParameter(i,tmp->GetParameter(i));
  
  tmp = (TF1*)infile.Get("fcnPFu2sigma2");
  tmp->GetRange(xmin,xmax);
  fcnPFu2sigma2 = new TF1("fcnPFu2sigma2",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
  for(int i=0; i<tmp->GetNpar(); i++) 
    fcnPFu2sigma2->SetParameter(i,tmp->GetParameter(i));
  
  tmp = (TF1*)infile.Get("fcnPFu2sigma0");
  tmp->GetRange(xmin,xmax);
  fcnPFu2sigma0 = new TF1("fcnPFu2sigma0",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
  for(int i=0; i<tmp->GetNpar(); i++) 
    fcnPFu2sigma0->SetParameter(i,tmp->GetParameter(i));        

  fitresPFu1mean  = (TFitResult*)(infile.Get("fitresPFu1mean")->Clone("fitresPFu1mean"));    
  fitresPFu2mean  = (TFitResult*)(infile.Get("fitresPFu2mean")->Clone("fitresPFu2mean"));
  // (!) explicitly storing covariance matrix, something was screwy with the TFitResult objects
  for(int i=0; i<2; i++) {
    for(int j=0; j<2; j++) {
      pfu1meanCov[i][j] = fitresPFu1mean->GetCovarianceMatrix()[i][j];
      pfu2meanCov[i][j] = fitresPFu2mean->GetCovarianceMatrix()[i][j];
    }
  }

  fitresPFu1sigma0  = (TFitResult*)(infile.Get("fitresPFu1sigma0")->Clone("fitresPFu1sigma0"));
  fitresPFu2sigma0  = (TFitResult*)(infile.Get("fitresPFu2sigma0")->Clone("fitresPFu2sigma0"));
  
  fitresPFu1sigma1  = (TFitResult*)(infile.Get("fitresPFu1sigma1")->Clone("fitresPFu1sigma1"));
  fitresPFu2sigma1  = (TFitResult*)(infile.Get("fitresPFu2sigma1")->Clone("fitresPFu2sigma1"));
  
  fitresPFu1sigma2  = (TFitResult*)(infile.Get("fitresPFu1sigma2")->Clone("fitresPFu1sigma2"));
  fitresPFu2sigma2  = (TFitResult*)(infile.Get("fitresPFu2sigma2")->Clone("fitresPFu2sigma2"));
  // (!) explicitly storing covariance matrix, something was screwy with the TFitResult objects
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      pfu1sigma0Cov[i][j] = fitresPFu1sigma0->GetCovarianceMatrix()[i][j];
      pfu2sigma0Cov[i][j] = fitresPFu2sigma0->GetCovarianceMatrix()[i][j];
      pfu1sigma1Cov[i][j] = fitresPFu1sigma1->GetCovarianceMatrix()[i][j];
      pfu2sigma1Cov[i][j] = fitresPFu2sigma1->GetCovarianceMatrix()[i][j];
      pfu1sigma2Cov[i][j] = fitresPFu1sigma2->GetCovarianceMatrix()[i][j];
      pfu2sigma2Cov[i][j] = fitresPFu2sigma2->GetCovarianceMatrix()[i][j];
    }
  }
  
  infile.Close();


  fcnPFu1mean_Wp=0, fcnPFu1sigma1_Wp=0, fcnPFu1sigma2_Wp=0, fcnPFu1sigma0_Wp=0;
  fcnPFu2mean_Wp=0, fcnPFu2sigma1_Wp=0, fcnPFu2sigma2_Wp=0, fcnPFu2sigma0_Wp=0;
  if(fname_Wp.Length()>0) {
    TFile infile_Wp(fname_Wp);

    fcnPFu1mean_Wp  = (TF1*)(infile_Wp.Get("fcnPFu1mean")->Clone("fcnPFu1mean"));
    fcnPFu2mean_Wp  = (TF1*)(infile_Wp.Get("fcnPFu2mean")->Clone("fcnPFu2mean"));

    //
    // Load PF-recoil functions (W+)
    //
    tmp = (TF1*)infile_Wp.Get("fcnPFu1sigma1");
    tmp->GetRange(xmin,xmax);
    fcnPFu1sigma1_Wp = new TF1("fcnPFu1sigma1",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
    for(int i=0; i<tmp->GetNpar(); i++) 
      fcnPFu1sigma1_Wp->SetParameter(i,tmp->GetParameter(i));
  
    tmp = (TF1*)infile_Wp.Get("fcnPFu1sigma2");
    tmp->GetRange(xmin,xmax);
    fcnPFu1sigma2_Wp = new TF1("fcnPFu1sigma2",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
    for(int i=0; i<tmp->GetNpar(); i++) 
      fcnPFu1sigma2_Wp->SetParameter(i,tmp->GetParameter(i));
  
    tmp = (TF1*)infile_Wp.Get("fcnPFu1sigma0");
    tmp->GetRange(xmin,xmax);
    fcnPFu1sigma0_Wp = new TF1("fcnPFu1sigma0",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
    for(int i=0; i<tmp->GetNpar(); i++) 
      fcnPFu1sigma0_Wp->SetParameter(i,tmp->GetParameter(i));
  
    tmp = (TF1*)infile_Wp.Get("fcnPFu2sigma1");
    tmp->GetRange(xmin,xmax);
    fcnPFu2sigma1_Wp = new TF1("fcnPFu2sigma1",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
    for(int i=0; i<tmp->GetNpar(); i++) 
      fcnPFu2sigma1_Wp->SetParameter(i,tmp->GetParameter(i));
  
    tmp = (TF1*)infile_Wp.Get("fcnPFu2sigma2");
    tmp->GetRange(xmin,xmax);
    fcnPFu2sigma2_Wp = new TF1("fcnPFu2sigma2",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
    for(int i=0; i<tmp->GetNpar(); i++) 
      fcnPFu2sigma2_Wp->SetParameter(i,tmp->GetParameter(i));
  
    tmp = (TF1*)infile_Wp.Get("fcnPFu2sigma0");
    tmp->GetRange(xmin,xmax);
    fcnPFu2sigma0_Wp = new TF1("fcnPFu2sigma0",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
    for(int i=0; i<tmp->GetNpar(); i++) 
      fcnPFu2sigma0_Wp->SetParameter(i,tmp->GetParameter(i));   

    infile_Wp.Close();
  }


  fcnPFu1mean_Wm=0, fcnPFu1sigma1_Wm=0, fcnPFu1sigma2_Wm=0, fcnPFu1sigma0_Wm=0;
  fcnPFu2mean_Wm=0, fcnPFu2sigma1_Wm=0, fcnPFu2sigma2_Wm=0, fcnPFu2sigma0_Wm=0;
  if(fname_Wm.Length()>0) {
    TFile infile_Wm(fname_Wm);

    fcnPFu1mean_Wm  = (TF1*)(infile_Wm.Get("fcnPFu1mean")->Clone("fcnPFu1mean"));
    fcnPFu2mean_Wm  = (TF1*)(infile_Wm.Get("fcnPFu2mean")->Clone("fcnPFu2mean"));

    //
    // Load PF-recoil functions (W-)
    //
    tmp = (TF1*)infile_Wm.Get("fcnPFu1sigma1");
    tmp->GetRange(xmin,xmax);
    fcnPFu1sigma1_Wm = new TF1("fcnPFu1sigma1",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
    for(int i=0; i<tmp->GetNpar(); i++) 
      fcnPFu1sigma1_Wm->SetParameter(i,tmp->GetParameter(i));
  
    tmp = (TF1*)infile_Wm.Get("fcnPFu1sigma2");
    tmp->GetRange(xmin,xmax);
    fcnPFu1sigma2_Wm = new TF1("fcnPFu1sigma2",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
    for(int i=0; i<tmp->GetNpar(); i++) 
      fcnPFu1sigma2_Wm->SetParameter(i,tmp->GetParameter(i));
  
    tmp = (TF1*)infile_Wm.Get("fcnPFu1sigma0");
    tmp->GetRange(xmin,xmax);
    fcnPFu1sigma0_Wm = new TF1("fcnPFu1sigma0",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
    for(int i=0; i<tmp->GetNpar(); i++) 
      fcnPFu1sigma0_Wm->SetParameter(i,tmp->GetParameter(i));
  
    tmp = (TF1*)infile_Wm.Get("fcnPFu2sigma1");
    tmp->GetRange(xmin,xmax);
    fcnPFu2sigma1_Wm = new TF1("fcnPFu2sigma1",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
    for(int i=0; i<tmp->GetNpar(); i++) 
      fcnPFu2sigma1_Wm->SetParameter(i,tmp->GetParameter(i));
  
    tmp = (TF1*)infile_Wm.Get("fcnPFu2sigma2");
    tmp->GetRange(xmin,xmax);
    fcnPFu2sigma2_Wm = new TF1("fcnPFu2sigma2",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
    for(int i=0; i<tmp->GetNpar(); i++) 
      fcnPFu2sigma2_Wm->SetParameter(i,tmp->GetParameter(i));
  
    tmp = (TF1*)infile_Wm.Get("fcnPFu2sigma0");
    tmp->GetRange(xmin,xmax);
    fcnPFu2sigma0_Wm = new TF1("fcnPFu2sigma0",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
    for(int i=0; i<tmp->GetNpar(); i++) 
      fcnPFu2sigma0_Wm->SetParameter(i,tmp->GetParameter(i)); 
  
    infile_Wm.Close();
  }
  

  fcnPFu1mean_Z=0, fcnPFu1sigma1_Z=0, fcnPFu1sigma2_Z=0, fcnPFu1sigma0_Z=0;
  fcnPFu2mean_Z=0, fcnPFu2sigma1_Z=0, fcnPFu2sigma2_Z=0, fcnPFu2sigma0_Z=0;
  if(fname_Z.Length()>0) {
    TFile infile_Z(fname_Z);

    fcnPFu1mean_Z  = (TF1*)(infile_Z.Get("fcnPFu1mean")->Clone("fcnPFu1mean"));
    fcnPFu2mean_Z  = (TF1*)(infile_Z.Get("fcnPFu2mean")->Clone("fcnPFu2mean"));

    //
    // Load PF-recoil functions (Z)
    //
    tmp = (TF1*)infile_Z.Get("fcnPFu1sigma1");
    tmp->GetRange(xmin,xmax);
    fcnPFu1sigma1_Z = new TF1("fcnPFu1sigma1",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
    for(int i=0; i<tmp->GetNpar(); i++) 
      fcnPFu1sigma1_Z->SetParameter(i,tmp->GetParameter(i));
  
    tmp = (TF1*)infile_Z.Get("fcnPFu1sigma2");
    tmp->GetRange(xmin,xmax);
    fcnPFu1sigma2_Z = new TF1("fcnPFu1sigma2",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
    for(int i=0; i<tmp->GetNpar(); i++) 
      fcnPFu1sigma2_Z->SetParameter(i,tmp->GetParameter(i));
  
    tmp = (TF1*)infile_Z.Get("fcnPFu1sigma0");
    tmp->GetRange(xmin,xmax);
    fcnPFu1sigma0_Z = new TF1("fcnPFu1sigma0",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
    for(int i=0; i<tmp->GetNpar(); i++) 
      fcnPFu1sigma0_Z->SetParameter(i,tmp->GetParameter(i));
  
    tmp = (TF1*)infile_Z.Get("fcnPFu2sigma1");
    tmp->GetRange(xmin,xmax);
    fcnPFu2sigma1_Z = new TF1("fcnPFu2sigma1",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
    for(int i=0; i<tmp->GetNpar(); i++) 
      fcnPFu2sigma1_Z->SetParameter(i,tmp->GetParameter(i));
  
    tmp = (TF1*)infile_Z.Get("fcnPFu2sigma2");
    tmp->GetRange(xmin,xmax);
    fcnPFu2sigma2_Z = new TF1("fcnPFu2sigma2",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
    for(int i=0; i<tmp->GetNpar(); i++) 
      fcnPFu2sigma2_Z->SetParameter(i,tmp->GetParameter(i));
  
    tmp = (TF1*)infile_Z.Get("fcnPFu2sigma0");
    tmp->GetRange(xmin,xmax);
    fcnPFu2sigma0_Z = new TF1("fcnPFu2sigma0",sigmaFunc,xmin,xmax,tmp->GetNpar()); 
    for(int i=0; i<tmp->GetNpar(); i++) 
      fcnPFu2sigma0_Z->SetParameter(i,tmp->GetParameter(i)); 
  
    infile_Z.Close();
  }
}

//--------------------------------------------------------------------------------------------------
RecoilCorrector::~RecoilCorrector()
{
  delete fcnPFu1mean;
  delete fcnPFu1sigma1;
  delete fcnPFu1sigma2;
  delete fcnPFu1sigma0;
  delete fcnPFu2mean;
  delete fcnPFu2sigma1;
  delete fcnPFu2sigma2;
  delete fcnPFu2sigma0;

  delete fitresPFu1mean;
  delete fitresPFu1sigma1;
  delete fitresPFu1sigma2;
  delete fitresPFu1sigma0;
  delete fitresPFu2mean;
  delete fitresPFu2sigma1;
  delete fitresPFu2sigma2;
  delete fitresPFu2sigma0;

  delete fcnPFu1mean_Wp;
  delete fcnPFu1sigma1_Wp;
  delete fcnPFu1sigma2_Wp;
  delete fcnPFu1sigma0_Wp;
  delete fcnPFu2mean_Wp;
  delete fcnPFu2sigma1_Wp;
  delete fcnPFu2sigma2_Wp;
  delete fcnPFu2sigma0_Wp;

  delete fcnPFu1mean_Wm;
  delete fcnPFu1sigma1_Wm;
  delete fcnPFu1sigma2_Wm;
  delete fcnPFu1sigma0_Wm;
  delete fcnPFu2mean_Wm;
  delete fcnPFu2sigma1_Wm;
  delete fcnPFu2sigma2_Wm;
  delete fcnPFu2sigma0_Wm;

  delete fcnPFu1mean_Z;
  delete fcnPFu1sigma1_Z;
  delete fcnPFu1sigma2_Z;
  delete fcnPFu1sigma0_Z;
  delete fcnPFu2mean_Z;
  delete fcnPFu2sigma1_Z;
  delete fcnPFu2sigma2_Z;
  delete fcnPFu2sigma0_Z;
}

//--------------------------------------------------------------------------------------------------
void RecoilCorrector::Correct(Double_t &pfmet, Double_t &pfmetphi,
                              Double_t genWPt, Double_t genWPhi,
			      Double_t lepPt,  Double_t lepPhi,
			      Double_t nsigma, Int_t charge
) {  
  
  //
  // Model for PF u1
  //
  Double_t pfu1mean   = fcnPFu1mean->Eval(genWPt);
  Double_t pfu1sigma1 = fcnPFu1sigma1->Eval(genWPt);
  Double_t pfu1sigma2 = fcnPFu1sigma2->Eval(genWPt);
  Double_t pfu1sigma0 = fcnPFu1sigma0->Eval(genWPt);
  if(nsigma!=0) {
    pfu1mean   += nsigma*dMean(fcnPFu1mean,genWPt,pfu1meanCov);//fitresPFu1mean);
    pfu1sigma1 += nsigma*dSigma(fcnPFu1sigma1,genWPt,pfu1sigma1Cov);//fitresPFu1sigma1);
    pfu1sigma2 += nsigma*dSigma(fcnPFu1sigma2,genWPt,pfu1sigma2Cov);//fitresPFu1sigma2);
    pfu1sigma0 += nsigma*dSigma(fcnPFu1sigma0,genWPt,pfu1sigma0Cov);//fitresPFu1sigma0);
  }
  
  //
  // Model for PF u2
  //
  Double_t pfu2mean   = fcnPFu2mean->Eval(genWPt);
  Double_t pfu2sigma1 = fcnPFu2sigma1->Eval(genWPt);
  Double_t pfu2sigma2 = fcnPFu2sigma2->Eval(genWPt);
  Double_t pfu2sigma0 = fcnPFu2sigma0->Eval(genWPt);
  if(nsigma!=0) {
    pfu2mean   += nsigma*dMean(fcnPFu2mean,genWPt,pfu2meanCov);//fitresPFu2mean);
    pfu2sigma1 += nsigma*dSigma(fcnPFu2sigma1,genWPt,pfu2sigma1Cov);//fitresPFu2sigma1);
    pfu2sigma2 += nsigma*dSigma(fcnPFu2sigma2,genWPt,pfu2sigma2Cov);//fitresPFu2sigma2);
    pfu2sigma0 += nsigma*dSigma(fcnPFu2sigma0,genWPt,pfu2sigma0Cov);//fitresPFu2sigma0);
  }
  
  //
  // Apply W/Z corrections if available
  //
  if(fcnPFu1mean_Wp && fcnPFu1mean_Wm && fcnPFu1mean_Z) {
    pfu1mean   *= (charge>0) ? fcnPFu1mean_Wp->Eval(genWPt) / fcnPFu1mean_Z->Eval(genWPt)     : fcnPFu1mean_Wm->Eval(genWPt) / fcnPFu1mean_Z->Eval(genWPt);
    pfu1sigma1 *= (charge>0) ? fcnPFu1sigma1_Wp->Eval(genWPt) / fcnPFu1sigma1_Z->Eval(genWPt) : fcnPFu1sigma1_Wm->Eval(genWPt) / fcnPFu1sigma1_Z->Eval(genWPt);
    pfu1sigma2 *= (charge>0) ? fcnPFu1sigma2_Wp->Eval(genWPt) / fcnPFu1sigma2_Z->Eval(genWPt) : fcnPFu1sigma2_Wm->Eval(genWPt) / fcnPFu1sigma2_Z->Eval(genWPt);
    pfu1sigma0 *= (charge>0) ? fcnPFu1sigma0_Wp->Eval(genWPt) / fcnPFu1sigma0_Z->Eval(genWPt) : fcnPFu1sigma0_Wm->Eval(genWPt) / fcnPFu1sigma0_Z->Eval(genWPt);
    
    pfu2mean   *= (charge>0) ? fcnPFu2mean_Wp->Eval(genWPt) / fcnPFu2mean_Z->Eval(genWPt)     : fcnPFu2mean_Wm->Eval(genWPt) / fcnPFu2mean_Z->Eval(genWPt);
    pfu2sigma1 *= (charge>0) ? fcnPFu2sigma1_Wp->Eval(genWPt) / fcnPFu2sigma1_Z->Eval(genWPt) : fcnPFu2sigma1_Wm->Eval(genWPt) / fcnPFu2sigma1_Z->Eval(genWPt);
    pfu2sigma2 *= (charge>0) ? fcnPFu2sigma2_Wp->Eval(genWPt) / fcnPFu2sigma2_Z->Eval(genWPt) : fcnPFu2sigma2_Wm->Eval(genWPt) / fcnPFu2sigma2_Z->Eval(genWPt);
    pfu2sigma0 *= (charge>0) ? fcnPFu2sigma0_Wp->Eval(genWPt) / fcnPFu2sigma0_Z->Eval(genWPt) : fcnPFu2sigma0_Wm->Eval(genWPt) / fcnPFu2sigma0_Z->Eval(genWPt);
  }
  
  Double_t pfu1frac2  = (pfu1sigma0 - pfu1sigma1)/(pfu1sigma2 - pfu1sigma1);
  Double_t pfu2frac2  = (pfu2sigma0 - pfu2sigma1)/(pfu2sigma2 - pfu2sigma1);
  
  
  Double_t z1 = gRandom->Gaus(0,1);
  Double_t z2 = gRandom->Gaus(0,1);

  Double_t pfu1 = (gRandom->Uniform(0,1) < pfu1frac2) ? z1*pfu1sigma2+pfu1mean : z1*pfu1sigma1+pfu1mean;
  Double_t pfu2 = (gRandom->Uniform(0,1) < pfu2frac2) ? z2*pfu2sigma2+pfu2mean : z2*pfu2sigma1+pfu2mean;
 
  TVector2 vpfmet(-pfu1*cos(genWPhi)+pfu2*sin(genWPhi)-lepPt*cos(lepPhi), -pfu1*sin(genWPhi)-pfu2*cos(genWPhi)-lepPt*sin(lepPhi));  
  pfmet    = vpfmet.Mod();
  pfmetphi = TVector2::Phi_mpi_pi(vpfmet.Phi());
}

