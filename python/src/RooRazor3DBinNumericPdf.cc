//---------------------------------------------------------------------------
#include "RooFit.h"

#include "Riostream.h"
#include <TMath.h>
#include <cassert>
#include <cmath>
#include <math.h>

#include "RooRazor3DBinNumericPdf.h"
#include "RooRealVar.h"
#include "Math/Functor.h"
#include "Math/WrappedFunction.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"
#include "Math/GSLIntegrator.h"

using namespace std;

ClassImp(RooRazor3DBinNumericPdf)
//---------------------------------------------------------------------------
RooRazor3DBinNumericPdf::RooRazor3DBinNumericPdf(const char *name, const char *title,
				   RooAbsReal& _th1x,  
				   RooAbsReal& _x0, RooAbsReal& _y0, 
				   RooAbsReal& _b, RooAbsReal& _n,
				   RooAbsReal& _y1, RooAbsReal& _y2,
				   RooAbsReal& _xCut, RooAbsReal& _yCut, RooAbsReal& _zCut) : RooAbsPdf(name, title), 
//TH3* _Hnominal) : RooAbsPdf(name, title), 
  th1x("th1x", "th1x Observable", this, _th1x),
  X0("X0", "X Offset", this, _x0),
  Y0("Y0", "Y Offset", this, _y0),
  B("B", "B Shape parameter", this, _b),
  N("N", "N Shape parameter", this, _n),
  Y1("Y1", "Y turn-off offset", this, _y1),
  Y2("Y2", "Y turn-off width", this, _y2),
  xCut("xCut", "X Cut parameter",this, _xCut),
  yCut("yCut", "Y Cut parameter",this, _yCut),
  zCut("zCut", "Z Cut parameter",this, _zCut),
  xBins(0),
  yBins(0),
  zBins(0),
  xMax(0),
  yMax(0),
  zMax(0),
  xMin(0),
  yMin(0),
  zMin(0),
  relTol(1E-12),
  absTol(1E-12)
{
  memset(&xArray, 0, sizeof(xArray));
  memset(&yArray, 0, sizeof(yArray));
  memset(&zArray, 0, sizeof(zArray));
}
//---------------------------------------------------------------------------
RooRazor3DBinNumericPdf::RooRazor3DBinNumericPdf(const RooRazor3DBinNumericPdf& other, const char* name) :
   RooAbsPdf(other, name), 
   th1x("th1x", this, other.th1x),  
   X0("X0", this, other.X0),
   Y0("Y0", this, other.Y0),
   B("B", this, other.B),
   N("N", this, other.N),
   Y1("Y1", this, other.Y1),
   Y2("Y2", this, other.Y2),
   xCut("xCut", this, other.xCut),
   yCut("yCut", this, other.yCut),
   zCut("zCut", this, other.zCut),
   xBins(other.xBins),
   yBins(other.yBins),
   zBins(other.zBins),
   xMax(other.xMax),
   yMax(other.yMax),
   zMax(other.zMax),
   xMin(other.xMin),
   yMin(other.yMin),
   zMin(other.zMin),
   relTol(other.relTol),
   absTol(other.absTol)
{
  //memset(&xArray, 0, sizeof(xArray));
  //memset(&yArray, 0, sizeof(yArray));
  //memset(&zArray, 0, sizeof(zArray));
  for (Int_t i=0; i<xBins+1; i++){
    xArray[i] = other.xArray[i];
  }
  for (Int_t j=0; j<yBins+1; j++){
    yArray[j] =  other.yArray[j];
  }
  for (Int_t k=0; k<zBins+1; k++){
    zArray[k] =  other.zArray[k];
  }
}
//---------------------------------------------------------------------------
void RooRazor3DBinNumericPdf::setTH3Binning(TH3* _Hnominal){
  xBins = _Hnominal->GetXaxis()->GetNbins();
  yBins = _Hnominal->GetYaxis()->GetNbins();
  zBins = _Hnominal->GetZaxis()->GetNbins();
  xMin = _Hnominal->GetXaxis()->GetBinLowEdge(1);
  yMin = _Hnominal->GetYaxis()->GetBinLowEdge(1);
  zMin = _Hnominal->GetZaxis()->GetBinLowEdge(1);
  xMax = _Hnominal->GetXaxis()->GetBinUpEdge(xBins);
  yMax = _Hnominal->GetYaxis()->GetBinUpEdge(yBins);
  zMax = _Hnominal->GetZaxis()->GetBinUpEdge(zBins);
  memset(&xArray, 0, sizeof(xArray));
  memset(&yArray, 0, sizeof(yArray));
  memset(&zArray, 0, sizeof(zArray));
  for (Int_t i=0; i<xBins+1; i++){
    xArray[i] =  _Hnominal->GetXaxis()->GetBinLowEdge(i+1);
  }
  for (Int_t j=0; j<yBins+1; j++){
    yArray[j] =  _Hnominal->GetYaxis()->GetBinLowEdge(j+1);
  }
  for (Int_t k=0; k<zBins+1; k++){
    zArray[k] =  _Hnominal->GetZaxis()->GetBinLowEdge(k+1);
  }
}
//---------------------------------------------------------------------------
void RooRazor3DBinNumericPdf::setRelTol(double _relTol){
  relTol = _relTol;
}
//---------------------------------------------------------------------------
void RooRazor3DBinNumericPdf::setAbsTol(double _absTol){
  absTol = _absTol;
}
//---------------------------------------------------------------------------
Double_t RooRazor3DBinNumericPdf::evaluate() const
{
  Double_t integral = 0.0;
  Double_t total_integral = 1.0;
  
  if(B <= 0. || N <= 0. || X0 >= xMin || Y0 >= yMin) return 0.0;

  Int_t nBins = xBins*yBins*zBins;

  Int_t iBin = (Int_t) th1x;
  if(iBin < 0 || iBin >= nBins) {
    //cout << "in bin " << iBin << " which is outside of range" << endl;
    return 0.0;
  }

  
  Int_t zBin = iBin % zBins;
  Int_t yBin = ( (iBin - zBin)/(zBins) ) % (yBins);
  Int_t xBin =  (iBin - zBin - yBin*zBins ) / (zBins*yBins);

  //cout << "in bin " << iBin << " which is in range" << endl;
  //cout << "(" << xBin+1 << "," << yBin+1 << "," << zBin+1 << ")" << endl;

  Double_t zLow = zArray[zBin];
  Double_t zHigh = zArray[zBin+1];
  if (zCut >= zLow and zCut < zHigh){
    Double_t xLow = xArray[xBin];
    Double_t xHigh = xArray[xBin+1];
    Double_t yLow = yArray[yBin];
    Double_t yHigh = yArray[yBin+1];

    
    // define the function to be integrated numerically
    RazorFunctionErf func;
    double params[8];
    params[0] = X0;    params[1] = Y0;
    params[2] = B;     params[3] = N;
    params[4] = Y1;    params[5] = Y2;
    func.SetParameters(params);
    ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,absTol,relTol);
    ig.SetFunction(func,false);
      
    total_integral = (N/pow(B*N,N))*(-Gfun(xMin,yMax)-Gfun(xMax,yMin)+Gfun(xMax,yMax)+Gfun(xMin,yCut)+Gfun(xCut,yMin)-Gfun(xCut,yCut)) ;

    if(xHigh <= xCut && yHigh <= yCut) {
      return 0.0;
    }
    else if(xLow < xCut && xHigh > xCut && yHigh <= yCut) {
      //integral = (N/pow(B*N,N))*(Gfun(xCut,yLow)-Gfun(xCut,yHigh)-Gfun(xHigh,yLow)+Gfun(xHigh,yHigh));
      params[6] = xCut;  params[7] = xHigh;
      integral = ig.Integral(yLow,yHigh);
    }
    else if(yLow < yCut && yHigh > yCut && xHigh <= xCut) {
      //integral = (N/pow(B*N,N))*(Gfun(xLow,yCut)-Gfun(xLow,yHigh)-Gfun(xHigh,yCut)+Gfun(xHigh,yHigh));
      params[6] = xLow;  params[7] = xHigh;
      integral = ig.Integral(yCut,yHigh);
    }
    else if(xLow < xCut && xHigh > xCut && yLow < yCut && yHigh > yCut) {
      //integral = (N/pow(B*N,N))*(-Gfun(xLow,yHigh)-Gfun(xHigh,yLow)+Gfun(xHigh,yHigh)+Gfun(xLow,yCut)+Gfun(xCut,yLow)-Gfun(xCut,yCut));
      params[6] = xLow;  params[7] = xHigh;
      integral = ig.Integral(yLow,yHigh);
      params[6] = xLow;  params[7] = xCut;
      integral -= ig.Integral(yLow,yCut);
    }
    else {
      //integral = (N/pow(B*N,N))*(Gfun(xLow,yLow)-Gfun(xLow,yHigh)-Gfun(xHigh,yLow)+Gfun(xHigh,yHigh));
      params[6] = xLow;  params[7] = xHigh;
      integral = ig.Integral(yLow,yHigh);
    }
      
  }

  if (total_integral>0.0) {
    return integral;
  } else return 0;

}

// //---------------------------------------------------------------------------
Int_t RooRazor3DBinNumericPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const{
  if (matchArgs(allVars, analVars, th1x)) return 1;
  return 0;
}

// //---------------------------------------------------------------------------
Double_t RooRazor3DBinNumericPdf::analyticalIntegral(Int_t code, const char* rangeName) const{

   Double_t th1xMin = th1x.min(rangeName); Double_t th1xMax = th1x.max(rangeName);
   Int_t iBinMin = (Int_t) th1xMin; Int_t iBinMax = (Int_t) th1xMax;

   
   if(B <= 0. || N <= 0. || X0 >= xMin || Y0 >= yMin) return 1.;

   Double_t integral = 0.0;
   Double_t total_integral =  1.0;
      
   //cout <<  "iBinMin = " << iBinMin << ",iBinMax = " << iBinMax << endl;
   Int_t nBins =  xBins*yBins*zBins;

   
   // define the function to be integrated numerically
   RazorFunctionErf func;
   double params[8];
   params[0] = X0;    params[1] = Y0;
   params[2] = B;     params[3] = N;
   params[4] = Y1;    params[5] = Y2;
   func.SetParameters(params);
   ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kNONADAPTIVE,absTol,relTol);
   ig.SetFunction(func,false);
    

   if (code==1 && iBinMin<=0 && iBinMax>=nBins){
     integral = (N/pow(B*N,N))*(-Gfun(xMin,yMax)-Gfun(xMax,yMin)+Gfun(xMax,yMax)+Gfun(xMin,yCut)+Gfun(xCut,yMin)-Gfun(xCut,yCut));
     
   }
   else if(code==1) { 
     total_integral = Gfun(xMin,yMin)-Gfun(xMin,yMax)-Gfun(xMax,yMin)+Gfun(xMax,yMax);
     for (Int_t iBin=iBinMin; iBin<iBinMax; iBin++){
       Int_t zBin = iBin % zBins;
       Int_t yBin = ( (iBin - zBin)/(zBins) ) % (yBins);
       Int_t xBin =  (iBin - zBin - yBin*zBins ) / (zBins*yBins);
 
       Double_t zLow = zArray[zBin];
       Double_t zHigh = zArray[zBin+1];
       
       if(iBin < 0 || iBin >= nBins) {
	 integral += 0.0;
       }
       else{
	 if (zCut >= zLow and zCut < zHigh){
	   Double_t xLow = xArray[xBin];
	   Double_t xHigh = xArray[xBin+1];
	   Double_t yLow = yArray[yBin];
	   Double_t yHigh = yArray[yBin+1];
	   if(xHigh <= xCut && yHigh <= yCut) integral += 0.0;
	   else if(xLow < xCut && xHigh > xCut && yHigh <= yCut) {
	     //integral += (N/pow(B*N,N))*(Gfun(xCut,yLow)-Gfun(xCut,yHigh)-Gfun(xHigh,yLow)+Gfun(xHigh,yHigh));
	     params[6] = xCut;  params[7] = xHigh;
	     integral += ig.Integral(yLow,yHigh);
	   }
	   else if(yLow < yCut && yHigh > yCut && xHigh <= xCut) {
	     //integral += (N/pow(B*N,N))*(Gfun(xLow,yCut)-Gfun(xLow,yHigh)-Gfun(xHigh,yCut)+Gfun(xHigh,yHigh));
	     params[6] = xLow;  params[7] = xHigh;
	     integral += ig.Integral(yCut,yHigh);
	   }
	   else if(xLow < xCut && xHigh > xCut && yLow < yCut && yHigh > yCut) {
	     //integral += (N/pow(B*N,N))*(-Gfun(xLow,yHigh)-Gfun(xHigh,yLow)+Gfun(xHigh,yHigh)+Gfun(xLow,yCut)+Gfun(xCut,yLow)-Gfun(xCut,yCut));
	     params[6] = xLow;  params[7] = xHigh;
	     integral += ig.Integral(yLow,yHigh);
	     params[6] = xLow;  params[7] = xCut;
	     integral -= ig.Integral(yLow,yCut);
	   }
	   else {
	     //integral += (N/pow(B*N,N))*(Gfun(xLow,yLow)-Gfun(xLow,yHigh)-Gfun(xHigh,yLow)+Gfun(xHigh,yHigh));
	     params[6] = xLow;  params[7] = xHigh;
	     integral += ig.Integral(yLow,yHigh);
	   }
	 }
       }
     }
   } else {
     cout << "WARNING IN RooRazor3DBinNumericPdf: integration code is not correct" << endl;
     cout << "                           what are you integrating on?" << endl;
     return 1.0;
   }

   if (total_integral>0.0) {
     
     return integral;
   } else return 1.0;
}
// //---------------------------------------------------------------------------

