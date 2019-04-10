//---------------------------------------------------------------------------
#include "RooFit.h"

#include "Riostream.h"
#include <TMath.h>
#include <cassert>
#include <cmath>
#include <math.h>

#include "RooRazor3DBinNumericMRSlopePdf.h"
#include "RooRealVar.h"
#include "Math/Functor.h"
#include "Math/WrappedFunction.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"
#include "Math/GSLIntegrator.h"

using namespace std;

ClassImp(RooRazor3DBinNumericMRSlopePdf)
//---------------------------------------------------------------------------
RooRazor3DBinNumericMRSlopePdf::RooRazor3DBinNumericMRSlopePdf(const char *name, const char *title,
				   RooAbsReal& _th1x,  
				   RooAbsReal& _x0, RooAbsReal& _y0, 
				   RooAbsReal& _b, RooAbsReal& _n,
				   RooAbsReal& _x1,
				   RooAbsReal& _xCut, RooAbsReal& _yCut, RooAbsReal& _zCut) : RooAbsPdf(name, title), 
//TH3* _Hnominal) : RooAbsPdf(name, title), 
  th1x("th1x", "th1x Observable", this, _th1x),
  X0("X0", "X Offset", this, _x0),
  Y0("Y0", "Y Offset", this, _y0),
  B("B", "B Shape parameter", this, _b),
  N("N", "N Shape parameter", this, _n),
  X1("X1", "X slope parameter", this, _x1),
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
RooRazor3DBinNumericMRSlopePdf::RooRazor3DBinNumericMRSlopePdf(const RooRazor3DBinNumericMRSlopePdf& other, const char* name) :
   RooAbsPdf(other, name), 
   th1x("th1x", this, other.th1x),  
   X0("X0", this, other.X0),
   Y0("Y0", this, other.Y0),
   B("B", this, other.B),
   N("N", this, other.N),
   X1("X1", this, other.X1),
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
void RooRazor3DBinNumericMRSlopePdf::setTH3Binning(TH3* _Hnominal){
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
void RooRazor3DBinNumericMRSlopePdf::setRelTol(double _relTol){
  relTol = _relTol;
}
//---------------------------------------------------------------------------
void RooRazor3DBinNumericMRSlopePdf::setAbsTol(double _absTol){
  absTol = _absTol;
}
//---------------------------------------------------------------------------
Double_t RooRazor3DBinNumericMRSlopePdf::evaluate() const
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
    RazorFunctionSlope func;
    double params[8];
    params[0] = X0;    params[1] = Y0;
    params[2] = B;     params[3] = N;
    params[4] = X1;    params[5] = xMin;
    params[6] = yMin;  params[7] = yMax;
    func.SetParameters(params);
    ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,absTol,relTol);
    ig.SetFunction(func,false);
      
    total_integral = (N/pow(B*N,N))*(-Gfun(xMin,yMax)-Gfun(xMax,yMin)+Gfun(xMax,yMax)+Gfun(xMin,yCut)+Gfun(xCut,yMin)-Gfun(xCut,yCut)) ;

    if(xHigh <= xCut && yHigh <= yCut) {
      return 0.0;
    }
    else if(yLow < yCut && yHigh > yCut && xHigh <= xCut) {
      params[6] = yCut;  params[7] = yHigh;
      integral = ig.Integral(xLow,xHigh);
    }
    else if(xLow < xCut && xHigh > xCut && yHigh <= yCut) {
      params[6] = yLow;  params[7] = yHigh;
      integral = ig.Integral(xCut,xHigh);
    }
    else if(yLow < yCut && yHigh > yCut && xLow < xCut && xHigh > xCut) {
      params[6] = yLow;  params[7] = yHigh;
      integral = ig.Integral(xLow,xHigh);
      params[6] = yLow;  params[7] = yCut;
      integral -= ig.Integral(xLow,xCut);
    }
    else {
      params[6] = yLow;  params[7] = yHigh;
      integral = ig.Integral(xLow,xHigh);
      //cout << endl;
      //cout << "X0 = " << X0 << ", Y0 = " << Y0 << endl;
      //cout << "B = " << B << ", N = " << N << endl;
      //cout << "X1 = " << X1 << ", xMin = " << xMin << endl;
      //cout << "yLow = " << yLow << ", yHigh = " << yHigh << endl;
      //cout << "xLow = " << xLow << ", xHigh = " << xHigh << endl;
      //cout << "integral = " << integral << endl;
    }
      
  }

  if (total_integral>0.0) {
    return integral;
  } else return 0;

}

// //---------------------------------------------------------------------------
Int_t RooRazor3DBinNumericMRSlopePdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const{
  if (matchArgs(allVars, analVars, th1x)) return 1;
  return 0;
}

// //---------------------------------------------------------------------------
Double_t RooRazor3DBinNumericMRSlopePdf::analyticalIntegral(Int_t code, const char* rangeName) const{

   Double_t th1xMin = th1x.min(rangeName); Double_t th1xMax = th1x.max(rangeName);
   Int_t iBinMin = (Int_t) th1xMin; Int_t iBinMax = (Int_t) th1xMax;

   
   if(B <= 0. || N <= 0. || X0 >= xMin || Y0 >= yMin) return 1.;

   Double_t integral = 0.0;
   Double_t total_integral =  1.0;
      
   //cout <<  "iBinMin = " << iBinMin << ",iBinMax = " << iBinMax << endl;
   Int_t nBins =  xBins*yBins*zBins;

   
   // define the function to be integrated numerically
   RazorFunctionSlope func;
   double params[8];
   params[0] = X0;    params[1] = Y0;
   params[2] = B;     params[3] = N;
   params[4] = X1;    params[5] = xMin;
   params[6] = yMin;  params[7] = yMax;
   func.SetParameters(params);
   ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,absTol,relTol);
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
	   else if(yLow < yCut && yHigh > yCut && xHigh <= xCut) {
	     params[6] = yCut;  params[7] = yHigh;
	     integral += ig.Integral(xLow,xHigh);
	   }
	   else if(xLow < xCut && xHigh > xCut && yHigh <= yCut) {
	     params[6] = yLow;  params[7] = yHigh;
	     integral += ig.Integral(xCut,xHigh);
	   }
	   else if(yLow < yCut && yHigh > yCut && xLow < xCut && xHigh > xCut) {
	     params[6] = yLow;  params[7] = yHigh;
	     integral += ig.Integral(xLow,xHigh);
	     params[6] = yLow;  params[7] = yCut;
	     integral -= ig.Integral(xLow,xCut);
	   }
	   else {
	     params[6] = yLow;  params[7] = yHigh;
	     integral += ig.Integral(xLow,xHigh);
	   }
	 }
       }
     }
   } else {
     cout << "WARNING IN RooRazor3DBinNumericMRSlopePdf: integration code is not correct" << endl;
     cout << "                           what are you integrating on?" << endl;
     return 1.0;
   }

   if (total_integral>0.0) {
     
     return integral;
   } else return 1.0;
}
// //---------------------------------------------------------------------------

