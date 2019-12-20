//---------------------------------------------------------------------------
#include "RooFit.h"

#include "Riostream.h"
#include <TMath.h>
#include <cassert>
#include <cmath>
#include <math.h>

#include "RooRazor3DBinMRSlopeErrorPdf.h"
#include "RooRealVar.h"

using namespace std;

ClassImp(RooRazor3DBinMRSlopeErrorPdf)
//---------------------------------------------------------------------------
RooRazor3DBinMRSlopeErrorPdf::RooRazor3DBinMRSlopeErrorPdf(const char *name, const char *title,
				   RooAbsReal& _th1x,  
				   RooAbsReal& _x0, RooAbsReal& _y0, 
				   RooAbsReal& _b, RooAbsReal& _n, RooAbsReal& _x1, 
						 RooAbsReal& _xCut, RooAbsReal& _yCut, RooAbsReal& _zCut, RooArgList& _pars) : RooAbsPdf(name, title), 
//TH3* _Hnominal) : RooAbsPdf(name, title), 
  th1x("th1x", "th1x Observable", this, _th1x),
  X0("X0", "X Offset", this, _x0),
  Y0("Y0", "Y Offset", this, _y0),
  B("B", "B Shape parameter", this, _b),
  N("N", "N Shape parameter", this, _n),
  X1("X1", "X Slope parameter", this, _x1),
  xCut("xCut", "X Cut parameter",this, _xCut),
  yCut("yCut", "Y Cut parameter",this, _yCut),
  zCut("zCut", "Z Cut parameter",this, _zCut),
  pars("pars","pars",this),
  xBins(0),
  yBins(0),
  zBins(0),
  xMax(0),
  yMax(0),
  zMax(0),
  xMin(0),
  yMin(0),
  zMin(0)
{
  memset(&xArray, 0, sizeof(xArray));
  memset(&yArray, 0, sizeof(yArray));
  memset(&zArray, 0, sizeof(zArray));
  memset(&histArray, 0, sizeof(histArray));
  
  TIterator *varIter=_pars.createIterator(); 
  RooAbsReal *fVar;
  while ( (fVar = (RooAbsReal*)varIter->Next()) ){
	pars.add(*fVar);
  }
}
//---------------------------------------------------------------------------
RooRazor3DBinMRSlopeErrorPdf::RooRazor3DBinMRSlopeErrorPdf(const RooRazor3DBinMRSlopeErrorPdf& other, const char* name) :
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
   pars("_pars",this,RooListProxy()),
   xBins(other.xBins),
   yBins(other.yBins),
   zBins(other.zBins),
   xMax(other.xMax),
   yMax(other.yMax),
   zMax(other.zMax),
   xMin(other.xMin),
   yMin(other.yMin),
   zMin(other.zMin)
{
  //memset(&xArray, 0, sizeof(xArray));
  //memset(&yArray, 0, sizeof(yArray));
  //memset(&zArray, 0, sizeof(zArray));
  //memset(&histArray, 0, sizeof(histArray));
  for (Int_t i=0; i<xBins+1; i++){
    xArray[i] = other.xArray[i];
  }
  for (Int_t j=0; j<yBins+1; j++){
    yArray[j] =  other.yArray[j];
  }
  for (Int_t k=0; k<zBins+1; k++){
    zArray[k] =  other.zArray[k];
  }

  Int_t th1xBin = 0;
  for (Int_t i=0; i<xBins+1; i++){
    for (Int_t j=0; j<yBins+1; j++){
      for (Int_t k=0; k<zBins+1; k++){
	histArray[th1xBin] =  other.histArray[th1xBin];
	th1xBin++;
      }
    }
  }
  
  TIterator *varIter=other.pars.createIterator(); 
  RooAbsReal *fVar;
  while ( (fVar = (RooAbsReal*) varIter->Next()) ){
	pars.add(*fVar);
  }
}
//---------------------------------------------------------------------------
void RooRazor3DBinMRSlopeErrorPdf::setTH3Binning(TH3* _Hnominal){
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
  memset(&histArray, 0, sizeof(histArray));
  for (Int_t i=0; i<xBins+1; i++){
    xArray[i] =  _Hnominal->GetXaxis()->GetBinLowEdge(i+1);
  }
  for (Int_t j=0; j<yBins+1; j++){
    yArray[j] =  _Hnominal->GetYaxis()->GetBinLowEdge(j+1);
  }
  for (Int_t k=0; k<zBins+1; k++){
    zArray[k] =  _Hnominal->GetZaxis()->GetBinLowEdge(k+1);
  }
  Int_t th1xBin = 0;
  for (Int_t i=0; i<xBins; i++){
    for (Int_t j=0; j<yBins; j++){
      for (Int_t k=0; k<zBins; k++){
	histArray[th1xBin] = _Hnominal->GetBinContent(i+1,j+1,k+1);
	//cout << xArray[i] << ", " << xArray[i+1] << endl;
	//cout << yArray[j] << ", " << yArray[j+1] << endl;
	//cout << zArray[k] << ", " << zArray[k+1] << endl;
	//cout << histArray[th1xBin] << endl;
	th1xBin++;
      }
    }
  }
  
}
//---------------------------------------------------------------------------
Double_t RooRazor3DBinMRSlopeErrorPdf::evaluate() const
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

  RooAbsReal *retVar = (RooAbsReal*)pars.at(iBin);    
  double ret = retVar->getVal();
  
  Int_t zBin = iBin % zBins;
  Int_t yBin = ( (iBin - zBin)/(zBins) ) % (yBins);
  Int_t xBin =  (iBin - zBin - yBin*zBins ) / (zBins*yBins);

  //cout << "in bin " << iBin << " which is in range" << endl;
  //cout << "(" << xBin+1 << "," << yBin+1 << "," << zBin+1 << ")" << endl;

  Double_t zLow = zArray[zBin];
  Double_t zHigh = zArray[zBin+1];
  Double_t suppress = 1.;
  if (zCut >= zLow and zCut < zHigh){
    Double_t xLow = xArray[xBin];
    Double_t xHigh = xArray[xBin+1];
    Double_t yLow = yArray[yBin];
    Double_t yHigh = yArray[yBin+1];
    suppress = X1*((xLow+xHigh)/2.-(xMin+xArray[1])/2.) + 1.0;
    
    total_integral = -Gfun(xMin,yMax)-Gfun(xMax,yMin)+Gfun(xMax,yMax)+Gfun(xMin,yCut)+Gfun(xCut,yMin)-Gfun(xCut,yCut);

    if(xHigh <= xCut && yHigh <= yCut) {
      return 0.0;
    }
    else if(xLow < xCut && xHigh > xCut && yHigh <= yCut) {
      integral = Gfun(xCut,yLow)-Gfun(xCut,yHigh)-Gfun(xHigh,yLow)+Gfun(xHigh,yHigh);
    }
    else if(yLow < yCut && yHigh > yCut && xHigh <= xCut) {
      integral = Gfun(xLow,yCut)-Gfun(xLow,yHigh)-Gfun(xHigh,yCut)+Gfun(xHigh,yHigh);
    }
    else if(xLow < xCut && xHigh > xCut && yLow < yCut && yHigh > yCut) {
      integral = -Gfun(xLow,yHigh)-Gfun(xHigh,yLow)+Gfun(xHigh,yHigh)+Gfun(xLow,yCut)+Gfun(xCut,yLow)-Gfun(xCut,yCut);
    }
    else {
      integral = Gfun(xLow,yLow)-Gfun(xLow,yHigh)-Gfun(xHigh,yLow)+Gfun(xHigh,yHigh);
    }

  }

  if (total_integral>0.0) {
    if (histArray[iBin] > 0)
      return pow(histArray[iBin],ret) * suppress * integral;
    else
      return suppress*integral;
  } else return 0;

}

// //---------------------------------------------------------------------------
Int_t RooRazor3DBinMRSlopeErrorPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const{
  if (matchArgs(allVars, analVars, th1x)) return 1;
  return 0;
}

// //---------------------------------------------------------------------------
Double_t RooRazor3DBinMRSlopeErrorPdf::analyticalIntegral(Int_t code, const char* rangeName) const{

   Double_t th1xMin = th1x.min(rangeName); Double_t th1xMax = th1x.max(rangeName);
   Int_t iBinMin = (Int_t) th1xMin; Int_t iBinMax = (Int_t) th1xMax;

   
   if(B <= 0. || N <= 0. || X0 >= xMin || Y0 >= yMin) return 1.;

   Double_t integral = 0.0;
   Double_t total_integral =  1.0;
      
   //cout <<  "iBinMin = " << iBinMin << ",iBinMax = " << iBinMax << endl;
   Int_t nBins =  xBins*yBins*zBins;

   //   if (code==1 && iBinMin==0 && iBinMax>=nBins){
   //  integral = -Gfun(xMin,yMax)-Gfun(xMax,yMin)+Gfun(xMax,yMax)+Gfun(xMin,yCut)+Gfun(xCut,yMin)-Gfun(xCut,yCut);
   //}
   if(code==1) { 
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
	 
	 RooAbsReal *retVar = (RooAbsReal*)pars.at(iBin);    
	 double ret = retVar->getVal();
	 
	 if (zCut >= zLow and zCut < zHigh){
	   Double_t xLow = xArray[xBin];
	   Double_t xHigh = xArray[xBin+1];
	   Double_t yLow = yArray[yBin];
	   Double_t yHigh = yArray[yBin+1];	   
	   Double_t suppress = X1*((xLow+xHigh)/2.-(xMin+xArray[1])/2.) + 1.0;
	   if(xHigh <= xCut && yHigh <= yCut) integral += 0.0;
	   else if(xLow < xCut && xHigh > xCut && yHigh <= yCut) {	     
	     if (histArray[iBin] > 0)
	       integral += pow(histArray[iBin],ret) * suppress * ( Gfun(xCut,yLow)-Gfun(xCut,yHigh)-Gfun(xHigh,yLow)+Gfun(xHigh,yHigh) );
	     else
	       integral += suppress * ( Gfun(xCut,yLow)-Gfun(xCut,yHigh)-Gfun(xHigh,yLow)+Gfun(xHigh,yHigh) );	       
	   }
	   else if(yLow < yCut && yHigh > yCut && xHigh <= xCut) {     
	     if (histArray[iBin] > 0)
	       integral += pow(histArray[iBin],ret) * suppress * ( Gfun(xLow,yCut)-Gfun(xLow,yHigh)-Gfun(xHigh,yCut)+Gfun(xHigh,yHigh) );
	     else
	       integral += suppress * ( Gfun(xLow,yCut)-Gfun(xLow,yHigh)-Gfun(xHigh,yCut)+Gfun(xHigh,yHigh) );	       
	   }
	   else if(xLow < xCut && xHigh > xCut && yLow < yCut && yHigh > yCut) {
	     if (histArray[iBin] > 0)
	       integral += pow(histArray[iBin],ret) * suppress * ( -Gfun(xLow,yHigh)-Gfun(xHigh,yLow)+Gfun(xHigh,yHigh)+Gfun(xLow,yCut)+Gfun(xCut,yLow)-Gfun(xCut,yCut) );
	     else
	       integral += suppress * ( -Gfun(xLow,yHigh)-Gfun(xHigh,yLow)+Gfun(xHigh,yHigh)+Gfun(xLow,yCut)+Gfun(xCut,yLow)-Gfun(xCut,yCut) );
	   }
	   else {
	     if (histArray[iBin] > 0)
	       integral += pow(histArray[iBin],ret) * suppress * ( Gfun(xLow,yLow)-Gfun(xLow,yHigh)-Gfun(xHigh,yLow)+Gfun(xHigh,yHigh) );
	     else
	       integral += suppress * ( Gfun(xLow,yLow)-Gfun(xLow,yHigh)-Gfun(xHigh,yLow)+Gfun(xHigh,yHigh) );
	   }
	 }
       }
     }
   } else {
     cout << "WARNING IN RooRazor3DBinMRSlopeErrorPdf: integration code is not correct" << endl;
     cout << "                           what are you integrating on?" << endl;
     return 1.0;
   }

   if (total_integral>0.0) {
     
     return integral;
   } else return 1.0;
}
// //---------------------------------------------------------------------------

