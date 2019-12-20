//---------------------------------------------------------------------------
#ifndef ROO_Razor3DBinMRSlopeErrorPdf
#define ROO_Razor3DBinMRSlopeErrorPdf
//---------------------------------------------------------------------------
#include "RooAbsPdf.h"
#include "RooConstVar.h"
#include "RooRealProxy.h"
//---------------------------------------------------------------------------
class RooRealVar;
class RooAbsReal;

#include "Riostream.h"
#include "TMath.h"
#include <TH3.h>
#include "Math/SpecFuncMathCore.h"
#include "Math/SpecFuncMathMore.h"

//---------------------------------------------------------------------------
class RooRazor3DBinMRSlopeErrorPdf : public RooAbsPdf
{
public:
   RooRazor3DBinMRSlopeErrorPdf() {} ;
   RooRazor3DBinMRSlopeErrorPdf(const char *name, const char *title,
		    RooAbsReal& _th1x,
		    RooAbsReal& _x0, RooAbsReal& _y0,
		    RooAbsReal& _b, RooAbsReal& _n, RooAbsReal& _x1,
			   RooAbsReal& _xCut, RooAbsReal& _yCut, RooAbsReal& _zCut, RooArgList& _pars);
     //TH3* _Hnominal);
   RooRazor3DBinMRSlopeErrorPdf(const RooRazor3DBinMRSlopeErrorPdf& other,
      const char* name = 0);
   void setTH3Binning(TH3* _Hnominal);
   virtual TObject* clone(const char* newname) const { return new RooRazor3DBinMRSlopeErrorPdf(*this,newname); }
   inline virtual ~RooRazor3DBinMRSlopeErrorPdf() { }

   Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
   Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:

   Double_t Chop(const Double_t x) const{
           return (TMath::Abs(x - 0) < 1e-30) ? TMath::Sign(0.0,x) : x;
   }
   Double_t Power(const Double_t x, const Double_t y) const{
	   return Chop(TMath::Power(x,y));
   }
   Double_t Gamma(const Double_t a, const Double_t x) const{
	   return TMath::Gamma(a)*ROOT::Math::inc_gamma_c(a,x);
   }
   Double_t ExpIntegralEi(const Double_t z) const{
	   return Chop(ROOT::Math::expint(z));
   }

   Double_t Gfun(const Double_t x, const Double_t y) const{
     //std::cout << "Gamma(N=" << N << ",BN[(x-X0)(y-Y0)]^(1/N)=" << B*N*pow(x-X0,1/N)*pow(y-Y0,1/N) <<") = " <<  Gamma(N,B*N*pow(x-X0,1/N)*pow(y-Y0,1/N)) << std::endl;
     return Gamma(N,B*N*pow((x-X0)*(y-Y0),1/N));
   }

   RooRealProxy th1x;        // dependent variable
   RooRealProxy X0;       // X offset
   RooRealProxy Y0;       // Y offset
   RooRealProxy B;        // shape parameter
   RooRealProxy N;        // shape parameter
   RooRealProxy X1;        // X slope parameter
   RooRealProxy xCut;        // X cut constant (set pdf to 0 for X < Xcut && Y < Ycut) 
   RooRealProxy yCut;        // Y cut constant (set pdf to 0 for X < Xcut && Y < Ycut) 
   RooRealProxy zCut;        // Z cut constant (set pdf to 0 unless Zut <= Z < Zcut)
   RooListProxy pars;
   Int_t xBins;        // X bins
   Int_t yBins;        // Y bins
   Int_t zBins;        // Z bins
   Double_t xArray[21]; // xArray[xBins+1]
   Double_t yArray[21]; // yArray[yBins+1]
   Double_t zArray[5]; // zArray[zBins+1]
   Double_t histArray[1600]; // histArray[xBins*yBins*zBins]
   //Double_t xArray[21]; // xArray[xBins+1]
   //Double_t yArray[21]; // yArray[yBins+1]
   //Double_t zArray[5]; // zArray[zBins+1]
   //Double_t histArray[1600]; // histArray[xBins*yBins*zBins]
   Double_t xMax;        // X max
   Double_t yMax;        // Y max
   Double_t zMax;        // Z max
   Double_t xMin;        // X min
   Double_t yMin;        // Y min
   Double_t zMin;        // Z min

   Double_t evaluate() const;
private:
  ClassDef(RooRazor3DBinMRSlopeErrorPdf,1) // Razor2DTail_SYS function
    
};
//---------------------------------------------------------------------------
#endif
