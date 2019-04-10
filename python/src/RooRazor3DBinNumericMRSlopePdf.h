//---------------------------------------------------------------------------
#ifndef ROO_Razor3DBinNumericMRSlopePdf
#define ROO_Razor3DBinNumericMRSlopePdf
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
#include "Math/Functor.h"
#include "Math/WrappedFunction.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"

//---------------------------------------------------------------------------
class RooRazor3DBinNumericMRSlopePdf : public RooAbsPdf
{
public:
   RooRazor3DBinNumericMRSlopePdf() {} ;
   RooRazor3DBinNumericMRSlopePdf(const char *name, const char *title,
		    RooAbsReal& _th1x,
		    RooAbsReal& _x0, RooAbsReal& _y0,
		    RooAbsReal& _b, RooAbsReal& _n,
		    RooAbsReal& _x1,
		    RooAbsReal& _xCut, RooAbsReal& _yCut, RooAbsReal& _zCut);
     //TH3* _Hnominal);
   RooRazor3DBinNumericMRSlopePdf(const RooRazor3DBinNumericMRSlopePdf& other,
      const char* name = 0);
   void setTH3Binning(TH3* _Hnominal);
   void setAbsTol(double _absTol);
   void setRelTol(double _relTol);
   virtual TObject* clone(const char* newname) const { return new RooRazor3DBinNumericMRSlopePdf(*this,newname); }
   inline virtual ~RooRazor3DBinNumericMRSlopePdf() { }

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
   Int_t xBins;        // X bins
   Int_t yBins;        // Y bins
   Int_t zBins;        // Z bins
   Double_t xArray[21]; // xArray[xBins+1]
   Double_t yArray[21]; // yArray[yBins+1]
   Double_t zArray[5]; // zArray[zBins+1]
   Double_t xMax;        // X max
   Double_t yMax;        // Y max
   Double_t zMax;        // Z max
   Double_t xMin;        // X min
   Double_t yMin;        // Y min
   Double_t zMin;        // Z min
   Double_t relTol;      //relative tolerance for numerical integration
   Double_t absTol;      //absolute tolerance for numerical integration

   Double_t evaluate() const;
private:
  ClassDef(RooRazor3DBinNumericMRSlopePdf,1) // Razor3DBinNumericMRSlopePdf function
    
};
//---------------------------------------------------------------------------
#endif

#include "Math/IFunction.h"
#include "Math/IParamFunction.h"
 
class RazorFunctionSlope: public ROOT::Math::IParametricFunctionOneDim
{
private:
   const double *pars;
 
public:
   double DoEvalPar(double x,const double* p) const
   {     
     //params[0] = X0;    params[1] = Y0;
     //params[2] = B;     params[3] = N;
     //params[4] = X1;    params[5] = xMin;
     //params[6] = yLow;  params[7] = yHigh;
     
     double integral = ( (p[6]-p[1])*exp(p[2]*p[3]*pow(p[7]-p[1],1/p[3])*pow(x-p[0],1/p[3])) - (p[7]-p[1])*exp(p[2]*p[3]*pow(p[6]-p[1],1/p[3])*pow(x-p[0],1/p[3])) )*exp(-p[2]*p[3]*(pow(p[6]-p[1],1/p[3])+pow(p[7]-p[1],1/p[3]))*pow(x-p[0],1/p[3]));
     double suppress = p[4]*(x-p[5]) + 1.0;
     return suppress*integral;
   }
   
   double DoEval(double x) const
   {     
     //params[0] = X0;    params[1] = Y0;
     //params[2] = B;     params[3] = N;
     //params[4] = X1;    params[5] = xMin;
     //params[6] = yLow;  params[7] = yHigh;
     
     double integral = ( (pars[6]-pars[1])*exp(pars[2]*pars[3]*pow(pars[7]-pars[1],1/pars[3])*pow(x-pars[0],1/pars[3])) - (pars[7]-pars[1])*exp(pars[2]*pars[3]*pow(pars[6]-pars[1],1/pars[3])*pow(x-pars[0],1/pars[3])) )*exp(-pars[2]*pars[3]*(pow(pars[6]-pars[1],1/pars[3])+pow(pars[7]-pars[1],1/pars[3]))*pow(x-pars[0],1/pars[3]));
     double suppress = pars[4]*(x-pars[5]) + 1.0;
     return suppress*integral;
   }
 
   ROOT::Math::IBaseFunctionOneDim* Clone() const
   {
      return new RazorFunctionSlope();
   }
 
   const double* Parameters() const
   {
      return pars;
   }
 
   void SetParameters(const double* p)
   {
      pars = p;
   }
 
   unsigned int NPar() const
   {
      return 8;
   }
};
