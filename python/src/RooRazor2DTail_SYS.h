//---------------------------------------------------------------------------
#ifndef ROO_Razor2DTail_SYS
#define ROO_Razor2DTail_SYS
//---------------------------------------------------------------------------
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
//---------------------------------------------------------------------------
class RooRealVar;
class RooAbsReal;

#include "Riostream.h"
#include "TMath.h"
#include "Math/SpecFuncMathCore.h"
#include "Math/SpecFuncMathMore.h"

//---------------------------------------------------------------------------
class RooRazor2DTail_SYS : public RooAbsPdf
{
public:
   RooRazor2DTail_SYS() {} ;
   RooRazor2DTail_SYS(const char *name, const char *title,
                  RooAbsReal &_x, RooAbsReal &_y, 
                  RooAbsReal &_x0, RooAbsReal &_y0,
                      RooAbsReal &_b, RooAbsReal &_n);
   RooRazor2DTail_SYS(const RooRazor2DTail_SYS& other,
      const char* name = 0);
   virtual TObject* clone(const char* newname) const { return new RooRazor2DTail_SYS(*this,newname); }
   inline virtual ~RooRazor2DTail_SYS() { }

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

   RooRealProxy X;        // dependent variable
   RooRealProxy Y;        // dependent variable
   RooRealProxy X0;       // X offset
   RooRealProxy Y0;       // Y offset
   RooRealProxy B;        // shape parameter
   RooRealProxy N;        // shape parameter

   Double_t evaluate() const;
private:
  ClassDef(RooRazor2DTail_SYS,1) // Razor2DTail_SYS function
    
};
//---------------------------------------------------------------------------
#endif
