//---------------------------------------------------------------------------
#ifndef ROO_BTagMult
#define ROO_BTagMult
//---------------------------------------------------------------------------
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
//---------------------------------------------------------------------------
class RooRealVar;
class RooAbsReal;

#include "TMath.h"
#include "Math/SpecFuncMathCore.h"
#include "Math/SpecFuncMathMore.h"

//---------------------------------------------------------------------------
class RooBTagMult : public RooAbsPdf
{
public:
   RooBTagMult() {} ;
   RooBTagMult(const char *name, const char *title,
	       RooAbsReal &_x, 
	       RooAbsReal &_f0, RooAbsReal &_f1, RooAbsReal &_f2, RooAbsReal &_f3);
   RooBTagMult(const RooBTagMult& other,
      const char* name = 0);
   virtual TObject* clone(const char* newname) const { return new RooBTagMult(*this,newname); }
   inline virtual ~RooBTagMult() { }

   Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
   Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:

   RooRealProxy X;        // dependent variable
   RooRealProxy f0;       // nBtag=0 fraction
   RooRealProxy f1;       // nBtag=1 fraction
   RooRealProxy f2;       // nBtag=2 fraction
   RooRealProxy f3;       // nBtag>=3 fraction

   Double_t evaluate() const;
private:
  ClassDef(RooBTagMult,1) // BTagMult function
};
//---------------------------------------------------------------------------
#endif
