//---------------------------------------------------------------------------
#ifndef ROO_RazorGenTail
#define ROO_RazorGenTail
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
class RooRazorGenTail : public RooAbsPdf
{
public:
   RooRazorGenTail() {} ;
   RooRazorGenTail(const char *name, const char *title,
                  RooAbsReal &_x, RooAbsReal &_y, 
                  RooAbsReal &_bx, RooAbsReal &_by, RooAbsReal &_bxx, 
                  RooAbsReal &_bxy, RooAbsReal &_byy, RooAbsReal &_alpha);
   RooRazorGenTail(const RooRazorGenTail& other,
      const char* name = 0);
   virtual TObject* clone(const char* newname) const { return new RooRazorGenTail(*this,newname); }
   inline virtual ~RooRazorGenTail() { }

protected:

   RooRealProxy X;        // dependent variable
   RooRealProxy Y;        // dependent variable
   RooRealProxy BX;       // X offset
   RooRealProxy BY;       // Y offset
   RooRealProxy BXX;        // shape parameter
   RooRealProxy BXY;        // shape parameter
   RooRealProxy BYY;        // shape parameter
   RooRealProxy ALPHA;        // shape parameter

   Double_t evaluate() const;
private:
  ClassDef(RooRazorGenTail,1) // Razor2DTail_SYS function
    
};
//---------------------------------------------------------------------------
#endif
