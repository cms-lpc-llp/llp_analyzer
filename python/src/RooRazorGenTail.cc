//---------------------------------------------------------------------------
#include "RooFit.h"

#include "Riostream.h"
#include <TMath.h>
#include <math.h>

#include "RooRazorGenTail.h"
#include "RooRealVar.h"

using namespace std;

ClassImp(RooRazorGenTail)
//---------------------------------------------------------------------------
RooRazorGenTail::RooRazorGenTail(const char *name, const char *title,
                               RooAbsReal &_x,  RooAbsReal &_y, 
                               RooAbsReal &_bx, RooAbsReal &_by, RooAbsReal &_bxx, 
                               RooAbsReal &_bxy, RooAbsReal &_byy, RooAbsReal &_alpha) : RooAbsPdf(name, title), 
  X("X", "X Observable", this, _x),
  Y("Y", "Y Observable", this, _y),
  BX("BX", "B X Coefficient", this, _bx),
  BY("BY", "B Y Coefficient", this, _by),
  BXX("BXX", "B XX Coefficient", this, _bxx),
  BXY("BXY", "B XY Coefficient", this, _bxy),
  BYY("BYY", "B YY Coefficient", this, _byy),
  ALPHA("ALPHA", "Shape parameter", this, _alpha)
{
}
//---------------------------------------------------------------------------
RooRazorGenTail::RooRazorGenTail(const RooRazorGenTail& other, const char* name) :
   RooAbsPdf(other, name), 
   X("X", this, other.X), 
   Y("Y", this, other.Y), 
   BX("BX", this, other.BX),
   BY("BY", this, other.BY),
   BXX("BXX", this, other.BXX),
   BXY("BXY", this, other.BXY),
   BYY("BYY", this, other.BYY),
   ALPHA("ALPHA", this, other.ALPHA)
{
}
//---------------------------------------------------------------------------
Double_t RooRazorGenTail::evaluate() const
{
  double myexp = pow(BX*X/1.E+3 + BY*Y + BXX*X*X/1.E+6 + BXY*X*Y/1.E+3 + BYY*Y*Y, ALPHA);
  return exp(-myexp);
}

// //---------------------------------------------------------------------------
