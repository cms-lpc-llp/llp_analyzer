//---------------------------------------------------------------------------
#include "RooFit.h"

#include "Riostream.h"
#include <TMath.h>
#include <math.h>

#include "RooBTagMult.h"
#include "RooRealVar.h"

using namespace std;

ClassImp(RooBTagMult)
//---------------------------------------------------------------------------
RooBTagMult::RooBTagMult(const char *name, const char *title,
			       RooAbsReal &_x,
			       RooAbsReal &_f0,
			       RooAbsReal &_f1,
			       RooAbsReal &_f2,
			       RooAbsReal &_f3) : RooAbsPdf(name, title), 
  X("X", "X Observable", this, _x),
  f0("f0", "nBtag=0 fraction", this, _f0),
  f1("f1", "nBtag=1 fraction", this, _f1),
  f2("f2", "nBtag=2 fraction", this, _f2),
  f3("f3", "nBtag=3 fraction", this, _f3)
{
}
//---------------------------------------------------------------------------
RooBTagMult::RooBTagMult(const RooBTagMult& other, const char* name) :
   RooAbsPdf(other, name), 
   X("X", this, other.X), 
   f0("f0", this, other.f0),
   f1("f1", this, other.f1),
   f2("f2", this, other.f2),
   f3("f3", this, other.f3)
{
}
//---------------------------------------------------------------------------
Double_t RooBTagMult::evaluate() const
{
  double thisf0 = f0;
  double thisf1 = f1;
  double thisf2 = f2;
  double thisf3 = f3;

  if (X<0.) return 1.7e-308;
  else if(X>=0. && X<1.) return thisf0;
  else if(X>=1. && X<2.) return thisf1;
  else if(X>=2. && X<3.) return thisf2;
  else if(X>=3. && X<4.) return thisf3;
  else return  1.7e-308;
}

// //---------------------------------------------------------------------------
Int_t RooBTagMult::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const{
  // integral over X
  if (matchArgs(allVars, analVars, X)) return 1;
  return 0;
}

// //---------------------------------------------------------------------------
Double_t RooBTagMult::analyticalIntegral(Int_t code, const char* rangeName) const{
  const Double_t xmin = X.min(rangeName); const Double_t xmax = X.max(rangeName);

  double thisf0 = f0;
  double thisf1 = f1;
  double thisf2 = f2;
  double thisf3 = f3;

  double binVol0 = 0;
  double binVol1 = 0;
  double binVol2 = 0;
  double binVol3 = 0;

  double integral = 0;
 
  if(code == 1) {
    if (xmin==0. and xmax==4.) 
      binVol0 = 1.; binVol1 = 1.; binVol2 = 1.; binVol3 = 1.;
    if (xmin>=0. && xmin <=1.){
      if (xmax>0. && xmax <=1.){
	binVol0 = xmax-xmin; binVol1 = 0; binVol2 = 0.; binVol3 = 0.;
	  }
      else if (xmax>1. && xmax<=2.){
	binVol0 = 1.-xmin; binVol1 = xmax-1.; binVol2 = 0.; binVol3 = 0; 
	  }
      else if (xmax>2. && xmax<=3.){
	binVol0 = 1.-xmin; binVol1 = 1.; binVol2 = xmax-2.; binVol3 = 0; 
	  }
      else if (xmax>3. && xmax<=4.){
	binVol0 = 1.-xmin; binVol1 = 1.; binVol2 = 1.; binVol3 = xmax-3.;
	  }
    }
    else if (xmin>=1. && xmin <=2.){
      if (xmax>1. && xmax <=2.){
	binVol0 = 0.; binVol1 = xmax-xmin; binVol2 = 0.; binVol3 = 0.;
	  }
      else if (xmax>2. && xmax<=3.){
	binVol0 = 0.; binVol1 = 2.-xmin; binVol2 = xmax-2.; binVol3 = 0.; 
	  }
      else if (xmax>3. && xmax<=4.){
	binVol0 = 0.; binVol1 = 2.-xmin; binVol2 = 1.; binVol3 = xmax-3.;
	  }
    }

    else if (xmin>=2. && xmin <=3.){
      if (xmax>2. && xmax <=3.){
	binVol0 = 0.; binVol1 = 0.; binVol2 = xmax-xmin; binVol3 = 0.;
	  }
      else if (xmax>3. && xmax<=4.){
	binVol0 = 0.; binVol1 = 0.; binVol2 = 3.-xmin; binVol3 = xmax-3.;
	  }
    }

    else if (xmin>=3. && xmin <=4.){
      if (xmax>3. && xmax<=4.){
	binVol0 = 0.; binVol1 = 0.; binVol2 = 0.; binVol3 = xmax-xmin;
	  }
    }

    integral = thisf0*binVol0 + thisf1*binVol1 + thisf2*binVol2 + thisf3*binVol3;
    
    return integral;
  }
  else {
     cout << "WARNING IN RooBTagMult: integration code is not correct" << endl;
     cout << "                           what are you integrating on?" << endl;
     return 1.;
  }
  return 1.;
}
// //---------------------------------------------------------------------------

