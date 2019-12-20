//---------------------------------------------------------------------------
#include "RooFit.h"

#include "Riostream.h"
#include <TMath.h>
#include <math.h>

#include "RooRazor2DTail_SYS.h"
#include "RooRealVar.h"

using namespace std;

ClassImp(RooRazor2DTail_SYS)
//---------------------------------------------------------------------------
RooRazor2DTail_SYS::RooRazor2DTail_SYS(const char *name, const char *title,
                               RooAbsReal &_x,  RooAbsReal &_y, 
                               RooAbsReal &_x0, RooAbsReal &_y0, 
                               RooAbsReal &_b, RooAbsReal &_n) : RooAbsPdf(name, title), 
  X("X", "X Observable", this, _x),
  Y("Y", "Y Observable", this, _y),
  X0("X0", "X Offset", this, _x0),
  Y0("Y0", "Y Offset", this, _y0),
  B("B", "Shape parameter", this, _b),
  N("N", "Shape parameter", this, _n)
{
}
//---------------------------------------------------------------------------
RooRazor2DTail_SYS::RooRazor2DTail_SYS(const RooRazor2DTail_SYS& other, const char* name) :
   RooAbsPdf(other, name), 
   X("X", this, other.X), 
   Y("Y", this, other.Y), 
   X0("X0", this, other.X0),
   Y0("Y0", this, other.Y0),
   B("B", this, other.B),
   N("N", this, other.N)
{
}
//---------------------------------------------------------------------------
Double_t RooRazor2DTail_SYS::evaluate() const
{
  double myexp = B*N*pow(fabs(X-X0)*fabs(Y-Y0),1./N);
  double mycoeff = B*pow(fabs(X-X0)*fabs(Y-Y0),1./N) - 1.;
  if(myexp < -700.) {
    //std::cout << "MYEXP = "<< myexp << ", < -700 -- BAD" << std::endl;
    return  1.7e-308;}
  if(mycoeff <= 0.) {
    //std::cout << "MYCOEFF = " << mycoeff <<", IS NEGATIVE -- BAD" << std::endl;
    return  1.7e-308;}
  if(X0 >= X.min() || Y0 >= Y.min() || B <= 0. || N <= 0.) {
    //std::cout << "PARAMETERS OUT OF PHYSICAL, INTEGRABLE RANGES -- BAD" << std::endl;
    return  1.7e-308;}
  return mycoeff*exp(-myexp);
}

// //---------------------------------------------------------------------------
Int_t RooRazor2DTail_SYS::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const{
  if(rangeName != 0) cout << "Warning: rangeName not used in calculation of integral" << endl;
  // integral on both X and Y
  if (matchArgs(allVars, analVars, X, Y)) return 1;
  // integral over X
  if (matchArgs(allVars, analVars, X)) return 2;
  // integral over Y
  if (matchArgs(allVars, analVars, Y)) return 3;
  // integrating nothing
  return 0;
}

// //---------------------------------------------------------------------------
Double_t RooRazor2DTail_SYS::analyticalIntegral(Int_t code, const char* rangeName) const{

   const Double_t xmin = X.min(rangeName); const Double_t xmax = X.max(rangeName);
   const Double_t ymin = Y.min(rangeName); const Double_t ymax = Y.max(rangeName);

   if(B <= 0. || N <= 0. || X0 >= X.min() || Y0 >= Y.min()) return 1.;

   double integral = 0.;
   if(code ==1) { // integral on both X and Y
     integral = N/pow(B*N,N)*(Gfun(xmin,ymin)-Gfun(xmin,ymax)-Gfun(xmax,ymin)+Gfun(xmax,ymax));
   } else if(code == 2) { // integral on X
     integral = ( (xmin-X0)*exp(B*N*pow(xmax-X0,1/N)*pow(Y-Y0,1/N)) - (xmax-X0)*exp(B*N*pow(xmin-X0,1/N)*pow(Y-Y0,1/N)) )*exp(-B*N*(pow(xmin-X0,1/N)+pow(xmax-X0,1/N))*pow(Y-Y0,1/N));
   } else if(code == 3) { // integral on Y
     integral = ( (ymin-Y0)*exp(B*N*pow(ymax-Y0,1/N)*pow(X-X0,1/N)) - (ymax-Y0)*exp(B*N*pow(ymin-Y0,1/N)*pow(X-X0,1/N)) )*exp(-B*N*(pow(ymin-Y0,1/N)+pow(ymax-Y0,1/N))*pow(X-X0,1/N));
   } else {
     cout << "WARNING IN RooRazor2DTail_SYS: integration code is not correct" << endl;
     cout << "                           what are you integrating on?" << endl;
     return 1.;
   }

   return integral;
}
// //---------------------------------------------------------------------------
