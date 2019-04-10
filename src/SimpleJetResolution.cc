#include "SimpleJetResolution.h"
#include "JetCorrectorParameters.h"
#include "JetCorrUtilities.cc"
#include <iostream>
#include <sstream>
#include <cmath>

//------------------------------------------------------------------------ 
//--- Default SimpleJetResolution constructor -----------------------------
//------------------------------------------------------------------------
SimpleJetResolution::SimpleJetResolution() 
{ 
  mParameters      = new JetCorrectorParameters();
  simpleJetCorrector_  = new SimpleJetCorrector();
}
//------------------------------------------------------------------------ 
//--- SimpleJetResolution constructor -------------------------------------
//--- reads arguments from a file ----------------------------------------
//------------------------------------------------------------------------
SimpleJetResolution::SimpleJetResolution(const std::string& fDataFile, const std::string& fOption) 
{
  mParameters      = new JetCorrectorParameters(fDataFile,fOption);
  simpleJetCorrector_= new SimpleJetCorrector(*mParameters);
}
//------------------------------------------------------------------------
//--- SimpleJetResolution constructor -------------------------------------
//--- reads arguments from a file ----------------------------------------
//------------------------------------------------------------------------
SimpleJetResolution::SimpleJetResolution(const JetCorrectorParameters& fParameters)
{
  mParameters      = new JetCorrectorParameters(fParameters);
  simpleJetCorrector_= new SimpleJetCorrector(*mParameters);
}
//------------------------------------------------------------------------ 
//--- SimpleJetResolution destructor --------------------------------------
//------------------------------------------------------------------------
SimpleJetResolution::~SimpleJetResolution() 
{
  delete simpleJetCorrector_;
}
//------------------------------------------------------------------------ 
//--- calculates the correction ------------------------------------------
//------------------------------------------------------------------------
float SimpleJetResolution::resolution(std::vector<float>& fX,const std::vector<float>& fY) const 
{
   fX[0]=std::abs(fX[0]); // use only positive eta values 
   return simpleJetCorrector_->correction(fX,fY);

}









