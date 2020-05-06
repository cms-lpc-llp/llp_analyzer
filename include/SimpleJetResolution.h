#ifndef SimpleJetResolution_h
#define SimpleJetResolution_h

#include <string>
#include <vector>
#include "SimpleJetCorrector.h"
#include <TFormula.h>


class JetCorrectorParameters;

class SimpleJetResolution 
{
 public:
  //-------- Constructors --------------
  SimpleJetResolution();
  SimpleJetResolution(const JetCorrectorParameters& fParameters);
  SimpleJetResolution(const std::string& fDataFile, const std::string& fOption);
  //-------- Destructor -----------------
  ~SimpleJetResolution();
  //-------- Member functions -----------
  float  resolution(std::vector<float>& fX,const std::vector<float>& fY) const;  
  const  JetCorrectorParameters& parameters() const {return *mParameters;} 
  SimpleJetCorrector*  getSimpleJetCorrector() const{return simpleJetCorrector_;}

 private:
  //-------- Member functions -----------
  SimpleJetResolution(const SimpleJetResolution&);
  SimpleJetResolution& operator= (const SimpleJetResolution&);
  SimpleJetCorrector*     simpleJetCorrector_;
  //-------- Member variables -----------
  JetCorrectorParameters* mParameters;
};

#endif


