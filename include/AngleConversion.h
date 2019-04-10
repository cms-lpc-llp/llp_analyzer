//----------------------------------------------------------------------------
#ifndef ANGLECONVERSION_H_11223_JASKDJGHASHDIKFAHNWEIKHFIEHWFI
#define ANGLECONVERSION_H_11223_JASKDJGHASHDIKFAHNWEIKHFIEHWFI
//----------------------------------------------------------------------------
// Conversion of the HZZ angles <==> lepton vectors
// Author: Yi Chen (11223)
//----------------------------------------------------------------------------
#include <ostream>
//----------------------------------------------------------------------------
#include "TauHelperFunctions2.h"
//----------------------------------------------------------------------------
struct EventParameters;
struct LeptonVectors;
LeptonVectors ConvertAnglesToVectors(const EventParameters &Angles, double HiggsPT, double HiggsEta);
EventParameters ConvertVectorsToAngles(const LeptonVectors &Leptons);
std::ostream &operator <<(std::ostream &out, const EventParameters &Parameters);
std::ostream &operator <<(std::ostream &out, const LeptonVectors &Leptons);
//----------------------------------------------------------------------------
struct EventParameters
{
public:
   double Phi0;
   double Theta0;
   double Phi;
   double Theta1;
   double Theta2;
   double HMass;
   double ZMass;
   double Z2Mass;
   double PhiH;
public:
   EventParameters()
      : Phi0(0), Theta0(0), Phi(0), Theta1(0), Theta2(0), HMass(125),
      ZMass(91), Z2Mass(25), PhiH(0) {}
   ~EventParameters() {}
   bool operator ==(const EventParameters &other);
};
//----------------------------------------------------------------------------
struct LeptonVectors
{
public:
   FourVector Lepton11;
   FourVector Lepton12;   // these two add up to on-shell Z
   FourVector Lepton21;
   FourVector Lepton22;   // these two add up to off-shell Z
public:
   LeptonVectors() {}
   ~LeptonVectors() {}
};
//----------------------------------------------------------------------------
LeptonVectors ConvertAnglesToVectors(const EventParameters &Angles, double HiggsPT, double HiggsEta)
{
   // Collect information
   double HMass = Angles.HMass;
   double ZMass = Angles.ZMass;
   double Z2Mass = Angles.Z2Mass;
   double PhiOffset = Angles.PhiH;
   double Phi0 = Angles.Phi0;
   double Theta0 = Angles.Theta0;
   double Phi = Angles.Phi;
   double Theta1 = Angles.Theta1;
   double Theta2 = Angles.Theta2;

   double Gamma1 = HMass / (2 * ZMass) * (1 + (ZMass * ZMass - Z2Mass * Z2Mass) / (HMass * HMass));
   double Gamma2 = HMass / (2 * Z2Mass) * (1 - (ZMass * ZMass - Z2Mass * Z2Mass) / (HMass * HMass));
   double Beta1 = GammaToBeta(Gamma1);
   double Beta2 = GammaToBeta(Gamma2);

   if(HMass < ZMass + Z2Mass)
      std::cerr << "[ConvertAnglesToVectors] WTF!  HMass < ZMass + Z2Mass.  Expect errors" << std::endl;

   // Calculate Higgs direction and boost in the ZZ frame (in that famous illustration used everywhere)
   FourVector BoostDirection;
   BoostDirection[0] = 1;
   BoostDirection[1] = -sin(Theta0) * sin(Phi0 - PI);
   BoostDirection[2] = sin(Theta0) * cos(Phi0 - PI);
   BoostDirection[3] = cos(Theta0);

   double HiggsBoostGamma =
      sqrt(HiggsPT * HiggsPT / HMass / HMass * cosh(HiggsEta) * cosh(HiggsEta) + 1);
   double HiggsBoostBeta = GammaToBeta(HiggsBoostGamma);

   // Write down the four-vectors in the ZZ frame
   FourVector Lepton11, Lepton12, Lepton21, Lepton22;

   Lepton11[0] = Gamma1 * (1 + Beta1 * cos(Theta1));
   Lepton11[1] = -sin(Phi) * sin(Theta1);
   Lepton11[2] = -cos(Phi) * sin(Theta1);
   Lepton11[3] = -Gamma1 * (Beta1 + cos(Theta1));
   Lepton11 = Lepton11 * ZMass / 2;

   Lepton12[0] = Gamma1 * (1 - Beta1 * cos(Theta1));
   Lepton12[1] = sin(Phi) * sin(Theta1);
   Lepton12[2] = cos(Phi) * sin(Theta1);
   Lepton12[3] = -Gamma1 * (Beta1 - cos(Theta1));
   Lepton12 = Lepton12 * ZMass / 2;

   Lepton21[0] = Gamma2 * (1 + Beta2 * cos(Theta2));
   Lepton21[1] = 0;
   Lepton21[2] = sin(Theta2);
   Lepton21[3] = Gamma2 * (Beta2 + cos(Theta2));
   Lepton21 = Lepton21 * Z2Mass / 2;

   Lepton22[0] = Gamma2 * (1 - Beta2 * cos(Theta2));
   Lepton22[1] = 0;
   Lepton22[2] = -sin(Theta2);
   Lepton22[3] = Gamma2 * (Beta2 - cos(Theta2));
   Lepton22 = Lepton22 * Z2Mass / 2;

   // Write down final higgs direction in the lab frame
   FourVector HiggsDirection;
   HiggsDirection.SetPtEtaPhiMass(HiggsPT, HiggsEta, 0, HMass);
   if(HiggsPT < 1e-5)
      HiggsDirection = FourVector(1, 0, 0, 1);

   // Calculate what is the rotation axis needed (and amount) to go from boosted-frame to lab frame
   FourVector RotationAxis = BoostDirection.SpatialCross(HiggsDirection).SpatialNormalize();
   double RotationAngle = BoostDirection.SpatialDot(HiggsDirection)
      / HiggsDirection.GetP() / BoostDirection.GetP();
   RotationAngle = acos(RotationAngle);

   // Boost to rotated lab frame
   Lepton11 = Lepton11.Boost(BoostDirection, HiggsBoostBeta);
   Lepton12 = Lepton12.Boost(BoostDirection, HiggsBoostBeta);
   Lepton21 = Lepton21.Boost(BoostDirection, HiggsBoostBeta);
   Lepton22 = Lepton22.Boost(BoostDirection, HiggsBoostBeta);

   // Rotate to lab frame
   Lepton11 = Lepton11.Rotate(RotationAxis, RotationAngle);
   Lepton12 = Lepton12.Rotate(RotationAxis, RotationAngle);
   Lepton21 = Lepton21.Rotate(RotationAxis, RotationAngle);
   Lepton22 = Lepton22.Rotate(RotationAxis, RotationAngle);

   // The extra degree of freedom that is not discussed anywhere - will become relevant for S != 0
   Lepton11 = Lepton11.Rotate(HiggsDirection, PhiOffset);
   Lepton12 = Lepton12.Rotate(HiggsDirection, PhiOffset);
   Lepton21 = Lepton21.Rotate(HiggsDirection, PhiOffset);
   Lepton22 = Lepton22.Rotate(HiggsDirection, PhiOffset);

   // Collect output
   LeptonVectors Result;
   Result.Lepton11 = Lepton11;
   Result.Lepton12 = Lepton12;
   Result.Lepton21 = Lepton21;
   Result.Lepton22 = Lepton22;

   return Result;
}
//----------------------------------------------------------------------------
EventParameters ConvertVectorsToAngles(const LeptonVectors &Leptons)
{
   // Some short-hand notation
   FourVector L11 = Leptons.Lepton11;
   FourVector L12 = Leptons.Lepton12;
   FourVector L21 = Leptons.Lepton21;
   FourVector L22 = Leptons.Lepton22;

   // Calculate basic information
   double HMass = (L11 + L12 + L21 + L22).GetMass();
   double ZMass = (L11 + L12).GetMass();
   double Z2Mass = (L21 + L22).GetMass();

   // Calculate higgs direction in lab frame, and boost to go to CM frame
   FourVector HiggsLab = L11 + L12 + L21 + L22;
   double HiggsBoostGamma = HiggsLab[0] / HMass;
   double HiggsBoostBeta = GammaToBeta(HiggsBoostGamma);

   // Boost everything back to the higgs rest frame
   L11 = L11.Boost(HiggsLab, -HiggsBoostBeta);
   L12 = L12.Boost(HiggsLab, -HiggsBoostBeta);
   L21 = L21.Boost(HiggsLab, -HiggsBoostBeta);
   L22 = L22.Boost(HiggsLab, -HiggsBoostBeta);

   double Theta0 = PI - GetAngle(L11 + L12, HiggsLab);

   // Z directions
   FourVector Z1 = L11 + L12;
   FourVector Z2 = L21 + L22;

   // Subtract out Z1-projection from L11, and same thing for L21 and Higgs
   // From these we can cook up the phi's
   FourVector L11Perpendicular = L11 - Z1 * (L11.SpatialDot(Z1)) / Z1.GetP2();
   FourVector L21Perpendicular = L21 - Z1 * (L21.SpatialDot(Z1)) / Z1.GetP2();
   FourVector HPerpendicular = HiggsLab - Z1 * (HiggsLab.SpatialDot(Z1)) / Z1.GetP2();

   double Phi0 = GetAngle(-HPerpendicular, L21Perpendicular);
   if(Z2.SpatialDot(HPerpendicular.SpatialCross(L21Perpendicular)) < 0)
      Phi0 = 2 * PI - Phi0;

   double Phi = GetAngle(-L21Perpendicular, L11Perpendicular);
   if(Z2.SpatialDot(L21Perpendicular.SpatialCross(L11Perpendicular)) < 0)
      Phi = 2 * PI - Phi;
   while(Phi < 0)         Phi = Phi + 2 * PI;
   while(Phi >= 2 * PI)   Phi = Phi - 2 * PI;

   // now, finally, the remaining theta's
   double Z1BoostBeta = GammaToBeta(Z1[0] / ZMass);
   double Z2BoostBeta = GammaToBeta(Z2[0] / Z2Mass);

   double Theta1 = GetAngle(Z1, L11.Boost(Z1, -Z1BoostBeta));
   double Theta2 = GetAngle(Z2, L21.Boost(Z2, -Z2BoostBeta));

   // collect output
   EventParameters Result;
   Result.Phi0 = Phi0;
   Result.Theta0 = Theta0;
   Result.Phi = Phi;
   Result.Theta1 = Theta1;
   Result.Theta2 = Theta2;
   Result.HMass = HMass;
   Result.ZMass = ZMass;
   Result.Z2Mass = Z2Mass;
   Result.PhiH = 0;
   return Result;
}
//----------------------------------------------------------------------------
std::ostream &operator <<(std::ostream &out, EventParameters &Value)
{
   out << "Event parameters: " << std::endl;
   out << "   Phi0   = " << Value.Phi0 << std::endl;
   out << "   Theta0 = " << Value.Theta0 << std::endl;
   out << "   Phi    = " << Value.Phi << std::endl;
   out << "   Theta1 = " << Value.Theta1 << std::endl;
   out << "   Theta2 = " << Value.Theta2 << std::endl;
   out << "   PhiH   = " << Value.PhiH << std::endl;
   out << "   HMass  = " << Value.HMass << std::endl;
   out << "   ZMass  = " << Value.ZMass << std::endl;
   out << "   Z2Mass = " << Value.Z2Mass << std::endl;

   return out;
}
//----------------------------------------------------------------------------
std::ostream &operator <<(std::ostream &out, LeptonVectors &Leptons)
{
   out << "Vectors of the 4 leptons = " << std::endl;
   out << "   L11 = " << Leptons.Lepton11 << std::endl;
   out << "   L12 = " << Leptons.Lepton12 << std::endl;
   out << "   L21 = " << Leptons.Lepton21 << std::endl;
   out << "   L22 = " << Leptons.Lepton22 << std::endl;

   return out;
}
//----------------------------------------------------------------------------
bool EventParameters::operator ==(const EventParameters &other)
{
   if(std::min(fabs(Phi0 - other.Phi0), 2 * PI - std::fabs(Phi0 - other.Phi0)) > 1e-5)   return false;
   if(std::fabs(Theta0 - other.Theta0) > 1e-5)                                           return false;
   if(std::min(fabs(Phi - other.Phi), 2 * PI - std::fabs(Phi - other.Phi)) > 1e-5)       return false;
   if(std::fabs(Theta1 - other.Theta1) > 1e-5)                                           return false;
   if(std::fabs(Theta2 - other.Theta2) > 1e-5)                                           return false;

   return true;
}
//----------------------------------------------------------------------------
#endif



