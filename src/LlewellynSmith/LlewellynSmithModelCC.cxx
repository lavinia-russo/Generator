//____________________________________________________________________________
/*!

\class    genie::LlewellynSmithModelCC

\brief    Form Factors for Quasi Elastic CC vN scattering according to
          Llewellyn-Smith model.

          Is a concrete implementation of the QELFormFactorsModelI :

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory
          
\created  May 03, 2004

*/ 
//____________________________________________________________________________

#include "Conventions/Constants.h"
#include "LlewellynSmith/LlewellynSmithModelCC.h"
#include "Messenger/Messenger.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
LlewellynSmithModelCC::LlewellynSmithModelCC() :
LlewellynSmithModel()
{
  fName = "genie::LlewellynSmithModelCC";
}
//____________________________________________________________________________
LlewellynSmithModelCC::LlewellynSmithModelCC(const char * param_set) :
LlewellynSmithModel(param_set)
{
  fName = "genie::LlewellynSmithModelCC";

  FindConfig();
}
//____________________________________________________________________________
LlewellynSmithModelCC::~LlewellynSmithModelCC()
{

}
//____________________________________________________________________________
double LlewellynSmithModelCC::F1V(const Interaction * interaction) const
{
  return LlewellynSmithModel::F1V(interaction);
}  
//____________________________________________________________________________
double LlewellynSmithModelCC::xiF2V(const Interaction * interaction) const
{
  return LlewellynSmithModel::xiF2V(interaction);
}
//____________________________________________________________________________
double LlewellynSmithModelCC::FA(const Interaction * interaction) const
{
  return LlewellynSmithModel::FA(interaction);
}
//____________________________________________________________________________
double LlewellynSmithModelCC::Fp(const Interaction * interaction) const
{
  return LlewellynSmithModel::Fp(interaction);
}
//____________________________________________________________________________



