//____________________________________________________________________________
/*!

\class    genie::LlewellynSmithModelNC

\brief    Concrete implementation of the QELFormFactorsModelI :
          Form Factors for Quasi Elastic NC vN scattering according to
          Llewellyn-Smith model

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#ifndef _LLEWELLYN_SMITH_MODEL_NC_H_
#define _LLEWELLYN_SMITH_MODEL_NC_H_

#include "LlewellynSmith/LlewellynSmithModel.h"

namespace genie {

class LlewellynSmithModelNC : public LlewellynSmithModel {

public:

  LlewellynSmithModelNC();
  LlewellynSmithModelNC(const char * param_set);
  virtual ~LlewellynSmithModelNC();

  //-- QELFormFactorModelI interface implementation

  double F1V     (const Interaction * interaction) const;
  double xiF2V   (const Interaction * interaction) const;
  double FA      (const Interaction * interaction) const;
  double Fp      (const Interaction * interaction) const;  
};

}       // genie namespace

#endif  // _LLEWELLYN_SMITH_MODEL_NC_H_

