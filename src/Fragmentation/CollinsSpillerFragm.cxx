//____________________________________________________________________________
/*!

\class    genie::CollinsSpillerFragm

\brief    The Collins-Spiller fragmentation function.

          Is a concrete implementation of the FragmentationFunctionI interface.
          
\ref      P.D.B.Collins and T.P.Spiller, J.Phys.G11, 1289 (1984)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 15, 2004

*/ 
//____________________________________________________________________________

#include "Fragmentation/CollinsSpillerFragm.h"
#include "Fragmentation/fragmentation_functions.h"

using namespace genie;

//___________________________________________________________________________
CollinsSpillerFragm::CollinsSpillerFragm() :
FragmentationFunctionI()
{
  fName     = "genie::CollinsSpillerFragm";
  fParamSet = "Default";

  FindConfig();
}
//___________________________________________________________________________
CollinsSpillerFragm::CollinsSpillerFragm(const char * param_set) :
FragmentationFunctionI(param_set)
{
  fName = "genie::CollinsSpillerFragm";

  FindConfig();
}
//___________________________________________________________________________
CollinsSpillerFragm::~CollinsSpillerFragm()
{
  delete fFunc;
}
//___________________________________________________________________________
double CollinsSpillerFragm::Value(double z) const
{
  //-- get fragmentation function parameters from the config. registry

  double N = (fConfig->Exists("norm"))    ? fConfig->GetDouble("norm")    : 0;
  double e = (fConfig->Exists("epsilon")) ? fConfig->GetDouble("epsilon") : 0;

  //-- evaluate the fragmentation function

  assert( z > 0 && z < 1);

  fFunc->SetParameters(N,e);

  double D = fFunc->Eval(z);

  return D;
}
//___________________________________________________________________________
double CollinsSpillerFragm::GenerateZ(void) const
{
  //-- get fragmentation function parameters from the config. registry

  double N = (fConfig->Exists("norm"))    ? fConfig->GetDouble("norm")    : 0;
  double e = (fConfig->Exists("epsilon")) ? fConfig->GetDouble("epsilon") : 0;

  fFunc->SetParameters(N,e);

  double Rndm = fFunc->GetRandom();

  return Rndm;
}
//___________________________________________________________________________
void CollinsSpillerFragm::BuildFunction(void)
{
  fFunc = new TF1("fFunc",genie::collins_spiller_fragmentation_function,0,1,2);

  fFunc->SetParNames("Norm","Epsilon");
}
//___________________________________________________________________________



