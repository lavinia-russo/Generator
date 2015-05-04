//____________________________________________________________________________
/*!

\class    genie::BergerSehgalRESPXSec2014

\brief    Computes the double differential cross section for resonance 
          electro- or neutrino-production according to the Rein-Sehgal model.

          The computed cross section is the d^2 xsec/ dQ^2 dW \n

          where \n
            \li \c Q^2 : momentum transfer ^ 2
            \li \c W   : invariant mass of the final state hadronic system

          Is a concrete implementation of the XSecAlgorithmI interface.

          Modifications based on a MiniBooNE tune courtesy of J. Nowak
          (http://www.physics.lancs.ac.uk/people/jaroslaw-nowak) and 
          S. Dytman.

\ref      main model: Berger, Sehgal Phys. Rev. D76, 113004 (2007) \n
          alternate: Kuzmin, Lyubushkin, Naumov Mod. Phys. Lett. A19 (2004) 2815 \n
          modifications within original format of 
	  D.Rein and L.M.Sehgal, Neutrino Excitation of Baryon Resonances
          and Single Pion Production, Ann.Phys.133, 79 (1981)

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 05, 2004

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _BERGER_SEHGAL_RES_PXSEC_2014_H_
#define _BERGER_SEHGAL_RES_PXSEC_2014_H_

#include "Base/XSecAlgorithmI.h"
#include "BaryonResonance/BaryonResonance.h"
#include "ReinSehgal/FKR.h"

namespace genie {

  class RSHelicityAmplModelI;
  class Spline;
  class XSecIntegratorI;

  class BergerSehgalRESPXSec2014 : public XSecAlgorithmI {

    public:
      BergerSehgalRESPXSec2014();
      BergerSehgalRESPXSec2014(string config);
      virtual ~BergerSehgalRESPXSec2014();

      // implement the XSecAlgorithmI interface 
      double XSec         (const Interaction * i, KinePhaseSpace_t k) const;
      double Integral     (const Interaction * i) const;
      bool   ValidProcess (const Interaction * i) const;

      // overload the Algorithm::Configure() methods to load private data
      // members from configuration options
      void Configure(const Registry & config);
      void Configure(string config);

    private:

      void LoadConfig (void);

      mutable FKR fFKR;

      const RSHelicityAmplModelI * fHAmplModelCC;
      const RSHelicityAmplModelI * fHAmplModelNCp;
      const RSHelicityAmplModelI * fHAmplModelNCn;
      const RSHelicityAmplModelI * fHAmplModelEMp;
      const RSHelicityAmplModelI * fHAmplModelEMn;

      // configuration data
      bool     fWghtBW;            ///< weight with resonance breit-wigner?
      double   fZeta;              ///< FKR parameter Zeta
      double   fOmega;             ///< FKR parameter Omega
      double   fMa2;               ///< (axial mass)^2
      double   fMv2;               ///< (vector mass)^2
      double   fSin48w;            ///< sin^4(Weingberg angle)
      bool     fUsingDisResJoin;   ///< use a DIS/RES joining scheme?
      bool     fUsingNuTauScaling; ///< use NeuGEN nutau xsec reduction factors?
      double   fWcut;              ///< apply DIS/RES joining scheme < Wcut
      double   fN2ResMaxNWidths;   ///< limits allowed phase space for n=2 res
      double   fN0ResMaxNWidths;   ///< limits allowed phase space for n=0 res
      double   fGnResMaxNWidths;   ///< limits allowed phase space for other res
      Spline * fNuTauRdSpl;        ///< xsec reduction spline for nu_tau
      Spline * fNuTauBarRdSpl;     ///< xsec reduction spline for nu_tau_bar


      bool fKNL;
      bool fBRS;
      bool fGA;
      bool fGV;

      const XSecIntegratorI * fXSecIntegrator;
  };

}       // genie namespace

#endif  // _BERGER_SEHGAL_RES_PXSEC_2014_H_