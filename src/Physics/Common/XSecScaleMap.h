//____________________________________________________________________________
/*!

\class    genie::XSecScaleMap

\brief    This class is responsible to compute a scaling factor for the XSec

\author   Code contributed by J.Tena Vidal and M.Roda

\created  June, 2020

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _XSEC_SCALE_MAP_H_
#define _XSEC_SCALE_MAP_H_

#include "Physics/Common/XSecScaleI.h"
#include "Framework/Interaction/Interaction.h"
#include <map>
#include <TSpline.h>

using std::map; 

namespace genie {
  
  class XSecScaleMap: public XSecScaleI {
    
  public:
    XSecScaleMap();
    XSecScaleMap(string config);
    virtual ~XSecScaleMap();
    
    // This function returns the scaling value for a given interaction:
    virtual double GetScaling( const Interaction & ) const ; 
    
    void Configure (const Registry & config);
    void Configure (string config);
    
  protected:
    
    // Load algorithm configuration
    void LoadConfig (void);
    
 private: 
    const XSecScaleI * fXSecScaleDefault ; 
    std::map<double,XSecScaleI *> fXSecScaleMap ;

  };
  
}       // genie namespace
#endif  // _XSEC_SCALE_MAP_H_
