//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Anselmo Meregaglia <anselmo.meregaglia@cern.ch>, ETH Zurich
         Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>, STFC - Rutherford Lab
         May 24, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cassert>

#include <TGeoVolume.h>
#include <TGeoManager.h>
#include <TGeoShape.h>
#include <TGeoMedium.h>
#include <TGeoMaterial.h>
#include <TObjArray.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TSystem.h>
#include <TMath.h>

#include "Conventions/GBuild.h"
#include "Conventions/Units.h"
#include "EVGDrivers/PathLengthList.h"
#include "EVGDrivers/GFluxI.h"
#include "Geo/ROOTGeomAnalyzer.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGLibrary.h"
#include "Utils/PrintUtils.h"

#include <TPolyMarker3D.h>
#include <TGeoBBox.h>
#include "Numerical/RandomGen.h"

using namespace genie;
using namespace genie::geometry;

//___________________________________________________________________________
ROOTGeomAnalyzer::ROOTGeomAnalyzer(string geometry_filename) :
GeomAnalyzerI()
{
  this->Initialize();
  this->Load(geometry_filename);
}
//___________________________________________________________________________
ROOTGeomAnalyzer::ROOTGeomAnalyzer(TGeoManager * gm) :
GeomAnalyzerI()
{
  this->Initialize();
  this->Load(gm);
}
//___________________________________________________________________________
ROOTGeomAnalyzer::~ROOTGeomAnalyzer()
{
  this->CleanUp();
}
//___________________________________________________________________________
void ROOTGeomAnalyzer::SetLengthUnits(double u)
{
// Use the units of your input geometry, eg
//               geom.SetUnits(genie::units::centimeter)
// GENIE uses the physical system of units (hbar=c=1) almost throughtout so
// everything is expressed in GeV but when analyzing detector geometries we
// use meters. Setting your input geometry units will allow us to figure the
// conversion factor.
// As input, use one of the constants in $GENIE/src/Conventions/Units.h

  fLengthScale = u/units::meter;
  LOG("GROOTGeom", pNOTICE)
     << "Geometry length units scale factor (GU -> m): " << fLengthScale;
}
//___________________________________________________________________________
void ROOTGeomAnalyzer::SetDensityUnits(double u)
{
// Like above, but for density (def units = kgr/m3)

  fDensityScale = u / (units::kilogram / units::meter3);
  LOG("GROOTGeom", pNOTICE)
    << "Geometry density units scale factor (GU -> kgr/m3): " << fDensityScale;
}
//___________________________________________________________________________
void ROOTGeomAnalyzer::SetMaxPlSafetyFactor(double sf)
{
// Set a factor that can multiply the computed max path lengths.
// The maximum path lengths are computed by performing an MC scanning of the
// input geometry. If you configure the scanner with a low number of points
// or rays you might understimate the path lengths, so you might want to
// 'inflate' them a little bit using this method.
// Do not set this number too high, because the max interaction probability
// will be grossly overestimated and you would need lots of attempts before
// getting a flux neutrino to interact...

  fMaxPlSafetyFactor = sf;

  LOG("GROOTGeom", pNOTICE)
            << "Max path length safety factor: " << fMaxPlSafetyFactor;
}
//___________________________________________________________________________
void ROOTGeomAnalyzer::SetMixtureWeightsSum(double sum)
{
// Set it to x, if the relative weight proportions of elements in a mixture
// add up to x (eg x=1, 100, etc). Set it to a negative value if you want
// driver to be explicitly computing the correct weight normalization.

  fMixtWghtSum = sum;
}
//___________________________________________________________________________
void ROOTGeomAnalyzer::SetTopVolName(string name)
{
// Set the name of the top volume.
// This driver would ask the TGeoManager::GetTopVolume() for the top volume.
// Use this method for changing this if for example you want to set a smaller
// volume as the top one so as to generate events only in a specific part of
// your detector.

  fTopVolumeName = name;

  LOG("GROOTGeom",pNOTICE) << "Geometry Top Volume name: " << fTopVolumeName;

  TGeoVolume * gvol = fGeometry->GetVolume(fTopVolumeName.c_str());

  if(!gvol) {
     LOG("GROOTGeom",pWARN) << "Could not find volume: " << name.c_str();
     LOG("GROOTGeom",pWARN) << "Will not change the current top volume";
     fTopVolumeName = "";
     return;
  }
  fTopVolume = gvol;
}
//___________________________________________________________________________
const PathLengthList & ROOTGeomAnalyzer::ComputeMaxPathLengths(void)
{
// Computes the maximum path lengths for all materials in the input geometry
//
  LOG("GROOTGeom", pNOTICE)
                  << "Computing the maximum path lengths for all materials";

  if(!fGeometry) {
      LOG("GROOTGeom",pFATAL) << "No ROOT geometry is loaded!!";
      exit(1);
  }

  //-- initialize max path lengths
  fCurrMaxPathLengthList->SetAllToZero();

  //-- select maximum path length calculation method
  if(fFlux) this->MaxPathLengthsFluxMethod();
  else      this->MaxPathLengthsBoxMethod();

  return *fCurrMaxPathLengthList;
}
//___________________________________________________________________________
void ROOTGeomAnalyzer::MaxPathLengthsFluxMethod(void)
{
// Use the input flux driver to generate "rays", and then follow them through
// the detector and figure out the maximum path length for each material

  LOG("GROOTGeom", pNOTICE)
               << "Computing the maximum path lengths using the FLUX method";

  int iparticle = 0;
  PathLengthList::const_iterator pl_iter;

  while(iparticle < this->ScannerNParticles()) {

    bool ok = fFlux->GenerateNext();
    if(!ok) {
       LOG("GROOTGeom", pWARN) << "Couldn't generate a flux neutrino";
       continue;
    }

    const TLorentzVector & nup4  = fFlux -> Momentum ();
    const TLorentzVector & nux4  = fFlux -> Position ();

    LOG("GMCJDriver", pNOTICE)
       << "\n [-] Generated flux neutrino: "
       << "\n  |----o 4-momentum : " << utils::print::P4AsString(&nup4)
       << "\n  |----o 4-position : " << utils::print::X4AsString(&nux4);

    const PathLengthList & pl = this->ComputePathLengths(nux4, nup4);

    bool enters = false;

    for(pl_iter = pl.begin(); pl_iter != pl.end(); ++pl_iter) {
       int    pdgc = pl_iter->first;
       double pl   = pl_iter->second;

       if(pl>0) {
          pl *= (this->MaxPlSafetyFactor());

          pl = TMath::Max(pl, fCurrMaxPathLengthList->PathLength(pdgc));
          fCurrMaxPathLengthList->SetPathLength(pdgc,pl);
       }
    }
    if(enters) iparticle++;
  }
}
//___________________________________________________________________________
void ROOTGeomAnalyzer::MaxPathLengthsBoxMethod(void)
{
// Generate points in the geometry's bounding box and for each point generate
// random rays, follow them through the detector and figure out the maximum
// path length for each material

  LOG("GROOTGeom", pNOTICE)
               << "Computing the maximum path lengths using the BOX method";

  // get geometry's bounding box
  LOG("GROOTGeom", pNOTICE) << "Getting a TGeoBBox enclosing the detector";
  TGeoShape * TS  = fTopVolume->GetShape();
  TGeoBBox *  box = (TGeoBBox *)TS;

  //get box origin and dimensions (in the same units as the geometry)
  double dx = box->GetDX(); // half-length
  double dy = box->GetDY(); // half-length
  double dz = box->GetDZ(); // half-length
  double ox = (box->GetOrigin())[0];
  double oy = (box->GetOrigin())[1];
  double oz = (box->GetOrigin())[2];

  LOG("GROOTGeom",pNOTICE) << "Box size (GU)   :"
               << " x = " << 2*dx << ", y = " << 2*dy << ", z = " << 2*dz;
  LOG("GROOTGeom",pNOTICE) << "Box origin (GU) :"
               << " x = " << ox   << ", y = " << oy   << ", z = " <<   oz;

  // generate 200 random points on each surface, use 200 rays to
  // calculate maximum path for each material

  RandomGen * rnd = RandomGen::Instance();

  LOG("GROOTGeom",pNOTICE)
        << "Will generate [" << fNPoints << "] random points / box surface";
  LOG("GROOTGeom",pNOTICE)
        << "Will generate [" << fNRays   << "] rays / point";

  //loop on materials

  vector<int>::iterator itr;
  for(itr=fCurrPDGCodeList->begin();itr!=fCurrPDGCodeList->end();itr++) {

    int pdgc = *itr;
    LOG("GROOTGeom", pNOTICE)
                 <<"Calculating max path length for material: " << pdgc;

    int    ipoint    (0);
    int    iray      (0);
    int    maxPoints (fNPoints);
    int    maxRays   (fNRays);
    double maxPath   (0);

    TVector3 pos(0.,0.,0.); // position (GU)
    TVector3 dir(0.,0.,0.); // direction

    //top:
    LOG("GROOTGeom",pNOTICE) << "Box surface scanned: [TOP]";
    ipoint=0;
    while (ipoint++ < maxPoints) {
      iray=0;
      pos.SetXYZ(ox-dx+2*dx*rnd->RndGeom().Rndm(), oy+dy, oz-dz+2*dz*rnd->RndGeom().Rndm());
      while (iray++ < maxRays) {
        dir.SetXYZ(-0.5+rnd->RndGeom().Rndm(), -rnd->RndGeom().Rndm(), -0.5+rnd->RndGeom().Rndm());
        maxPath = TMath::Max(maxPath,
                   this->ComputePathLengthPDG(pos,dir.Unit(),pdgc));
      }
    }
/*
    //bottom:
    LOG("GROOTGeom",pNOTICE) << "Box surface scanned: [BOTTOM]";
    ipoint=0;
    while (ipoint++ < maxPoints) {
      iray=0;
      pos.SetXYZ(ox-dx+2*dx*rnd->RndGeom().Rndm(), oy-dy, oz-dz+2*dz*rnd->RndGeom().Rndm());
      while (iray++ < maxRays) {
        dir.SetXYZ(-0.5+rnd->RndGeom().Rndm(), rnd->RndGeom().Rndm(), -0.5+rnd->RndGeom().Rndm());
        maxPath = TMath::Max(maxPath,
                   this->ComputePathLengthPDG(pos,dir.Unit(),pdgc));
      }
    }
*/
    //left:
    LOG("GROOTGeom",pNOTICE) << "Box surface scanned: [LEFT]";
    ipoint=0;
    while (ipoint++ < maxPoints) {
      iray=0;
      pos.SetXYZ(ox-dx, oy-dy+2*dy*rnd->RndGeom().Rndm(), oz-dz+2*dz*rnd->RndGeom().Rndm());
      while (iray++ < maxRays) {
        dir.SetXYZ(rnd->RndGeom().Rndm(), -0.5+rnd->RndGeom().Rndm(), -0.5+rnd->RndGeom().Rndm());
        maxPath = TMath::Max(maxPath,
                   this->ComputePathLengthPDG(pos,dir.Unit(),pdgc));
      }
    }
/*
    //right:
    LOG("GROOTGeom",pNOTICE) << "Box surface scanned: [RIGHT]";
    ipoint=0;
    while (ipoint++ < maxPoints) {
      iray=0;
      pos.SetXYZ(ox+dx, oy-dy+2*dy*rnd->RndGeom().Rndm(), oz-dz+2*dz*rnd->RndGeom().Rndm());
      while (iray++ < maxRays) {
        dir.SetXYZ(-rnd->RndGeom().Rndm(), -0.5+rnd->RndGeom().Rndm(), -0.5+rnd->RndGeom().Rndm());
        maxPath = TMath::Max(maxPath,
                   this->ComputePathLengthPDG(pos,dir.Unit(),pdgc));
      }
    }
*/
    //front:
    LOG("GROOTGeom",pNOTICE) << "Box surface scanned: [FRONT]";
    ipoint=0;
    while (ipoint++ < maxPoints) {
      iray=0;
      pos.SetXYZ(ox-dx+2*dx*rnd->RndGeom().Rndm(), oy-dy+2*dy*rnd->RndGeom().Rndm(), oz+dz);
      while (iray++ < maxRays) {
        dir.SetXYZ(-0.5+rnd->RndGeom().Rndm(), -0.5+rnd->RndGeom().Rndm(), -rnd->RndGeom().Rndm());
        maxPath = TMath::Max(maxPath,
                   this->ComputePathLengthPDG(pos,dir.Unit(),pdgc));
      }
    }
/*
    //back:
    LOG("GROOTGeom",pNOTICE) << "Box surface scanned: [BACK]";
    ipoint=0;
    while (ipoint++ < maxPoints) {
      iray=0;
      pos.SetXYZ(ox-dx+2*dx*rnd->RndGeom().Rndm(), oy-dy+2*dy*rnd->RndGeom().Rndm(), oz-dz);
      while (iray++ < maxRays) {
        dir.SetXYZ(-0.5+rnd->RndGeom().Rndm(), -0.5+rnd->RndGeom().Rndm(), rnd->RndGeom().Rndm());
        maxPath = TMath::Max(maxPath,
                   this->ComputePathLengthPDG(pos,dir.Unit(),pdgc));
      }
    }
*/
    maxPath *= (this->MaxPlSafetyFactor());
    fCurrMaxPathLengthList->AddPathLength(pdgc, maxPath); // GU

    LOG("GROOTGeom", pNOTICE)
                       << "Max path length found = " << maxPath << " GU";
  }
  this->ScalePathLengths(*fCurrMaxPathLengthList); // GU -> m
}
//________________________________________________________________________
void ROOTGeomAnalyzer::Initialize(void)
{
  LOG("GROOTGeom", pNOTICE)
                << "Initializing ROOT geometry driver & setting defaults";

  fCurrMaxPathLengthList = 0;
  fCurrPathLengthList    = 0;
  fCurrPDGCodeList       = 0;
  fTopVolume             = 0;
  fTopVolumeName         = "";

  // some defaults:
  this -> SetScannerNPoints    (200);
  this -> SetScannerNRays      (200);
  this -> SetScannerNParticles (10000);
  this -> SetScannerFlux       (0);
  this -> SetMaxPlSafetyFactor (1.1);
  this -> SetLengthUnits       (genie::units::meter);
  this -> SetDensityUnits      (genie::units::kilogram/genie::units::meter3);
  this -> SetWeightWithDensity (true);
  this -> SetMixtureWeightsSum (-1.);
}
//___________________________________________________________________________
void ROOTGeomAnalyzer::CleanUp(void)
{
  LOG("GROOTGeom", pNOTICE) << "Cleaning up...";

  if( fCurrPathLengthList    ) delete fCurrPathLengthList;
  if( fCurrMaxPathLengthList ) delete fCurrMaxPathLengthList;
  if( fCurrPDGCodeList       ) delete fCurrPDGCodeList;
}
//___________________________________________________________________________
void ROOTGeomAnalyzer::Load(string filename)
{
  LOG("GROOTGeom", pNOTICE) << "Loading geometry from: " << filename;

  bool is_accessible = ! (gSystem->AccessPathName( filename.c_str() ));
  if (!is_accessible) {
     LOG("GROOTGeom", pERROR)
       << "The ROOT geometry doesn't exist! Initialization failed!";
     exit(1);
  }
  TGeoManager * gm = TGeoManager::Import(filename.c_str());

  this->Load(gm);
}
//___________________________________________________________________________
void ROOTGeomAnalyzer::Load(TGeoManager * gm)
{
  LOG("GROOTGeom", pNOTICE)
         << "A TGeoManager is being loaded to the geometry driver";
  fGeometry = gm;

  if(!fGeometry) {
    LOG("GROOTGeom", pFATAL) << "Null TGeoManager! Aborting";
  }
  assert(fGeometry);

  this->BuildListOfTargetNuclei();

  const PDGCodeList & pdglist = this->ListOfTargetNuclei();

  fTopVolume             = 0;
  fCurrPathLengthList    = new PathLengthList(pdglist);
  fCurrMaxPathLengthList = new PathLengthList(pdglist);
  fCurrVertex            = new TVector3(0.,0.,0.);

  // ask geometry manager for its top volume
  fTopVolume = fGeometry->GetTopVolume();
  if(!fTopVolume) {
      LOG("GROOTGeom", pFATAL) << "Could not get top volume!!!";
  }
  assert(fTopVolume);
}
//___________________________________________________________________________
const PDGCodeList & ROOTGeomAnalyzer::ListOfTargetNuclei(void)
{
  return *fCurrPDGCodeList;
}
//___________________________________________________________________________
const PathLengthList & ROOTGeomAnalyzer::ComputePathLengths(
                          const TLorentzVector & x, const TLorentzVector & p)
{
// Computes the path-length within each detector material for a neutrino
// starting from point x and travelling along the direction of p.

  LOG("GROOTGeom", pNOTICE)
       << "Computing path-lengths for the input neutrino";

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GROOTGeom", pDEBUG)
       << "\nInput nu: 4p (GeV) = " << utils::print::P4AsShortString(&p)
       << ", 4x (m,s) = " << utils::print::X4AsString(&x);
#endif

  // reset current list of path-lengths
  fCurrPathLengthList->SetAllToZero();

  //loop over materials & compute the path-length
  vector<int>::iterator itr;
  for(itr=fCurrPDGCodeList->begin();itr!=fCurrPDGCodeList->end();itr++) {

    int pdgc = *itr;
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("GROOTGeom", pINFO)
           <<"Calculating path length for material: " << pdgc;
#endif

    TVector3 pos  = x.Vect();        // initial position
    pos *= (1./this->LengthUnits()); // m -> GU
    TVector3 udir = p.Vect().Unit(); // unit vector along direction

    fCurrPathLengthList->AddPathLength(
                       pdgc, this->ComputePathLengthPDG(pos,udir,pdgc));
  }

  this->ScalePathLengths(*fCurrPathLengthList); // GU -> m

  return *fCurrPathLengthList;
}
//___________________________________________________________________________
const TVector3 & ROOTGeomAnalyzer::GenerateVertex(
              const TLorentzVector & x, const TLorentzVector & p, int tgtpdg)
{
// Generates a random vertex, within the detector material with the input
// PDG code, for a neutrino starting from point x and travelling along the
// direction of p

  LOG("GROOTGeom", pNOTICE)
         << "Generating vtx in material: " << tgtpdg
                                  << " along the input neutrino direction";
  // reset current interaction vertex
  fCurrVertex->SetXYZ(0.,0.,0.);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GROOTGeom", pDEBUG)
       << "\nInput nu: 4p (GeV) = " << utils::print::P4AsShortString(&p)
       << ", 4x (m,s) = " << utils::print::X4AsString(&x);
#endif

  if(!fGeometry) {
      LOG("GROOTGeom", pFATAL) << "No ROOT geometry is loaded!!";
      exit(1);
  }

  // calculate the event length for the selected material starting from
  // x and looking along the direction of p
  TVector3 r    = x.Vect();
  TVector3 dir  = p.Vect().Unit();
  r *= (1./this->LengthUnits());  // m -> GU
  double   dist = this->ComputePathLengthPDG(r, dir, tgtpdg);

  LOG("GROOTGeom", pNOTICE)
        << "Max {L(GU) x Density x Weight} given (init,dir) = " << dist;

  if(dist==0) {
    LOG("GROOTGeom", pERROR)
     << "The current trajectory does not cross the selected material!!";
    return *fCurrVertex;
  }

  // generate random number between 0 and dist
  RandomGen * rnd = RandomGen::Instance();
  double distVertex(dist * rnd->RndGeom().Rndm());
  LOG("GROOTGeom", pNOTICE)
       << "Generated 'distance' in selected material = " << distVertex;

  //-- generate the vertex

  TGeoVolume *   vol = 0;
  TGeoMedium *   med = 0;
  TGeoMaterial * mat = 0;

  int    FlagNotInYet(0);
  bool   condition(kTRUE);
  double StepIncrease(0.001/this->LengthUnits());
  double distToVtx(0);

  r.SetXYZ(x.X(), x.Y(), x.Z());
  r *= (1./this->LengthUnits());  // m -> GU

  fGeometry -> SetCurrentPoint (r[0],r[1],r[2]);

  while(((!FlagNotInYet) || condition) && distToVtx<distVertex) {

      condition=kTRUE;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
      LOG("GROOTGeom",pDEBUG)
           << "Position = " << utils::print::Vec3AsString(&r)
                             << ", flag(not in yet) = " << FlagNotInYet;
#endif
      r = r + StepIncrease * dir;
      fGeometry -> SetCurrentPoint (r[0],r[1],r[2]);
      fGeometry->FindNode();

      med = 0;
      mat = 0;
      vol = fGeometry->GetCurrentVolume();

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
      LOG("GROOTGeom", pDEBUG) << "Current volume: " << vol->GetName();
#endif
      if(fGeometry->IsOutside() || !vol) {
         condition=kFALSE;
         if(FlagNotInYet) break;
      }

      if(condition) {
         if(!FlagNotInYet) FlagNotInYet=1;
         mat = vol->GetMedium()->GetMaterial();
         double weight = this->GetWeight(mat,tgtpdg);
         distToVtx+=(StepIncrease*weight);
     }
  }

  r = r - StepIncrease * dir;
  fCurrVertex->SetXYZ(r[0],r[1],r[2]);
  (*fCurrVertex) *= (this->LengthUnits()); // GU -> m

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GROOTGeom", pDEBUG) << "Vtx (m) = " << utils::print::Vec3AsString(&r);
#endif

  return *fCurrVertex;
}
//___________________________________________________________________________
void ROOTGeomAnalyzer::BuildListOfTargetNuclei(void)
{
  fCurrPDGCodeList = new PDGCodeList;

  if(!fGeometry) {
    LOG("GROOTGeom", pFATAL) << "No ROOT geometry is loaded!!";
    exit(1);
  }

  TObjArray * volume_list = fGeometry->GetListOfVolumes();
  if(!volume_list) {
     LOG("GROOTGeom", pERROR)
        << "Null list of geometry volumes. Can not find build target list!";
     return;
  }

  int numVol = volume_list->GetEntries();
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GROOTGeom", pDEBUG) << "Number of volumes found: " << numVol;
#endif

  for(int ivol = 0; ivol < numVol; ivol++) {
      TGeoVolume * volume = dynamic_cast <TGeoVolume *>(volume_list->At(ivol));
      if(!volume) {
         LOG("GROOTGeom", pWARN)
           << "Got a null geometry volume!! Skiping current list element";
         continue;
      }
      TGeoMaterial * mat = volume->GetMedium()->GetMaterial();
      if(mat->IsMixture()) {
         TGeoMixture * mixt = dynamic_cast <TGeoMixture*> (mat);
         int Nelements = mixt->GetNelements();
         for(int i=0; i<Nelements; i++) {
            TGeoElement * ele = mixt->GetElement(i);
            int ion_pdgc = this->GetTargetPdgCode(ele);
            fCurrPDGCodeList->push_back(ion_pdgc);
         }
      } else {
          int ion_pdgc = this->GetTargetPdgCode(mat);
          fCurrPDGCodeList->push_back(ion_pdgc);
      }
  }
}
//___________________________________________________________________________
double ROOTGeomAnalyzer::ComputePathLengthPDG(
                        const TVector3 & r0, const TVector3 & udir, int pdgc)
{
// Compute the path length for the material with pdg-code = pdc, staring from
// the input position r and moving along the direction of the unit vector udir
//
  double pl = 0; // <-- path length (x density, if weight by density is ON)

  int    counterloop  (0);
  int    FlagNotInYet (0);
  bool   condition    (kTRUE);

  double step   = 0;
  double weight = 0;

  TGeoVolume *   vol = 0;
  TGeoMedium *   med = 0;
  TGeoMaterial * mat = 0;

  fGeometry -> SetCurrentDirection (udir[0],udir[1],udir[2]);
  fGeometry -> SetCurrentPoint     (r0[0],r0[1],r0[2]);

  while(((!FlagNotInYet) || condition) && counterloop <100) {
     counterloop++;
     condition=kTRUE;

     fGeometry->FindNode();

     med = 0;
     mat = 0;
     vol = fGeometry->GetCurrentVolume();

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("GROOTGeom", pDEBUG) << "Current volume: " << vol->GetName();
#endif
     if (fGeometry->IsOutside() || !vol) {
        condition=kFALSE;
        if(FlagNotInYet) break;

        step = this->StepToNextBoundary();
        while(!fGeometry->IsEntering()) {
          step = this->Step();
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
          LOG("GROOTGeom", pDEBUG) << "Stepping...dr = " << step;
#endif
          if(this->WillNeverEnter(step)) return 0.;
        }
      }

      if(condition) {
       if(!FlagNotInYet) FlagNotInYet=1;
       med = vol->GetMedium();
       mat = med->GetMaterial();

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
       LOG("GROOTGeom",pDEBUG)
         << "Cur med.: " << med->GetName() << ", mat.: " << mat->GetName();
#endif
       step   = this->StepUntilEntering();
       weight = this->GetWeight(mat, pdgc);

       pl += (step*weight);
     }//condition
  }

  LOG("GROOTGeom", pDEBUG) << "PathLength[" << pdgc << "] = " << pl << " GU";

  return pl;
}
//___________________________________________________________________________
double ROOTGeomAnalyzer::GetWeight(TGeoMaterial * mat, int pdgc)
{
// Get the weight of the input material.
// Return the weight only if the material's pdg code matches the input code.
// If the material is found to be a mixture, call the corresponding method
// for mixtures.

  if(!mat) {
    LOG("GROOTGeom", pERROR) << "Null input material. Return weight = 0.";
    return 0;
  }

  bool exists = fCurrPDGCodeList->ExistsInPDGCodeList(pdgc);
  if(!exists) {
    LOG("GROOTGeom", pERROR) << "Target doesn; exist. Return weight = 0.";
    return 0;
  }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GROOTGeom",pDEBUG)
       << "Curr. material: A/Z = " << mat->GetA() << " / " << mat->GetZ();
#endif

  // if the input material is a mixture, get a the sum of weights for
  // all matching elements
  double weight = 0.;
  if(mat->IsMixture()) {
    TGeoMixture * mixt = dynamic_cast <TGeoMixture*> (mat);

    if(!mixt) {
     LOG("GROOTGeom", pERROR) << "Null input mixture. Return weight = 0.";
     return 0;
    }
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("GROOTGeom", pDEBUG)
      << "Material : " << mat->GetName()
          << " is a mixture with " << mixt->GetNelements() << " elements";
#endif
    // loop over elements & sum weights of matching elements
    weight = this->GetWeight(mixt,pdgc);
    return weight;
  } // is mixture?

  // pure material
  int ion_pdgc = this->GetTargetPdgCode(mat);
  if(ion_pdgc != pdgc) return 0.;

  if (this->WeightWithDensity()) 
    weight = mat->GetDensity() / this->DensityUnits();
  else                           
    weight = 1.0;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GROOTGeom", pDEBUG)
       << "Weight[mat:" << mat->GetName() << "] = " << weight;
#endif

  return weight;
}
//___________________________________________________________________________
double ROOTGeomAnalyzer::GetWeight(TGeoMixture * mixt, int pdgc)
{
// Loop over the mixture elements, find the one matching the input pdgc
// and  return its weight

  double weight = 0;

  int nm = 0;
  for(int i = 0; i < mixt->GetNelements(); i++) {
     double dw = (this->GetWeight(mixt,i,pdgc));
     if(dw>0) nm++;
     weight += dw;
  }

  if(nm>1) {
     for(int j = 0; j < mixt->GetNelements(); j++) {
           TGeoElement * e = mixt->GetElement(j);
           LOG("GROOTGeom", pWARN)
              << "[" << j << "] Z = " << e->Z() << ", A = " << e->A()
                      << " (pdgc = " << this->GetTargetPdgCode(e)
                                       << "), w = " << mixt->GetWmixt()[j];
     }
     LOG("GROOTGeom", pERROR)
        << "Material pdgc = " << pdgc << " appears " << nm
                        << " times (>1) in mixture = " << mixt->GetName();
     LOG("GROOTGeom", pFATAL)
        << "Your geometry must be incorrect - Aborting";
     exit(1);
  }

  // if we are not weighting with the density then the weight=1 if the pdg
  // code was matched for any element of this mixture
  if( !this->WeightWithDensity() && weight>0. ) weight=1.0;

  return weight;
}
//___________________________________________________________________________
double ROOTGeomAnalyzer::GetWeight(TGeoMixture* mixt, int ielement, int pdgc)
{
// Get the weight of the input ith element of the input material.
// Return the weight only if the element's pdg code matches the input code
//
  int ion_pdgc = this->GetTargetPdgCode(mixt->GetElement(ielement));
  if(ion_pdgc != pdgc) return 0.;

  double d = mixt->GetDensity() / this->DensityUnits(); // mixture density
  double w = mixt->GetWmixt()[ielement]; // relative proportion by mass

  double wtot = this->MixtureWeightsSum();

  // <0 forces explicit calculation of relative proportion normalization
  if(wtot < 0) {
    wtot = 0;
    for(int i = 0; i < mixt->GetNelements(); i++) {
      wtot += (mixt->GetWmixt()[ielement]);
    }
  }
  assert(wtot>0);

  w /= wtot;
  double weight = d*w;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GROOTGeom", pDEBUG)
       << "Weight[mixt:" << mixt->GetName()
                      << ", iel = " << ielement << "] = " << weight;
#endif

  return weight;
}
//___________________________________________________________________________
double ROOTGeomAnalyzer::StepToNextBoundary(void)
{
  fGeometry->FindNextBoundary();
  double step=fGeometry->GetStep();
  return step;
}
//___________________________________________________________________________
double ROOTGeomAnalyzer::Step(void)
{
  fGeometry->Step();
  double step=fGeometry->GetStep();
  return step;
}
//___________________________________________________________________________
double ROOTGeomAnalyzer::StepUntilEntering(void)
{
  double step  = this->StepToNextBoundary();

  while(!fGeometry->IsEntering()) {
    step = this->Step();
  }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__

  bool isen = fGeometry->IsEntering();
  bool isob = fGeometry->IsOnBoundary();

  LOG("GROOTGeom",pDEBUG)
      << "IsEntering = "     << utils::print::BoolAsYNString(isen)
      << ", IsOnBoundary = " << utils::print::BoolAsYNString(isob)
      << ", Step = " << step;
#endif

  return step;
}
//___________________________________________________________________________
bool ROOTGeomAnalyzer::WillNeverEnter(double step)
{
// If the neutrino trajectory would never enter the detector, then the
// TGeoManager::GetStep returns the maximum step (1E30).
// Compare surrent step with max step and figure out whether the particle
// would never enter the detector

  if(step > 9.99E29) {

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("GROOTGeom", pINFO) << "Wow! Current step is dr = " << step;
     LOG("GROOTGeom", pINFO) << "This trajectory isn't entering the detector";
#endif
     return true;

  } else 
    return false;
}
//___________________________________________________________________________
void ROOTGeomAnalyzer::ScalePathLengths(PathLengthList & pl)
{
// convert path lenghts to default GENIE length scale
//
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GROOTGeom", pDEBUG)
    << "Scaling path-lengths -> meters (scale = " << fLengthScale << ")";
#endif

  PathLengthList::iterator pliter;
  for(pliter = pl.begin(); pliter != pl.end(); ++pliter)
  {
    int pdgc = pliter->first;
    pl.ScalePathLength(pdgc,fLengthScale);
  }
}
//___________________________________________________________________________
int ROOTGeomAnalyzer::GetTargetPdgCode(const TGeoMaterial * const m) const
{
  int A = TMath::Nint(m->GetA());
  int Z = TMath::Nint(m->GetZ());

  int pdgc = pdg::IonPdgCode(A,Z);

  return pdgc;
}
//___________________________________________________________________________
int ROOTGeomAnalyzer::GetTargetPdgCode(const TGeoElement * const e) const
{
  int A = TMath::Nint(e->A());
  int Z = e->Z();

  int pdgc = pdg::IonPdgCode(A,Z);

  return pdgc;
}
//___________________________________________________________________________

