/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/include/PhononDetectorConstruction.hh
/// \brief Definition of the RISQTutorialDetectorConstruction class
//
// Simplified to a single Si slab for phonon propagation studies.

#ifndef RISQTutorialDetectorConstruction_h
#define RISQTutorialDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Material;
class G4VPhysicalVolume;
class G4CMPSurfaceProperty;
class G4CMPElectrodeSensitivity;

class RISQTutorialDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  RISQTutorialDetectorConstruction();
  virtual ~RISQTutorialDetectorConstruction();

  virtual G4VPhysicalVolume* Construct();

private:
  void DefineMaterials();
  void SetupGeometry();
  void AttachPhononSensor(G4CMPSurfaceProperty* surfProp);

  // Materials
  G4Material* fLiquidHelium;
  G4Material* fGermanium;
  G4Material* fAluminum;
  G4Material* fTungsten;
  G4Material* fSilicon;
  G4Material* fNiobium;

  // World physical volume
  G4VPhysicalVolume* fWorldPhys;

  // Surface property for Si/vacuum boundary
  G4CMPSurfaceProperty* fSiVacuumInterface;

  // Sensitive detector
  G4CMPElectrodeSensitivity* fSuperconductorSensitivity;

  G4bool fConstructed;
};

#endif
