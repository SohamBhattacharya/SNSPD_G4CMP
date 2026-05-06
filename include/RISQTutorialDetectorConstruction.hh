/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/include/PhononDetectorConstruction.hh
/// \brief Definition of the RISQTutorialDetectorConstruction class
//
// $Id: 4c06153e9ea08f2a90b22c53e5c39bde4b847c07 $
//

#ifndef RISQTutorialDetectorConstruction_h
#define RISQTutorialDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4NistManager.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4UserLimits.hh"
#include "G4SubtractionSolid.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
//#include "event.hh"
#include "G4Trap.hh"
#include "G4UnionSolid.hh"
//#include "G4Args.hh"
#include "G4TriangularFacet.hh"
#include "G4TessellatedSolid.hh"
#include "G4QuadrangularFacet.hh"
#include "G4Tubs.hh"


#include "RISQTutorialSensitivity.hh"
#include "RISQTutorialQubitHousing.hh" 
#include "RISQTutorialPad.hh"
#include "RISQTutorialTransmissionLine.hh"
#include "RISQTutorialStraightFluxLine.hh"
#include "RISQTutorialCornerFluxLine.hh"
#include "RISQTutorialResonatorAssembly.hh"
#include "G4CMPPhononElectrode.hh"
#include "G4CMPElectrodeSensitivity.hh"
#include "G4CMPLogicalBorderSurface.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4FieldManager.hh"
#include "G4GeometryManager.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4CMPLogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4Tubs.hh"
#include "G4UniformMagField.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"


class G4Material;
class G4VPhysicalVolume;
class G4CMPSurfaceProperty;
class G4CMPElectrodeSensitivity;

class RISQTutorialDetectorConstruction : public G4VUserDetectorConstruction {
public:
  RISQTutorialDetectorConstruction();
  virtual ~RISQTutorialDetectorConstruction();
  
public:
  virtual G4VPhysicalVolume* Construct();
  
private:
  void DefineMaterials();
  void SetupGeometry();
  void AttachPhononSensor(G4CMPSurfaceProperty * surfProp);

  
private:
  G4Material* fLiquidHelium;
  G4Material* fGermanium;
  G4Material* fAluminum;
  G4Material* fTungsten;
  G4Material* fSilicon;
  G4Material* fNiobium;
  G4VPhysicalVolume* fWorldPhys;

  G4CMPSurfaceProperty* fSiNbInterface;
  G4CMPSurfaceProperty* fSiCopperInterface;
  G4CMPSurfaceProperty* fSiVacuumInterface;
  
      G4Box *solidWorld, *solidRadiator, *solidDetector, *solidScintillator, *solidAir, *solidCu1, *solidCu2, *solidAl1, *solidAl2, *solid_strip, *solid_aSistrip;
    G4LogicalVolume *logicWorld, *logicRadiator, *logicDetector, *logicScintillator, *logicAir[10], *logicCu1, *logicCu2, *logicAl1, *logicAl2, *logicsubstrate, *logicCompositeDetector, *logicComposite_aSiStrip;
    G4VPhysicalVolume *physWorld, *physDetector, *physRadiator, *physScintillator, *physAir[10], *physCu1, *physCu2, *physAl1, *physAl2;


    G4Material *SiO2, *H2O, *Aerogel, *worldMat, *NaI, *Air[10], *AirMat, *CuMat, *AlMat, *aSiMat, *WSiMat, *Si3N4Mat, *VacuumMat;
    G4Element *C, *Na, *I;

  
  G4CMPElectrodeSensitivity* fSuperconductorSensitivity;
  G4bool fConstructed;
  //G4bool fIfField;
  
  //public:
  //inline void Field(G4bool bl) { fIfField = bl; }
};

#endif

