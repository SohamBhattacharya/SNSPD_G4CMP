/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/src/RISQTutorialDetectorConstruction.cc
/// \brief Two-layer detector: crystalline Si slab + amorphous SiO2 overlayer.
//
// Geometry (Z axis, top = positive Z):
//
//   +Z
//    |   [ SiO2 layer ]   280 nm  (proton enters here first)
//    |   [ Si   layer ]   525 µm
//   -Z
//
// Each layer has its own lattice. Phonons cross the Si/SiO2 interface
// via a fully transmissive G4CMPLogicalBorderSurface (pAbsProb=0,
// pReflProb=0 → G4CMP defaults to transmission at unspecified interfaces).
//
// To change a layer's crystal: uncomment the desired block in the
// "LAYER SELECTION" sections below and comment the others.
// latticeMaterial MUST match the G4LogicalVolume material for that layer.
//
// Available crystals:
//   Cubic (fully supported):  Si, Ge, Al, Cu, Nb, GaAs
//   Approximated as cubic:    SiO2_alpha, SiO2_amorph, Al2O3, CaWO4
//
// Run command (local CrystalMaps with cubic SiO2 approximation):
//   singularity run \
//     --env G4LATTICEDATA="...CrystalMaps" \
//     ~/ubuntu-sandbox/ ./RISQTutorial ../G4Macros/throwProton_batch.mac
//
// Run command (official CrystalMaps — Si only; SiO2 needs charge params):
//   singularity run \
//     --env G4LATTICEDATA="/opt/G4CMP/G4CMP-install/share/G4CMP/CrystalMaps" \
//     ~/ubuntu-sandbox/ ./RISQTutorial ../G4Macros/throwProton_batch.mac

#include "RISQTutorialDetectorConstruction.hh"
#include "RISQTutorialSensitivity.hh"
#include "G4CMPPhononElectrode.hh"
#include "G4CMPElectrodeSensitivity.hh"
#include "G4CMPLogicalBorderSurface.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4GeometryManager.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RISQTutorialDetectorConstruction::RISQTutorialDetectorConstruction()
  : fLiquidHelium(0), fGermanium(0), fAluminum(0), fTungsten(0),
    fSilicon(0), fNiobium(0),
    fWorldPhys(0),
    fSiVacuumInterface(0),
    fSuperconductorSensitivity(0),
    fConstructed(false)
{;}

RISQTutorialDetectorConstruction::~RISQTutorialDetectorConstruction() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* RISQTutorialDetectorConstruction::Construct()
{
  if (fConstructed) {
    if (!G4RunManager::IfGeometryHasBeenDestroyed()) {
      G4GeometryManager::GetInstance()->OpenGeometry();
      G4PhysicalVolumeStore::GetInstance()->Clean();
      G4LogicalVolumeStore::GetInstance()->Clean();
      G4SolidStore::GetInstance()->Clean();
    }
    G4LatticeManager::GetLatticeManager()->Reset();
    G4CMPLogicalBorderSurface::CleanSurfaceTable();
  }

  DefineMaterials();
  SetupGeometry();
  fConstructed = true;
  return fWorldPhys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RISQTutorialDetectorConstruction::DefineMaterials()
{
  G4NistManager* nist = G4NistManager::Instance();
  fLiquidHelium = nist->FindOrBuildMaterial("G4_AIR");
  fGermanium    = nist->FindOrBuildMaterial("G4_Ge");
  fSilicon      = nist->FindOrBuildMaterial("G4_Si");
  fAluminum     = nist->FindOrBuildMaterial("G4_Al");
  fTungsten     = nist->FindOrBuildMaterial("G4_W");
  fNiobium      = nist->FindOrBuildMaterial("G4_Nb");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RISQTutorialDetectorConstruction::SetupGeometry()
{
  // ── Geometry parameters ───────────────────────────────────────────────────
  // Both layers share the same XY footprint.
  const G4double halfXY    = 5.000*mm;      // 10 x 10 mm lateral size

  const G4double siThick   = 525.0*um;      // Si substrate thickness
  const G4double sio2Thick =   0.280*um;    // SiO2 overlayer thickness (280 nm)

  const G4double halfSi    = siThick   / 2.;
  const G4double halfSiO2  = sio2Thick / 2.;

  // Z positions (Si bottom face at -halfSi, SiO2 top face at halfSi+sio2Thick)
  // Si   centre: z = 0
  // SiO2 centre: z = halfSi + halfSiO2
  const G4double zSi   = 0.;
  const G4double zSiO2 = halfSi + halfSiO2;

  // Proton start: just above the SiO2 top face
  // Print for macro reference:
  G4cout << "### SETUP GEOM: Si(" << siThick/um << " um) + SiO2("
         << sio2Thick/nm << " nm)" << G4endl;
  G4cout << "### Proton should start above z = "
         << (zSiO2 + halfSiO2)/mm << " mm" << G4endl;

  // ── Surface properties ────────────────────────────────────────────────────
  if (!fConstructed) {
    const G4double GHz        = 1e9 * hertz;
    const G4double anhCutoff  = 520.;
    const G4double reflCutoff = 350.;
    const std::vector<G4double> anhCoeffs  = {0., 0., 0., 0., 0., 1.51e-14};
    const std::vector<G4double> diffCoeffs = {};
    const std::vector<G4double> specCoeffs = {};

    // Outer (slab/world) surfaces: phonons reflect diffusely.
    fSiVacuumInterface = new G4CMPSurfaceProperty("SlabVacuumInterface",
        0.0, 1.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 1.0);
    fSiVacuumInterface->AddScatteringProperties(
        anhCutoff, reflCutoff,
        anhCoeffs, diffCoeffs, specCoeffs,
        GHz, GHz, GHz);
  }

  // ══════════════════════════════════════════════════════════════════════════
  // LAYER 1 — Si substrate (bottom layer, proton exits through it)
  // Uncomment exactly ONE block. latticeMaterialSi must match the volume.
  // ══════════════════════════════════════════════════════════════════════════

  // ── 1a. Silicon (default, recommended) ────────────────────────────────────
  const G4String latticeNameSi   = "Si";
  G4Material*    latticeMaterialSi = fSilicon;

  // ── 1b. Germanium ─────────────────────────────────────────────────────────
  /*
  const G4String latticeNameSi     = "Ge";
  G4Material*    latticeMaterialSi = fGermanium;
  */

  // ── 1c. Aluminium ─────────────────────────────────────────────────────────
  /*
  const G4String latticeNameSi     = "Al";
  G4Material*    latticeMaterialSi = fAluminum;
  */

  // ── 1d. Niobium ───────────────────────────────────────────────────────────
  /*
  const G4String latticeNameSi     = "Nb";
  G4Material*    latticeMaterialSi = fNiobium;
  */

  // ══════════════════════════════════════════════════════════════════════════
  // LAYER 2 — SiO2 overlayer (top layer, proton enters here first)
  // Uncomment exactly ONE block. latticeMaterialSiO2 must match the volume.
  // ══════════════════════════════════════════════════════════════════════════

  // ── 2a. Amorphous SiO2 (default — thermal oxide approximation) ────────────
  //    Cubic KV table approximation. Physically correct sound speeds.
  const G4String latticeNameSiO2     = "SiO2_amorph";
  G4Material*    latticeMaterialSiO2 = G4NistManager::Instance()
                                         ->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

  // ── 2b. Alpha-quartz SiO2 ─────────────────────────────────────────────────
  //    Cubic KV table approximation. Correct for crystalline quartz substrates.
  /*
  const G4String latticeNameSiO2     = "SiO2_alpha";
  G4Material*    latticeMaterialSiO2 = G4NistManager::Instance()
                                         ->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  */

  // ── 2c. No SiO2 layer — set thickness to 0 in geometry params above ───────
  //    (Comment out both volume and lattice blocks for SiO2 if not needed)

  // ══════════════════════════════════════════════════════════════════════════

  // ── World ─────────────────────────────────────────────────────────────────
  G4VSolid*        solid_world = new G4Box("World", 5.*cm, 5.*cm, 5.*cm);
  G4LogicalVolume* log_world   = new G4LogicalVolume(solid_world,
                                                      fLiquidHelium, "World");
  log_world->SetVisAttributes(G4VisAttributes::Invisible);
  fWorldPhys = new G4PVPlacement(0, G4ThreeVector(),
                                  log_world, "World", 0, false, 0);

  // ── Si substrate ──────────────────────────────────────────────────────────
  G4Box* solid_Si = new G4Box("SiSolid", halfXY, halfXY, halfSi);

  G4LogicalVolume* log_Si = new G4LogicalVolume(solid_Si,
                                                  latticeMaterialSi,
                                                  "SiLog");

  G4VPhysicalVolume* phys_Si = new G4PVPlacement(
      0, G4ThreeVector(0., 0., zSi),
      log_Si, "SiSlab", log_world, false, 0, true);

  log_Si->SetVisAttributes(
      new G4VisAttributes(G4Colour(0.5, 0.5, 1.0, 0.4)));  // blue-ish

  G4UserLimits* timeCutSi = new G4UserLimits();
  timeCutSi->SetUserMaxTime(100.0 * CLHEP::nanosecond);
  log_Si->SetUserLimits(timeCutSi);

  // ── SiO2 overlayer ────────────────────────────────────────────────────────
  G4Box* solid_SiO2 = new G4Box("SiO2Solid", halfXY, halfXY, halfSiO2);

  G4LogicalVolume* log_SiO2 = new G4LogicalVolume(solid_SiO2,
                                                    latticeMaterialSiO2,
                                                    "SiO2Log");

  G4VPhysicalVolume* phys_SiO2 = new G4PVPlacement(
      0, G4ThreeVector(0., 0., zSiO2),
      log_SiO2, "SiO2Layer", log_world, false, 0, true);

  log_SiO2->SetVisAttributes(
      new G4VisAttributes(G4Colour(1.0, 0.8, 0.3, 0.6)));  // gold

  G4UserLimits* timeCutSiO2 = new G4UserLimits();
  timeCutSiO2->SetUserMaxTime(100.0 * CLHEP::nanosecond);
  log_SiO2->SetUserLimits(timeCutSiO2);

  // ── Load and register Si lattice ──────────────────────────────────────────
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();

  G4LatticeLogical* log_lat_Si = LM->LoadLattice(latticeMaterialSi, latticeNameSi);
  if (!log_lat_Si) {
    G4Exception("RISQTutorialDetectorConstruction::SetupGeometry",
                "LatticeLoadSi", FatalException,
                ("Could not load Si lattice: " + latticeNameSi).c_str());
  }
  G4LatticePhysical* phys_lat_Si = new G4LatticePhysical(log_lat_Si);
  phys_lat_Si->SetMillerOrientation(1, 0, 0);
  LM->RegisterLattice(phys_Si, phys_lat_Si);
  G4cout << "### Si layer lattice:   " << latticeNameSi
         << "  material: " << latticeMaterialSi->GetName() << G4endl;

  // ── Load and register SiO2 lattice ────────────────────────────────────────
  G4LatticeLogical* log_lat_SiO2 = LM->LoadLattice(latticeMaterialSiO2,
                                                      latticeNameSiO2);
  if (!log_lat_SiO2) {
    G4Exception("RISQTutorialDetectorConstruction::SetupGeometry",
                "LatticeLoadSiO2", FatalException,
                ("Could not load SiO2 lattice: " + latticeNameSiO2).c_str());
  }
  G4LatticePhysical* phys_lat_SiO2 = new G4LatticePhysical(log_lat_SiO2);
  phys_lat_SiO2->SetMillerOrientation(1, 0, 0);
  LM->RegisterLattice(phys_SiO2, phys_lat_SiO2);
  G4cout << "### SiO2 layer lattice: " << latticeNameSiO2
         << "  material: " << latticeMaterialSiO2->GetName() << G4endl;

  // ── Boundary surfaces ─────────────────────────────────────────────────────
  // Si bottom face / world: diffuse reflection
  new G4CMPLogicalBorderSurface("border_Si_world_bot",
                                  phys_Si, fWorldPhys, fSiVacuumInterface);

  // SiO2 top face / world: diffuse reflection
  new G4CMPLogicalBorderSurface("border_SiO2_world_top",
                                  phys_SiO2, fWorldPhys, fSiVacuumInterface);

  // Si/SiO2 interface: fully transmissive — phonons cross freely.
  // G4CMP transmits phonons by default when no absorption/reflection is set.
  // pAbsProb=0, pReflProb=0, pSpecProb=0 → full transmission.
  G4CMPSurfaceProperty* siSio2Interface = new G4CMPSurfaceProperty(
      "Si_SiO2_Interface",
      0.0, 0.0, 0.0, 0.0,   // qAbsProb, qReflProb, eMinK, hMinK
      0.0, 0.0, 0.0, 0.0,   // pAbsProb, pReflProb, pSpecProb, pMinK
      0.0, 0.0);             // qpAbsProb, qpReflProb

  // Register both directions so phonons cross in either direction
  new G4CMPLogicalBorderSurface("border_Si_to_SiO2",
                                  phys_Si, phys_SiO2, siSio2Interface);
  new G4CMPLogicalBorderSurface("border_SiO2_to_Si",
                                  phys_SiO2, phys_Si, siSio2Interface);

  // ── Sensitive detector on Si bottom face ──────────────────────────────────
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  if (!fSuperconductorSensitivity)
    fSuperconductorSensitivity = new RISQTutorialSensitivity("PhononElectrode");
  SDman->AddNewDetector(fSuperconductorSensitivity);
  log_Si->SetSensitiveDetector(fSuperconductorSensitivity);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RISQTutorialDetectorConstruction::AttachPhononSensor(
    G4CMPSurfaceProperty* surfProp)
{
  if (!surfProp) return;

  auto sensorProp = surfProp->GetPhononMaterialPropertiesTablePointer();
  sensorProp->AddConstProperty("filmAbsorption",      0.0);
  sensorProp->AddConstProperty("filmThickness",       90.*CLHEP::nm);
  sensorProp->AddConstProperty("gapEnergy",           1.6e-3*CLHEP::eV);
  sensorProp->AddConstProperty("lowQPLimit",          3.);
  sensorProp->AddConstProperty("phononLifetime",      4.17*CLHEP::ps);
  sensorProp->AddConstProperty("phononLifetimeSlope", 0.29);
  sensorProp->AddConstProperty("vSound",              3.480*CLHEP::km/CLHEP::s);
  sensorProp->AddConstProperty("subgapAbsorption",    0.0);

  surfProp->SetPhononElectrode(new G4CMPPhononElectrode);
}
