/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/src/PhononDetectorConstruction.cc
/// \brief Implementation of the RISQTutorialDetectorConstruction class
//
// Simplified to a single 8x8x0.525 mm Si slab for phonon propagation studies.

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

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
  G4cout << "### SETUP GEOM: minimal Si slab (8x8x0.525 mm)" << G4endl;

  // ── Surface properties ────────────────────────────────────────────────────
  // Created once; G4CMP reuses them across geometry rebuilds.
  // New G4CMP API (post-2024):
  //   G4CMPSurfaceProperty(name,
  //     qAbsProb, qReflProb, eMinK, hMinK,       // charge parameters
  //     pAbsProb, pReflProb, pSpecProb, pMinK,   // phonon parameters
  //     qpAbsProb, qpReflProb)                   // quasiparticle parameters
  // Frequency-dependent scattering set via AddScatteringProperties().

  if (!fConstructed) {
    const G4double GHz = 1e9 * hertz;

    // Polynomial coefficients for frequency-dependent phonon scattering.
    // anhCoeffs: anharmonic decay (last term ~ f^5 Akhiezer-like scaling).
    // diffCoeffs / specCoeffs: empty → no frequency-dependent diffuse/specular.
    const std::vector<G4double> anhCoeffs  = {0., 0., 0., 0., 0., 1.51e-14};
    const std::vector<G4double> diffCoeffs = {};
    const std::vector<G4double> specCoeffs = {};

    // Cutoff frequencies in GHz (units supplied separately to AddScatteringProperties)
    const G4double anhCutoff  = 520.;  // anharmonic decay cutoff [GHz]
    const G4double reflCutoff = 350.;  // diffuse reflection cutoff [GHz]

    // Si / vacuum (world) interface:
    //   phonons reflect diffusely; nothing absorbed; QPs fully reflected.
    fSiVacuumInterface = new G4CMPSurfaceProperty("SiVacuumInterface",
        0.0, 1.0, 0.0, 0.0,   // qAbsProb, qReflProb, eMinK, hMinK
        0.0, 1.0, 0.0, 0.0,   // pAbsProb, pReflProb, pSpecProb, pMinK
        0.0, 1.0);             // qpAbsProb, qpReflProb

    fSiVacuumInterface->AddScatteringProperties(
        anhCutoff, reflCutoff,
        anhCoeffs, diffCoeffs, specCoeffs,
        GHz, GHz, GHz);
  }

  // ── World ─────────────────────────────────────────────────────────────────
  G4VSolid* solid_world = new G4Box("World", 5.*cm, 5.*cm, 5.*cm);
  G4LogicalVolume* log_world = new G4LogicalVolume(solid_world,
                                                    fLiquidHelium, "World");
  log_world->SetVisAttributes(G4VisAttributes::Invisible);
  fWorldPhys = new G4PVPlacement(0, G4ThreeVector(),
                                  log_world, "World", 0, false, 0);

  // ── Silicon slab: 10 x 10 x 0.525 mm ─────────────────────────────────────
  G4Box* solid_Si = new G4Box("SiSlab_solid",
                               5.000*mm,    // half X → 10 mm total
                               5.000*mm,    // half Y → 10 mm total
                               0.2625*mm);  // half Z → 0.525 mm total

  G4LogicalVolume* log_Si = new G4LogicalVolume(solid_Si, fSilicon,
                                                  "SiSlab_log");

  G4VPhysicalVolume* phys_Si = new G4PVPlacement(0, G4ThreeVector(),
                                                   log_Si, "SiSlab",
                                                   log_world, false, 0, true);

  G4VisAttributes* siVis = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.3));
  siVis->SetVisibility(true);
  log_Si->SetVisAttributes(siVis);

  // ── Time cut ──────────────────────────────────────────────────────────────
  // At ballistic phonon speeds in Si (~5000 m/s), 1 ns ≈ 5 µm.
  // Increase to e.g. 100*ns to follow phonons across the full 525 µm slab.
  G4UserLimits* timeCut = new G4UserLimits();
  timeCut->SetUserMaxTime(1.0 * CLHEP::nanosecond);
  log_Si->SetUserLimits(timeCut);

 // ── G4CMP lattice ─────────────────────────────────────────────────────────
  // Change fLatticeName + fLatticeMaterial to switch crystal.
  //
  // Available lattices in $G4LATTICEDATA (CrystalMaps/):
  //   "Si"         → fSilicon   (G4_Si)              ← default, use this
  //   "Ge"         → fGermanium (G4_Ge)
  //   "GaAs"       → GaAs material (build manually)
  //   "Al2O3"      → sapphire
  //   "SiO2_alpha" → quartz
  //   "SiO2_amorph"→ amorphous SiO2
  //   "Al"         → aluminium
  //   "Cu"         → copper
  //   "Nb"         → niobium
  //   "LiF"        → lithium fluoride
  //   "CaF2"       → calcium fluoride
  //   "CaWO4"      → calcium tungstate
  //
  // To switch: change latticeName and latticeMaterial to match each other.
  // The material must be the same object used for the physical volume (fSilicon here).

  const G4String latticeName     = "Si";   // ← change this to switch crystal
  G4Material*    latticeMaterial = fSilicon; // ← must match latticeName

  G4LatticeManager*  LM       = G4LatticeManager::GetLatticeManager();
  G4LatticeLogical*  log_lat  = LM->LoadLattice(latticeMaterial, latticeName);
  if (!log_lat) {
    G4Exception("RISQTutorialDetectorConstruction::SetupGeometry",
                "LatticeLoad", FatalException,
                ("Could not load lattice: " + latticeName).c_str());
  }
  G4LatticePhysical* phys_lat = new G4LatticePhysical(log_lat);
  phys_lat->SetMillerOrientation(1, 0, 0);
  LM->RegisterLattice(phys_Si, phys_lat);
  G4cout << "### Lattice loaded: " << latticeName
         << " for material " << latticeMaterial->GetName() << G4endl;
         
         
  // ── Si / world boundary: phonons reflect at all six faces ─────────────────
  new G4CMPLogicalBorderSurface("border_Si_world",
                                  phys_Si, fWorldPhys, fSiVacuumInterface);

  // ── Sensitive detector ────────────────────────────────────────────────────
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
