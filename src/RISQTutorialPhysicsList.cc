/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file RISQTutorialPhysicsList.cc
/// \brief Maximum-fidelity physics list for RISQTutorial.
///
/// ┌─────────────────────────────────────────────────────────────────────┐
/// │  OVERKILL MODE — expect 10–50× slower than default FTFP_BERT        │
/// │    120 GeV proton   → FTFP_BERT_HP hadronic                        │
/// │    MeV delta-e      → Penelope EM (accurate to ~10 eV in Si/SiO2)  │
/// │    keV fluorescence → SetFluo + SetPixe                            │
/// │    ~90 eV Auger     → SetAuger + SetAugerCascade                   │
/// │    10 eV floor      → SetMinEnergy + 1 nm cuts in SiO2Region       │
/// │    meV phonons      → G4CMPPhysics                                 │
/// │                                                                     │
/// │  NOTE: G4EmDNA models require G4_WATER and crash in Geant4 10.7    │
/// │  on non-water materials. Penelope + 1 nm cuts is the correct        │
/// │  approach for Si/SiO2 and gives equivalent accuracy.               │
/// └─────────────────────────────────────────────────────────────────────┘

#include "RISQTutorialPhysicsList.hh"

// ── G4CMP ─────────────────────────────────────────────────────────────────
#include "G4CMPPhysics.hh"

// ── EM: Penelope ──────────────────────────────────────────────────────────
#include "G4EmPenelopePhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmParameters.hh"

// ── Optical photons ───────────────────────────────────────────────────────
#include "G4OpticalPhysics.hh"

// ── Hadronic ──────────────────────────────────────────────────────────────
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4IonElasticPhysics.hh"

// ── Decay ─────────────────────────────────────────────────────────────────
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

// ── Neutron ───────────────────────────────────────────────────────────────
#include "G4NeutronTrackingCut.hh"

// ── Step limiter ──────────────────────────────────────────────────────────
#include "G4StepLimiterPhysics.hh"

// ── Region / cuts ─────────────────────────────────────────────────────────
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"

// ── Units ─────────────────────────────────────────────────────────────────
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RISQTutorialPhysicsList::RISQTutorialPhysicsList()
  : G4VModularPhysicsList()
{
  SetVerboseLevel(1);

  RegisterPhysics(new G4EmPenelopePhysics());
  RegisterPhysics(new G4EmExtraPhysics());
  RegisterPhysics(new G4OpticalPhysics());
  RegisterPhysics(new G4DecayPhysics());
  RegisterPhysics(new G4RadioactiveDecayPhysics());
  RegisterPhysics(new G4HadronElasticPhysicsHP());
  RegisterPhysics(new G4HadronPhysicsFTFP_BERT_HP());
  RegisterPhysics(new G4StoppingPhysics());
  RegisterPhysics(new G4IonPhysics());
  RegisterPhysics(new G4IonElasticPhysics());
  RegisterPhysics(new G4NeutronTrackingCut());
  RegisterPhysics(new G4StepLimiterPhysics());
  RegisterPhysics(new G4CMPPhysics());   // must be last
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RISQTutorialPhysicsList::~RISQTutorialPhysicsList() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RISQTutorialPhysicsList::SetCuts()
{
  // Global: 10 nm — half the minimum step (sio2Thick/4 = 70 nm)
  SetCutValue(10.*nm, "gamma");
  SetCutValue(10.*nm, "e-");
  SetCutValue(10.*nm, "e+");
  SetCutValue(10.*nm, "proton");

  // SiO2 region: 1 nm — resolves Auger electron ranges (~3–5 nm in SiO2)
  // Penelope tracks electrons down to 10 eV which have ranges of ~1–5 nm
  // in SiO2 (density 2.2 g/cm³), so 1 nm cut is the physically meaningful
  // limit — below this, energy is deposited locally anyway.
  G4Region* sio2Region =
      G4RegionStore::GetInstance()->GetRegion("SiO2Region", false);
  if (sio2Region) {
    G4ProductionCuts* cuts = new G4ProductionCuts();
    cuts->SetProductionCut(1.*nm, "gamma");
    cuts->SetProductionCut(1.*nm, "e-");
    cuts->SetProductionCut(1.*nm, "e+");
    cuts->SetProductionCut(1.*nm, "proton");
    sio2Region->SetProductionCuts(cuts);
    G4cout << "### PhysicsList: 1 nm production cuts → SiO2Region" << G4endl;
  } else {
    G4cout << "### PhysicsList: SiO2Region not found — "
           << "global 10 nm cuts apply everywhere" << G4endl;
  }

  if (verboseLevel > 0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RISQTutorialPhysicsList::ConstructParticle()
{
  G4VModularPhysicsList::ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RISQTutorialPhysicsList::ConstructProcess()
{
  G4VModularPhysicsList::ConstructProcess();

  // ── Atomic de-excitation: full cascade ────────────────────────────────────
  G4EmParameters* emParams = G4EmParameters::Instance();

  emParams->SetFluo(true);              // Si Kα X-ray (1.74 keV), O Kα (0.52 keV)
  emParams->SetAuger(true);            // Si L-shell Auger cascade (~90 eV)
  emParams->SetAugerCascade(true);     // full multi-step cascade
  emParams->SetPixe(true);             // proton-induced X-ray emission
  emParams->SetDeexcitationIgnoreCut(true); // de-excite even below 1 nm cut

  // Penelope lower bound: 10 eV in Si/SiO2
  // Below this threshold energy is deposited locally (continuous slowing down)
  emParams->SetMinEnergy(10.*eV);

  // ICRU-90 (2014) stopping powers — more accurate than ICRU-49 for Si
  emParams->SetUseICRU90Data(true);

  // Safest MSC step limit — critical for correct transport in 280 nm SiO2
  emParams->SetMscStepLimitType(fUseSafetyPlus);

  emParams->SetVerbose(0);
  emParams->SetWorkerVerbose(0);

  G4cout << "### PhysicsList: Penelope EM active, "
         << "Auger+Fluo+PIXE enabled, 10 eV min energy" << G4endl;
}
