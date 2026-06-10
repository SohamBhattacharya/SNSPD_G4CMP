/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// Run action — opens/closes phonon_steps.root using G4RootAnalysisManager.
// G4RootAnalysisManager writes ROOT binary format natively; no external
// ROOT installation required. Output is ~10x smaller than text and
// ~5x faster to write due to buffering and compression.
//
// NOTE: ntuple must be created AFTER OpenFile() to avoid uninitialized
// ROOT buffer crash (vector::reserve) on the first AddNtupleRow call.

#include "RISQTutorialRunAction.hh"
#include "G4RootAnalysisManager.hh"
#include "G4Run.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RISQTutorialRunAction::RISQTutorialRunAction()
{
  // Constructor intentionally empty.
  // Ntuple schema is defined in BeginOfRunAction, after OpenFile(),
  // so ROOT buffers are properly initialized before any Fill calls.
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RISQTutorialRunAction::~RISQTutorialRunAction() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RISQTutorialRunAction::BeginOfRunAction(const G4Run*)
{
  auto am = G4RootAnalysisManager::Instance();
  am->SetVerboseLevel(1);

  // Open file first — ntuple creation must follow OpenFile()
  am->OpenFile("phonon_steps");   // produces phonon_steps.root

  // ── Define TTree "Steps" (ntuple ID = 0) ──────────────────────────────────
  // ptype encoding:
  //   1=phononL  2=phononTF  3=phononTS
  //   4=G4CMPDriftElectron  5=G4CMPDriftHole
  //   6=e-  7=e+  8=gamma  9=proton  10=neutron
  //   11=pi+  12=pi-  13=pi0  14=kaon+  15=kaon-
  //   16=mu+  17=mu-  99=other
  //
  // procID encoding:
  //   1=phononScattering  2=phononDownconversion  3=phononReflection
  //   4=Transportation    5=UserMaxTime           0=other

  am->CreateNtuple("Steps", "Phonon and primary step data");

  // Integer columns
  am->CreateNtupleIColumn("run");         // 0
  am->CreateNtupleIColumn("event");       // 1
  am->CreateNtupleIColumn("trackID");     // 2
  am->CreateNtupleIColumn("ptype");       // 3

  // Pre-step point
  am->CreateNtupleDColumn("preX_mm");     // 4
  am->CreateNtupleDColumn("preY_mm");     // 5
  am->CreateNtupleDColumn("preZ_mm");     // 6
  am->CreateNtupleDColumn("preE_eV");     // 7
  am->CreateNtupleDColumn("preKE_eV");    // 8
  am->CreateNtupleDColumn("preT_ns");     // 9

  // Post-step point
  am->CreateNtupleDColumn("postX_mm");    // 10
  am->CreateNtupleDColumn("postY_mm");    // 11
  am->CreateNtupleDColumn("postZ_mm");    // 12
  am->CreateNtupleDColumn("postE_eV");    // 13
  am->CreateNtupleDColumn("postKE_eV");   // 14
  am->CreateNtupleDColumn("postT_ns");    // 15

  // Step quantities
  am->CreateNtupleDColumn("stepLen_mm");  // 16
  am->CreateNtupleDColumn("edep_eV");     // 17

  // Process identifier
  am->CreateNtupleIColumn("procID");      // 18

  am->FinishNtuple();

  G4cout << "### RunAction: phonon_steps.root opened" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RISQTutorialRunAction::EndOfRunAction(const G4Run*)
{
  auto am = G4RootAnalysisManager::Instance();
  am->Write();
  am->CloseFile();
  G4cout << "### RunAction: phonon_steps.root written and closed" << G4endl;
}
