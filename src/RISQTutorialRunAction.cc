/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// Run action — opens/closes phonon_steps.root using G4RootAnalysisManager.
// G4RootAnalysisManager writes ROOT binary format natively; no external
// ROOT installation required. Output is ~10x smaller than text and
// ~5x faster to write due to buffering and compression.

#include "RISQTutorialRunAction.hh"
#include "g4root.hh" // <-- CHANGED: Legacy Geant4 10 ROOT Analysis header
#include "G4Run.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RISQTutorialRunAction::RISQTutorialRunAction()
{
  auto am = G4RootAnalysisManager::Instance(); // <-- CHANGED: Root Instance
  am->SetVerboseLevel(1);
  
  // am->SetNtupleMerging(true); // <-- CHANGED: Commented out for Geant4 10 compatibility

  // ── Define TTree "Steps" (ntuple ID = 0) ──────────────────────────────────
  // Columns match the 17-column text format exactly, stored as binary.
  // Integer columns
  am->CreateNtuple("Steps", "Phonon and primary step data");
  am->CreateNtupleIColumn("run");        // 0
  am->CreateNtupleIColumn("event");      // 1
  am->CreateNtupleIColumn("trackID");    // 2
  // ptype: 1=phononL  2=phononTF  3=phononTS  0=other/primary
  am->CreateNtupleIColumn("ptype");      // 3

  // Pre-step
  am->CreateNtupleDColumn("preX_mm");    // 4
  am->CreateNtupleDColumn("preY_mm");    // 5
  am->CreateNtupleDColumn("preZ_mm");    // 6
  am->CreateNtupleDColumn("preE_eV");    // 7
  am->CreateNtupleDColumn("preKE_eV");   // 8
  am->CreateNtupleDColumn("preT_ns");    // 9

  // Post-step
  am->CreateNtupleDColumn("postX_mm");   // 10
  am->CreateNtupleDColumn("postY_mm");   // 11
  am->CreateNtupleDColumn("postZ_mm");   // 12
  am->CreateNtupleDColumn("postE_eV");   // 13
  am->CreateNtupleDColumn("postKE_eV");  // 14
  am->CreateNtupleDColumn("postT_ns");   // 15

  // Step quantities
  am->CreateNtupleDColumn("stepLen_mm"); // 16
  am->CreateNtupleDColumn("edep_eV");    // 17
  // procID: 1=scatter 2=downconv 3=reflect 4=transport 5=timekill
  am->CreateNtupleIColumn("procID");     // 18

  am->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RISQTutorialRunAction::~RISQTutorialRunAction() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RISQTutorialRunAction::BeginOfRunAction(const G4Run*)
{
  auto am = G4RootAnalysisManager::Instance(); // <-- CHANGED: Root Instance
  am->OpenFile("phonon_steps");   // produces phonon_steps.root
  G4cout << "### RunAction: phonon_steps.root opened" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RISQTutorialRunAction::EndOfRunAction(const G4Run*)
{
  auto am = G4RootAnalysisManager::Instance(); // <-- CHANGED: Root Instance
  am->Write();
  am->CloseFile();
  G4cout << "### RunAction: phonon_steps.root written and closed" << G4endl;
}
