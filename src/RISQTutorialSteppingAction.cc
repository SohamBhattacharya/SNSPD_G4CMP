/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// Stepping action — fills ROOT ntuple with every step and kills tracks
// whose post-step global time exceeds kMaxGlobalTime.

#include "RISQTutorialSteppingAction.hh"
#include "G4RootAnalysisManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "globals.hh"

// ── Kill threshold ────────────────────────────────────────────────────────────
// At ~5 mm/ns ballistic phonon speed in Si:
//   100 ns ≈ ~500 mm total path  (~1000 bounces across the 0.525 mm slab)
//    10 ns ≈ ~50  mm total path  (~100 bounces)
static constexpr G4double kMaxGlobalTime = 50.0 * CLHEP::ns;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RISQTutorialSteppingAction::RISQTutorialSteppingAction()
{
  G4cout << "### SteppingAction: writing steps to ROOT ntuple" << G4endl;
  G4cout << "### SteppingAction: killing tracks after "
         << kMaxGlobalTime / CLHEP::ns << " ns" << G4endl;
}

RISQTutorialSteppingAction::~RISQTutorialSteppingAction() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RISQTutorialSteppingAction::UserSteppingAction(const G4Step* step)
{
  G4String pname = step->GetTrack()
                       ->GetParticleDefinition()
                       ->GetParticleName();

  // Debug: print first occurrence of each particle type seen
  if (fSeenParticles.find(pname) == fSeenParticles.end()) {
    G4cout << "### SteppingAction: first step for particle: "
           << pname << G4endl;
    fSeenParticles.insert(pname);
  }

  // ── Map particle name → integer type ──────────────────────────────────────
  G4int ptype = 0;
  if      (pname == "phononL")             ptype = 1;
  else if (pname == "phononTF")            ptype = 2;
  else if (pname == "phononTS")            ptype = 3;
  else if (pname == "G4CMPDriftElectron")  ptype = 4;
  else if (pname == "G4CMPDriftHole")      ptype = 5;
  else if (pname == "e-")                  ptype = 6;
  else if (pname == "e+")                  ptype = 7;
  else if (pname == "gamma")               ptype = 8;
  else if (pname == "proton")              ptype = 9;
  else if (pname == "neutron")             ptype = 10;
  else if (pname == "pi+")                 ptype = 11;
  else if (pname == "pi-")                 ptype = 12;
  else if (pname == "pi0")                 ptype = 13;
  else if (pname == "kaon+")              ptype = 14;
  else if (pname == "kaon-")              ptype = 15;
  else if (pname == "mu+")                ptype = 16;
  else if (pname == "mu-")                ptype = 17;
  else                                     ptype = 99;

  // ── Collect step data ──────────────────────────────────────────────────────
  const G4StepPoint* pre  = step->GetPreStepPoint();
  const G4StepPoint* post = step->GetPostStepPoint();

  G4int runNo   = G4RunManager::GetRunManager()
                    ->GetCurrentRun()->GetRunID();
  G4int eventNo = G4RunManager::GetRunManager()
                    ->GetCurrentEvent()->GetEventID();
  G4int trackID = step->GetTrack()->GetTrackID();

  G4double preX  = pre->GetPosition().x()  / CLHEP::mm;
  G4double preY  = pre->GetPosition().y()  / CLHEP::mm;
  G4double preZ  = pre->GetPosition().z()  / CLHEP::mm;
  G4double preE  = pre->GetTotalEnergy()   / CLHEP::eV;
  G4double preKE = pre->GetKineticEnergy() / CLHEP::eV;
  G4double preT  = pre->GetGlobalTime()    / CLHEP::ns;

  G4double postX  = post->GetPosition().x()  / CLHEP::mm;
  G4double postY  = post->GetPosition().y()  / CLHEP::mm;
  G4double postZ  = post->GetPosition().z()  / CLHEP::mm;
  G4double postE  = post->GetTotalEnergy()   / CLHEP::eV;
  G4double postKE = post->GetKineticEnergy() / CLHEP::eV;
  G4double postT  = post->GetGlobalTime()    / CLHEP::ns;

  G4double stepLen = step->GetStepLength()         / CLHEP::mm;
  G4double edep    = step->GetTotalEnergyDeposit() / CLHEP::eV;

  // ── Map process name → integer ID ─────────────────────────────────────────
  G4int procID = 0;
  if (post->GetProcessDefinedStep()) {
    G4String proc = post->GetProcessDefinedStep()->GetProcessName();
    if      (proc.find("phononScattering")     != G4String::npos) procID = 1;
    else if (proc.find("phononDownconversion") != G4String::npos) procID = 2;
    else if (proc.find("phononReflection")     != G4String::npos) procID = 3;
    else if (proc.find("Transportation")       != G4String::npos) procID = 4;
    else if (proc.find("UserMaxTime")          != G4String::npos) procID = 5;
  }

  // ── Fill ntuple row (columns 0–18 match RunAction schema) ─────────────────
  auto am = G4RootAnalysisManager::Instance();
  const G4int id = 0;
  G4int col = 0;
  am->FillNtupleIColumn(id, col++, runNo);
  am->FillNtupleIColumn(id, col++, eventNo);
  am->FillNtupleIColumn(id, col++, trackID);
  am->FillNtupleIColumn(id, col++, ptype);
  am->FillNtupleDColumn(id, col++, preX);
  am->FillNtupleDColumn(id, col++, preY);
  am->FillNtupleDColumn(id, col++, preZ);
  am->FillNtupleDColumn(id, col++, preE);
  am->FillNtupleDColumn(id, col++, preKE);
  am->FillNtupleDColumn(id, col++, preT);
  am->FillNtupleDColumn(id, col++, postX);
  am->FillNtupleDColumn(id, col++, postY);
  am->FillNtupleDColumn(id, col++, postZ);
  am->FillNtupleDColumn(id, col++, postE);
  am->FillNtupleDColumn(id, col++, postKE);
  am->FillNtupleDColumn(id, col++, postT);
  am->FillNtupleDColumn(id, col++, stepLen);
  am->FillNtupleDColumn(id, col++, edep);
  am->FillNtupleIColumn(id, col++, procID);
  am->AddNtupleRow(id);

  // ── Time cut — kill track after threshold ──────────────────────────────────
  // Checked AFTER filling so the final step is still recorded.
  if (post->GetGlobalTime() >= kMaxGlobalTime) {
    step->GetTrack()->SetTrackStatus(fStopAndKill);
  }
}
