/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// Stepping action — exports full per-step information for all phonon tracks.

#include "RISQTutorialSteppingAction.hh"
#include <iostream>
#include "globals.hh"
#include "G4Run.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Threading.hh"
#include "G4RunManager.hh"
#include "G4StepPoint.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RISQTutorialSteppingAction::RISQTutorialSteppingAction()
{
  // Open output file and write header
  fOutputFile.open("phonon_steps.txt", std::ios::trunc);
  fOutputFile << "# run  event  track  particle  "
              << "preX_mm  preY_mm  preZ_mm  preE_eV  preKE_eV  "
              << "postX_mm  postY_mm  postZ_mm  postE_eV  postKE_eV  "
              << "preT_ns  postT_ns  "
              << "process\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RISQTutorialSteppingAction::~RISQTutorialSteppingAction()
{
  fOutputFile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RISQTutorialSteppingAction::UserSteppingAction(const G4Step* step)
{
  // Only record phonon tracks
  G4String pname = step->GetTrack()->GetParticleDefinition()->GetParticleName();
  if (pname != "phononL" && pname != "phononTF" && pname != "phononTS") return;

  ExportStepInformation(step);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RISQTutorialSteppingAction::ExportStepInformation(const G4Step* step)
{
  G4StepPoint* preSP  = step->GetPreStepPoint();
  G4StepPoint* postSP = step->GetPostStepPoint();

  int runNo    = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
  int eventNo  = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  int trackNo  = step->GetTrack()->GetTrackID();

  std::string particleName =
      step->GetTrack()->GetParticleDefinition()->GetParticleName();

  // Pre-step
  double preX  = preSP->GetPosition().x() / CLHEP::mm;
  double preY  = preSP->GetPosition().y() / CLHEP::mm;
  double preZ  = preSP->GetPosition().z() / CLHEP::mm;
  double preE  = preSP->GetTotalEnergy()  / CLHEP::eV;
  double preKE = preSP->GetKineticEnergy()/ CLHEP::eV;
  double preT  = preSP->GetGlobalTime()   / CLHEP::ns;

  // Post-step
  double postX  = postSP->GetPosition().x() / CLHEP::mm;
  double postY  = postSP->GetPosition().y() / CLHEP::mm;
  double postZ  = postSP->GetPosition().z() / CLHEP::mm;
  double postE  = postSP->GetTotalEnergy()  / CLHEP::eV;
  double postKE = postSP->GetKineticEnergy()/ CLHEP::eV;
  double postT  = postSP->GetGlobalTime()   / CLHEP::ns;

  // Process that ended this step
  std::string process = "unknown";
  if (postSP->GetProcessDefinedStep())
    process = postSP->GetProcessDefinedStep()->GetProcessName();

  fOutputFile
      << runNo      << " "
      << eventNo    << " "
      << trackNo    << " "
      << particleName << " "
      << preX       << " " << preY  << " " << preZ  << " "
      << preE       << " " << preKE << " "
      << postX      << " " << postY << " " << postZ << " "
      << postE      << " " << postKE << " "
      << preT       << " " << postT << " "
      << process    << "\n";
}
