/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/
/// \file RISQTutorial/RISQTutorial.cc
/// \brief Main program of the RISQTutorial example

#include "G4RunManager.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4CMPPhysicsList.hh"
#include "G4CMPPhysics.hh"
#include "G4CMPConfigManager.hh"
#include "RISQTutorialActionInitialization.hh"
#include "RISQTutorialConfigManager.hh"
#include "RISQTutorialDetectorConstruction.hh"
#include "RISQTutorialDetectorParameters.hh"
#include "FTFP_BERT.hh"

using namespace RISQTutorialDetectorParameters;

int main(int argc, char** argv)
{
  // ── Run manager ────────────────────────────────────────────────────────────
  G4RunManager* runManager = new G4RunManager;

  runManager->SetUserInitialization(new RISQTutorialDetectorConstruction);

  FTFP_BERT* physics = new FTFP_BERT;
  physics->RegisterPhysics(new G4CMPPhysics);
  physics->SetCuts();
  runManager->SetUserInitialization(physics);

  runManager->SetUserInitialization(new RISQTutorialActionInitialization);

  G4CMPConfigManager::Instance();
  RISQTutorialConfigManager::Instance();

  runManager->Initialize();

  // ── Visualisation (only initialised — not opened yet) ─────────────────────
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // ── Batch mode: macro passed as argument ───────────────────────────────────
  if (argc > 1) {
    G4String command  = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
  }
  // ── Interactive mode: no argument → open Qt session ───────────────────────
  else {
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);

    UImanager->ApplyCommand("/vis/open OGL");
    UImanager->ApplyCommand("/vis/viewer/set/upVector 0 1 0");
    UImanager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 70 20");
    UImanager->ApplyCommand("/vis/viewer/zoom 1.4");
    UImanager->ApplyCommand("/vis/drawVolume");
    UImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate");
    UImanager->ApplyCommand("/vis/scene/add/trajectories");
    UImanager->ApplyCommand("/tracking/storeTrajectory 1");

    ui->SessionStart();
    delete ui;
  }

  delete visManager;
  delete runManager;
  return 0;
}
