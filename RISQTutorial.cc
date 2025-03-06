/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file RISQTutorial/RISQTutorial.cc
/// \brief Main program of the RISQTutorial example (based on G4CMP's phonon example)
//
// $Id$
//
// 20140509  Add conditional code for Geant4 10.0 vs. earlier
// 20150112  Remove RM->Initialize() call to allow macro configuration
// 20160111  Remove Geant4 version check since we now hard depend on 10.2+
// 20170816  Add example-specific configuration manager
// 20220718  Remove obsolete pre-processor macros G4VIS_USE and G4UI_USE
// 20240521  Renamed for tutorial use

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

int main(int argc,char** argv)
{
	/*
    runManager->SetUserInitialization(new MyDetectorConstruction(G4Args));
	 G4VUserPhysicsList* physics = new G4CMPPhysicsList();
	 physics->SetCuts();
	 runManager->SetUserInitialization(physics);
    runManager->SetUserInitialization(new MyPhysicsList());
    runManager->SetUserInitialization(new MyActionInitialization(G4Args));
    runManager->Initialize();	
	*/
 // Construct the run manager
 //
 G4RunManager * runManager = new G4RunManager;

 // Set mandatory initialization classes
 //
 RISQTutorialDetectorConstruction* detector = new RISQTutorialDetectorConstruction();
 runManager->SetUserInitialization(detector);

 FTFP_BERT* physics = new FTFP_BERT;  
 physics->RegisterPhysics(new G4CMPPhysics);
 physics->SetCuts();
 runManager->SetUserInitialization(physics);
 
 // Set user action classes (different for Geant4 10.0)
 //
 runManager->SetUserInitialization(new RISQTutorialActionInitialization);

 // Create configuration managers to ensure macro commands exist
 G4CMPConfigManager::Instance();
 RISQTutorialConfigManager::Instance();
runManager->Initialize();	
 // Visualization manager
 //
     G4cout << "### Starting Vis: " << G4endl;

 G4VisManager* visManager = new G4VisExecutive;
 visManager->Initialize();
 
 
 // Get the pointer to the User Interface manager
 //
 G4UImanager* UImanager = G4UImanager::GetUIpointer();  


      G4UIExecutive * ui = new G4UIExecutive(argc,argv);

    // Open the visualization and set up the scene
    UImanager->ApplyCommand("/vis/open OGL");
    UImanager->ApplyCommand("/vis/viewer/set/viewpointVector -1 0 0");
    UImanager->ApplyCommand("/vis/drawVolume");
    UImanager->ApplyCommand("/vis/viewer/set/autoRefresh true");
    UImanager->ApplyCommand("/vis/scene/add/trajectories smooth");
    UImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate");
    UImanager->ApplyCommand("/vis/scene/add/axes");
	UImanager->ApplyCommand("/tracking/verbose 1");
	UImanager->ApplyCommand("/vis/scene/add/trajectories smooth");
	UImanager->ApplyCommand("/vis/modeling/trajectories/create smooth");
	UImanager->ApplyCommand("/vis/modeling/trajectories/drawByCharge");
	UImanager->ApplyCommand("/tracking/verbose 2");
	UImanager->ApplyCommand("/process/verbose 2");
	UImanager->ApplyCommand("/run/verbose 2");
      
      
      ui->SessionStart();

 return 0;
}


