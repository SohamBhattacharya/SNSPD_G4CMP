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

// Function to set up Geant4 visualization settings
void SetupVisualization(G4UImanager* UImanager) {
    // Open the OpenGL visualization driver
    UImanager->ApplyCommand("/vis/open OGL");

    // Set the up vector for the viewer (defines the orientation of the view)
    UImanager->ApplyCommand("/vis/viewer/set/upVector 0 1 0");

    // Set the viewpoint vector (defines the direction the camera is facing)
    UImanager->ApplyCommand("/vis/viewer/set/viewpointVector 0 0 1");

    // Set the zoom level for the viewer
    UImanager->ApplyCommand("/vis/viewer/zoom 1.4");

    // Draw the entire volume in the visualization
    UImanager->ApplyCommand("/vis/drawVolume");

    // Set the verbosity level for tracking (2 for detailed output)
    UImanager->ApplyCommand("/tracking/verbose 2");

    // Store trajectories during the tracking process
    UImanager->ApplyCommand("/tracking/storeTrajectory 1");

    // Set the action for the scene when the event ends (accumulate scene data)
    UImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate");

    // Add trajectory visualization to the scene
    UImanager->ApplyCommand("/vis/scene/add/trajectories");

    // Create trajectory visualization by particle ID (different trajectories for different particles)
    UImanager->ApplyCommand("/vis/modeling/trajectories/create/drawByParticleID");

    // Set specific colors for different particles in the trajectory visualization
    UImanager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set mu- White");  // Muon: White
    UImanager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set e- Yellow");  // Electron: Yellow
    UImanager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set G4CMPDriftElectron Violet");  // Drift Electron: Violet
    UImanager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set G4CMPDriftHole Orange");  // Drift Hole: Orange
    UImanager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set phononTS Red");  // Phonon (TS): Red
    UImanager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set phononTF Green");  // Phonon (TF): Green
    UImanager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set phononL Blue");  // Phonon (L): Blue
    UImanager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set proton Cyan");  // Proton: Cyan

    // Set the viewer's style to wireframe for clearer visualization
    UImanager->ApplyCommand("/vis/viewer/set/style wireframe");

    // Hide marker to simplify the view
    UImanager->ApplyCommand("/vis/viewer/set/hiddenMarker true");

    // Set the viewpoint angles (theta and phi) for the camera
    UImanager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 0 90");

    // Set a second zoom level for the viewer
    UImanager->ApplyCommand("/vis/viewer/zoom 1.6");
}


// Function to configure G4CMP settings
void SetupG4CMP(G4UImanager* UImanager) {
    UImanager->ApplyCommand("/g4cmp/producePhonons 1.0");
    UImanager->ApplyCommand("/g4cmp/sampleLuke 1.0");
    UImanager->ApplyCommand("/g4cmp/produceCharges 1.0");
}

// Function to configure particle generation settings
void SetupMuonGeneration(G4UImanager* UImanager) {
    UImanager->ApplyCommand("/gps/number 1");
    UImanager->ApplyCommand("/gps/particle mu-");
    UImanager->ApplyCommand("/gps/pos/type Point");
    UImanager->ApplyCommand("/gps/direction 0 -4 -1");
    UImanager->ApplyCommand("/gps/pos/centre 0.0 0.5 0.56 cm");  // A bit above the chip top
    // UImanager->ApplyCommand("/gps/pos/centre 0.0 0.1 0.481 cm");  // halfway between bottom and top of chip (Uncomment if needed)
    UImanager->ApplyCommand("/gps/ene/type Mono");
    UImanager->ApplyCommand("/gps/energy 4 GeV");
}


// Function to configure particle generation settings
void SetupPhononGeneration(G4UImanager* UImanager) {
	
    UImanager->ApplyCommand("/gps/number 1");
    UImanager->ApplyCommand("/gps/particle phononL");
    UImanager->ApplyCommand("/gps/pos/type Point");
    UImanager->ApplyCommand("/gps/direction 0 -4 -1");
    UImanager->ApplyCommand("/gps/pos/centre 0.0 0.1 0.481 cm");  // halfway between bottom and top of chip
    // UImanager->ApplyCommand("/gps/pos/centre 0.0 0.1 0.481 cm");  // halfway between bottom and top of chip (Uncomment if needed)
    UImanager->ApplyCommand("/gps/ene/type Mono");
    UImanager->ApplyCommand("/gps/energy 0.03 eV");
}

using namespace RISQTutorialDetectorParameters;

int main(int argc,char** argv)
{
	
G4Args = new MyG4Args(mainargc, mainargv);

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
    G4VisManager* visManager = new G4VisExecutive();
    visManager->Initialize();
    
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4UIExecutive* ui = nullptr;

 
    if (argc == 1)  // Define UI session for interactive mode
    {
        ui = new G4UIExecutive(argc, argv);
        SetupVisualization(UImanager);
        SetupG4CMP(UImanager);
        SetupPhononGeneration(UImanager);
        
        // Check for -muon argument -- not working
        if (std::find(argv, argv + argc, std::string("-muon")) != argv + argc)
        {
        SetupMuonGeneration(UImanager);
        }

        // Check for -phonon argument -- not working
        if (std::find(argv, argv + argc, std::string("-phonon")) != argv + argc)
        {
        SetupPhononGeneration(UImanager);
        }

        ui->SessionStart();
    }
    else  // Batch mode
    {
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    }
 

 delete visManager;
 delete runManager;

 return 0;
}





