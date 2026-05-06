/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "Physics.hh"
#include "G4CMPPhysicsList.hh"
#include "G4CMPPhysics.hh"
#include "G4CMPConfigManager.hh"
#include "RISQTutorialActionInitialization.hh"
#include "RISQTutorialConfigManager.hh"
#include "RISQTutorialDetectorConstruction.hh"
#include "RISQTutorialDetectorParameters.hh"
#include "FTFP_BERT.hh"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

// Function to execute a macro file
void ExecuteMacroFile(G4UImanager* UImanager, const std::string& macroFile) {
    std::ifstream file(macroFile);
    if (!file) {
        std::cerr << "Error: Could not open macro file " << macroFile << std::endl;
        return;
    }
    UImanager->ApplyCommand("/control/execute " + macroFile);
}

// Function to set up Geant4 visualization settings
void SetupVisualization(G4UImanager* UImanager) {
    UImanager->ApplyCommand("/vis/open OGL");
    UImanager->ApplyCommand("/vis/viewer/set/upVector 0 1 0");
    UImanager->ApplyCommand("/vis/viewer/set/viewpointVector 0 0 1");
    UImanager->ApplyCommand("/vis/viewer/zoom 1.4");
    UImanager->ApplyCommand("/vis/drawVolume");
    UImanager->ApplyCommand("/tracking/verbose 2");
    UImanager->ApplyCommand("/tracking/storeTrajectory 1");
    UImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate");
    UImanager->ApplyCommand("/vis/scene/add/trajectories");
    UImanager->ApplyCommand("/vis/modeling/trajectories/create/drawByParticleID");
    UImanager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set mu- White");
    UImanager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set e- Yellow");
    UImanager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set G4CMPDriftElectron Violet");
    UImanager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set G4CMPDriftHole Orange");
    UImanager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set phononTS Red");
    UImanager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set phononTF Green");
    UImanager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set phononL Blue");
UImanager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set proton Cyan"); // Using Cyan for protons
    UImanager->ApplyCommand("/vis/viewer/set/style wireframe");
    UImanager->ApplyCommand("/vis/viewer/set/hiddenMarker true");
    UImanager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 0 90");
    UImanager->ApplyCommand("/vis/viewer/zoom 1.6");
}

// Function to configure G4CMP settings
void SetupG4CMP(G4UImanager* UImanager) {
    UImanager->ApplyCommand("/g4cmp/producePhonons 1.0");
    UImanager->ApplyCommand("/g4cmp/sampleLuke 1.0");
    UImanager->ApplyCommand("/g4cmp/produceCharges 1.0");
}

int main(int argc, char** argv) {
    G4RunManager* runManager = new G4RunManager;
    runManager->SetUserInitialization(new RISQTutorialDetectorConstruction());
    
    FTFP_BERT* physics = new FTFP_BERT;
    physics->RegisterPhysics(new G4CMPPhysics);
    physics->SetCuts();
    runManager->SetUserInitialization(physics);
    runManager->SetUserInitialization(new RISQTutorialActionInitialization);
    
    G4CMPConfigManager::Instance();
    RISQTutorialConfigManager::Instance();
    runManager->Initialize();
    
    G4VisManager* visManager = new G4VisExecutive();
    visManager->Initialize();
    
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4UIExecutive* ui = nullptr;
    
    std::string macroFile;
    if (argc > 1) {
        macroFile = argv[1];
        ExecuteMacroFile(UImanager, macroFile);
    } else {
        ui = new G4UIExecutive(argc, argv);
        SetupVisualization(UImanager);
        SetupG4CMP(UImanager);
        ui->SessionStart();
    }
    
    delete ui;
    delete visManager;
    delete runManager;
    return 0;
}
