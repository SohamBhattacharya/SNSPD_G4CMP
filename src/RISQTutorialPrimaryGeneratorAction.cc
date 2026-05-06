/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/src/RISQTutorialPrimaryGeneratorAction.cc
/// \brief Implementation of the RISQTutorialPrimaryGeneratorAction class
//
// $Id: e75f788b103aef810361fad30f75077829192c13 $
//
// 20140519  Allow the user to specify phonon type by name in macro; if
//	     "geantino" is set, use random generator to select.

#include "RISQTutorialPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4Geantino.hh"
#include "G4ParticleGun.hh"
//#include "G4GeneralParticleSource.hh"
#include "G4RandomDirection.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhononLong.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

RISQTutorialPrimaryGeneratorAction::RISQTutorialPrimaryGeneratorAction() { 
  //fParticleGun  = new G4GeneralParticleSource();
  
  

  // default particle kinematics ("geantino" triggers random phonon choice)
  //  fParticleGun->SetParticleDefinition(G4Geantino::Definition());
  //  fParticleGun->SetParticleMomentumDirection(G4RandomDirection());
  //  fParticleGun->SetParticlePosition(G4ThreeVector(0.0,0.0,0.4998095*CLHEP::cm));
  //  fParticleGun->SetParticleEnergy(0.0075*eV);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


RISQTutorialPrimaryGeneratorAction::~RISQTutorialPrimaryGeneratorAction() {
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

 
void RISQTutorialPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

  /*
  if (fParticleGun->GetParticleDefinition() == G4Geantino::Definition()) {
    G4double selector = G4UniformRand();
    if (selector<0.53539) {
      fParticleGun->SetParticleDefinition(G4PhononTransSlow::Definition()); 
    } else if (selector<0.90217) {
      fParticleGun->SetParticleDefinition(G4PhononTransFast::Definition());
    } else {
      fParticleGun->SetParticleDefinition(G4PhononLong::Definition());
    }
  }
  */
  
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
  //  fParticleGun->SetParticleMomentumDirection(G4RandomDirection());
  //fParticleGun->GeneratePrimaryVertex(anEvent);
  
  
  	G4cout<< " ### Starting Generator  " <<G4endl;    

    fParticleGun = new G4ParticleGun(1); /*Number of particles*/	
    
    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();

    // Define the proton particle
    G4String particleName = "proton";
    //G4String particleName = "mu+";
    G4ParticleDefinition *particle_p = particleTable->FindParticle(particleName);
    
    // Set particle properties for 120 GeV proton
    fParticleGun->SetParticleDefinition(particle_p);
    fParticleGun->SetParticleMomentum(120. * GeV); // Set momentum to 120 GeV    
    
    //G4String particleName = "e-";
    //G4ParticleDefinition *particle_e = particleTable->FindParticle(particleName);
    
    // Set particle properties for 120 GeV proton
    //fParticleGun->SetParticleDefinition(particle_e);
    //fParticleGun->SetParticleMomentum(8. * GeV); // Set momentum to 120 GeV    
    
	// Declare pos outside the if-else blocks
	G4ThreeVector pos;
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1)); // Direction -Z

	// Check if randomGunLocation is true or false

		// Set position to -200 mm in Z with X and Y as 0 if randomGunLocation is false
		pos = G4ThreeVector(0. * mm, 0. * mm, -200 * mm);
	
	//PassArgs->StorePosition(pos);

	// Set the particle gun position
	fParticleGun->SetParticlePosition(pos);

	G4cout<< " ### Finshing Generator  " <<G4endl;    
    
    fParticleGun->GeneratePrimaryVertex(anEvent);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


