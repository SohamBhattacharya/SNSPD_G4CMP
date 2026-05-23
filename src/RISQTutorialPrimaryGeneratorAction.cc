/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/src/RISQTutorialPrimaryGeneratorAction.cc
/// \brief Implementation of the RISQTutorialPrimaryGeneratorAction class
//
// Restored to use G4GeneralParticleSource so macros control the particle.

#include "RISQTutorialPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"

RISQTutorialPrimaryGeneratorAction::RISQTutorialPrimaryGeneratorAction()
{
  fParticleGun = new G4GeneralParticleSource();
}

RISQTutorialPrimaryGeneratorAction::~RISQTutorialPrimaryGeneratorAction()
{
  delete fParticleGun;
}

void RISQTutorialPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
