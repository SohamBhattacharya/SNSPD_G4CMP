/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef RISQTutorialPrimaryGeneratorAction_h
#define RISQTutorialPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include "globals.hh"

class G4Event;

class RISQTutorialPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  RISQTutorialPrimaryGeneratorAction();
  virtual ~RISQTutorialPrimaryGeneratorAction();
  virtual void GeneratePrimaries(G4Event*);

private:
  G4GeneralParticleSource* fParticleGun;
};

#endif
