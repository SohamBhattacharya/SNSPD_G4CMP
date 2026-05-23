/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef RISQTutorialSteppingAction_hh
#define RISQTutorialSteppingAction_hh 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <set>

class G4Step;

class RISQTutorialSteppingAction : public G4UserSteppingAction
{
public:
  RISQTutorialSteppingAction();
  virtual ~RISQTutorialSteppingAction();
  virtual void UserSteppingAction(const G4Step* step);

private:
  std::set<G4String> fSeenParticles;
};

#endif
