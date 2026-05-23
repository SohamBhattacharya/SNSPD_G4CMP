/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "RISQTutorialActionInitialization.hh"
#include "RISQTutorialPrimaryGeneratorAction.hh"
#include "RISQTutorialRunAction.hh"
#include "RISQTutorialSteppingAction.hh"
#include "G4CMPStackingAction.hh"

void RISQTutorialActionInitialization::Build() const
{
  // RunAction MUST be registered first — it creates the ntuple
  // that SteppingAction fills. Wrong order → crash on first step.
  SetUserAction(new RISQTutorialRunAction);
  SetUserAction(new RISQTutorialPrimaryGeneratorAction);
  SetUserAction(new G4CMPStackingAction);
  SetUserAction(new RISQTutorialSteppingAction);
}
