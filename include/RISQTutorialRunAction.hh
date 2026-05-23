/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef RISQTutorialRunAction_hh
#define RISQTutorialRunAction_hh 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class RISQTutorialRunAction : public G4UserRunAction
{
public:
  RISQTutorialRunAction();
  virtual ~RISQTutorialRunAction();
  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);
};

#endif
