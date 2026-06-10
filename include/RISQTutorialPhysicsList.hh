/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file RISQTutorialPhysicsList.hh
/// \brief Maximum-fidelity physics list for RISQTutorial.
///
/// Combines:
///   - Penelope EM (most accurate below 1 keV in Si/SiO2)
///   - Full atomic de-excitation: fluorescence, Auger cascade, PIXE
///   - DNA physics on the SiO2 region (track-structure level, ~1 eV)
///   - FTFP_BERT hadronic for 120 GeV proton
///   - High-precision neutron (HP) for any thermal neutrons
///   - Optical photon physics (scintillation, Cherenkov, Rayleigh)
///   - Radioactive decay
///   - G4CMP phonon + charge carrier transport

#ifndef RISQTutorialPhysicsList_hh
#define RISQTutorialPhysicsList_hh 1

#include "G4VModularPhysicsList.hh"

class RISQTutorialPhysicsList : public G4VModularPhysicsList
{
public:
  RISQTutorialPhysicsList();
  virtual ~RISQTutorialPhysicsList();

  virtual void SetCuts()            override;
  virtual void ConstructParticle()  override;
  virtual void ConstructProcess()   override;
};

#endif
