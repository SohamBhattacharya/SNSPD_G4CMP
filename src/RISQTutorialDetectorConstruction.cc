/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/src/PhononDetectorConstruction.cc \brief
/// Implementation of the PhononDetectorConstruction class
//
// $Id: a2016d29cc7d1e75482bfc623a533d20b60390da $
//
// 20140321  Drop passing placement transform to G4LatticePhysical
// 20211207  Replace G4Logical*Surface with G4CMP-specific versions.
// 20220809  [ For M. Hui ] -- Add frequency dependent surface properties.

#include "RISQTutorialDetectorConstruction.hh"
#include "RISQTutorialSensitivity.hh"
#include "RISQTutorialQubitHousing.hh"
#include "RISQTutorialPad.hh"
#include "RISQTutorialTransmissionLine.hh"
#include "RISQTutorialStraightFluxLine.hh"
#include "RISQTutorialCornerFluxLine.hh"
#include "RISQTutorialResonatorAssembly.hh"
#include "G4CMPPhononElectrode.hh"
#include "G4CMPElectrodeSensitivity.hh"
#include "G4CMPLogicalBorderSurface.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4FieldManager.hh"
#include "G4GeometryManager.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4CMPLogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4Tubs.hh"
#include "G4UniformMagField.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"

using namespace RISQTutorialDetectorParameters;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RISQTutorialDetectorConstruction::RISQTutorialDetectorConstruction()
  : fLiquidHelium(0), fGermanium(0), fAluminum(0), fTungsten(0),
    fWorldPhys(0), 
    fSuperconductorSensitivity(0), fConstructed(false) {;}//, fIfField(true) {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RISQTutorialDetectorConstruction::~RISQTutorialDetectorConstruction() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* RISQTutorialDetectorConstruction::Construct()
{
  if (fConstructed) {
    if (!G4RunManager::IfGeometryHasBeenDestroyed()) {
      // Run manager hasn't cleaned volume stores. This code shouldn't execute
      G4GeometryManager::GetInstance()->OpenGeometry();
      G4PhysicalVolumeStore::GetInstance()->Clean();
      G4LogicalVolumeStore::GetInstance()->Clean();
      G4SolidStore::GetInstance()->Clean();
    }
    // Have to completely remove all lattices to avoid warning on reconstruction
    G4LatticeManager::GetLatticeManager()->Reset();
    // Clear all LogicalSurfaces
    // NOTE: No need to redefine the G4CMPSurfaceProperties
    G4CMPLogicalBorderSurface::CleanSurfaceTable();
  }

  DefineMaterials();
  SetupGeometry();
  fConstructed = true;

  return fWorldPhys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RISQTutorialDetectorConstruction::DefineMaterials()
{ 
  G4NistManager* nistManager = G4NistManager::Instance();

  fLiquidHelium = nistManager->FindOrBuildMaterial("G4_AIR"); // to be corrected
  fGermanium = nistManager->FindOrBuildMaterial("G4_Ge");
  fSilicon = nistManager->FindOrBuildMaterial("G4_Si");
  fAluminum = nistManager->FindOrBuildMaterial("G4_Al");
  fTungsten = nistManager->FindOrBuildMaterial("G4_W");
  fNiobium = nistManager->FindOrBuildMaterial("G4_Nb");
  
  
	//////////////////
	// Materials //
	//////////////////	
 	G4cout<< " ### Starting Material Definition " <<G4endl;    

    G4NistManager *nist = G4NistManager::Instance();

	G4cout<< " ### - Define SiO2" <<G4endl;    
    SiO2 = new G4Material("SiO2", 2.201*g/cm3, 2);
    SiO2->AddElement(nist->FindOrBuildElement("Si"), 1);
    SiO2->AddElement(nist->FindOrBuildElement("O"), 2);
		G4MaterialPropertiesTable *mptSiO2 = new G4MaterialPropertiesTable();
		G4double energySiO2[2] = {1.378*eV, 6.199*eV};
		G4double rindexSiO2[2] = {1.4585, 1.53};
		G4double ABSSiO2[2] = {0.01*mm, 0.01*mm};
		mptSiO2->AddProperty("RINDEX", energySiO2, rindexSiO2, 2);
		//mptSiO2->AddProperty("ABSLENGTH", energySiO2, ABSSiO2, 2);
		SiO2->SetMaterialPropertiesTable(mptSiO2);

	G4cout<< " ### - Define Air" <<G4endl;    
    AirMat = nist->FindOrBuildMaterial("G4_AIR");
		G4double energyWorld[2] = {1.378*eV, 6.199*eV};
		G4double rindexWorld[2] = {1.0, 1.0};
		G4MaterialPropertiesTable *mptWorld = new G4MaterialPropertiesTable();
		mptWorld->AddProperty("RINDEX", energyWorld, rindexWorld, 2);
		AirMat->SetMaterialPropertiesTable(mptWorld);
    
	G4cout << " ### - Define Vacuum" << G4endl;    
	VacuumMat = nist->FindOrBuildMaterial("G4_Galactic");

	// For a vacuum, you may want to define properties like the refractive index (which is typically 1 for vacuum).
	G4double energyVacuum[2] = {1.378*eV, 6.199*eV}; 
	G4double rindexVacuum[2] = {1.0, 1.0};  // Refractive index for vacuum is 1
	G4MaterialPropertiesTable* mptVacuum = new G4MaterialPropertiesTable();
	mptVacuum->AddProperty("RINDEX", energyVacuum, rindexVacuum, 2);

	// Set the material properties table for the vacuum material
	VacuumMat->SetMaterialPropertiesTable(mptVacuum);

    	G4cout << " ### - Define Cu" << G4endl;    
	CuMat = nist->FindOrBuildMaterial("G4_Cu");
	G4MaterialPropertiesTable *mptCu = new G4MaterialPropertiesTable();
	// Valores realistas para el índice de refracción del cobre en el rango de 2-6 eV
	G4double energyCu[2] = {2*eV, 6*eV};
	G4double realRindexCu[2] = {0.271, 0.238};  // Parte real del índice de refracción
	G4double imagRindexCu[2] = {3.56, 4.24};    // Parte imaginaria para la absorción
	mptCu->AddProperty("REALRINDEX", energyCu, realRindexCu, 2);
	mptCu->AddProperty("IMAGINARYRINDEX", energyCu, imagRindexCu, 2);
	CuMat->SetMaterialPropertiesTable(mptCu);
	
	G4cout<< " ### - Define Al" <<G4endl;    
    AlMat = nist->FindOrBuildMaterial("G4_Al");
		G4MaterialPropertiesTable *mptAl = new G4MaterialPropertiesTable();
		G4double energyAl[2] = {400*eV, 1000*eV};
		G4double rindexAl[2] = {0.99, 0.99};
		G4double ABSAl[2] = {6.6e-9*m, 6.63e-9*m};
		mptAl->AddProperty("RINDEX", energyAl, ABSAl, 2);
		mptAl->AddProperty("ABSLENGTH", energyAl, ABSAl, 2);
		AlMat->SetMaterialPropertiesTable(mptAl);
   
  // https://refractiveindex.info/?shelf=main&book=Si3N4&page=Kischkat
  // https://refractiveindex.info/?shelf=main&book=Si3N4&page=Luke
	G4cout << "### - Define Silicon Nitride (Si3N4)" << G4endl;
	// Define the material using the NIST database or creating a custom material
	Si3N4Mat = new G4Material("SiliconNitride", 3.44*g/cm3, 2);
	Si3N4Mat->AddElement(nist->FindOrBuildElement("Si"), 3);
	Si3N4Mat->AddElement(nist->FindOrBuildElement("N"), 4);
	// Define optical properties
	G4MaterialPropertiesTable* mptSi3N4 = new G4MaterialPropertiesTable();
	G4double energySi3N4[2] = {400*eV, 1000*eV}; // Define photon energy range
	//G4double rindexSi3N4[2] = {2.05, 2.05};      // Example refractive index values
	G4double absLengthSi3N4[2] = {4.08*mm, 4.08*mm}; // Absorption length in mm
	// Define photon energies in eV
	const G4int numEntries = 10;
	G4double photonEnergy[numEntries] = {
		4.465*eV, 3.79*eV, 3.115*eV, 2.647*eV, 2.18*eV,
		1.712*eV, 1.245*eV, 0.8294*eV, 0.5178*eV, 0.3619*eV};
	// Corresponding refractive indices
	G4double refractiveIndex[numEntries] = {
		1.8721213170359, 1.9110346794217, 1.943231277958, 1.9619062604624, 1.9778699738446,
		1.9917473810267, 2.0050993473208, 2.0222778887066, 2.0590472271639, 2.1244526847835};
	// Add properties to the material properties table
	mptSi3N4->AddProperty("RINDEX", photonEnergy, refractiveIndex, numEntries);
	mptSi3N4->AddProperty("ABSLENGTH", energySi3N4, absLengthSi3N4, 2);
	// Assign the properties table to the material
	Si3N4Mat->SetMaterialPropertiesTable(mptSi3N4);
	G4cout << "### - Silicon Nitride (Si3N4) defined" << G4endl;

    
	G4cout<< " ### - Define a-Si" <<G4endl;    
    // Define Amorphous Silicon (a-Si)
    G4double density_aSi = 2.32 * g/cm3;  // Typical density for Amorphous Silicon
    G4Element* Si = nist->FindOrBuildElement("Si");  // Silicon element
    aSiMat = new G4Material("AmorphousSi", density_aSi, 1);
    aSiMat->AddElement(Si, 1);  // 1 Silicon atom
		G4MaterialPropertiesTable *mptSi = new G4MaterialPropertiesTable();
		G4double energySi[2] = {1.378*eV, 6.199*eV};
		G4double rindexSi[2] = {4, 4};
		G4double ABSSi[2] = {0.01*mm, 0.01*mm};
		mptSi->AddProperty("RINDEX", energySi, rindexSi, 2);
		mptSi->AddProperty("ABSLENGTH", energySi, ABSSi, 2);
		aSiMat->SetMaterialPropertiesTable(mptSi);

	G4cout<< " ### - Define WSi" <<G4endl;    
    // Define Tungsten Silicide (WSi)
    G4double density_WSi = 9.3 * g/cm3;  // Approximate density of WSi
    G4Element* W = nist->FindOrBuildElement("W");  // Tungsten element
    WSiMat = new G4Material("WSi", density_WSi, 2);  // 2 elements in WSi
    WSiMat->AddElement(W, 1);  // 1 Tungsten atom
    WSiMat->AddElement(Si, 2); // 2 Silicon atoms
    
		G4MaterialPropertiesTable *mptWSi = new G4MaterialPropertiesTable();
		G4double energyWSi[2] = {1.5*eV, 3*eV};
		G4double rindexWSi[2] = {4, 4};
		mptWSi->AddProperty("RINDEX", energyWSi, rindexWSi, 2);
		WSiMat->SetMaterialPropertiesTable(mptWSi);
    
    G4cout<< " ### Finished Material Definition " <<G4endl;    
  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RISQTutorialDetectorConstruction::SetupGeometry()
{
     G4cout << "### SETUP GEOM: " << G4endl;



  //---------------------------------------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------------------------------------
  // First, define border surface properties that can be referenced later
  const G4double GHz = 1e9 * hertz; 

  //the following coefficients and cutoff values are not well-motivated
  //the code below is used only to demonstrate how to set these values.
  const std::vector<G4double> anhCoeffs = {0,0,0,0,0,0};//Turn this off temporarily
  const std::vector<G4double> diffCoeffs = {1,0,0,0,0,0};//Explicitly make this 1 for now
  const std::vector<G4double> specCoeffs = {0,0,0,0,0,0};//Turn this off temporarily
  const G4double anhCutoff = 520., reflCutoff = 350.;   // Units external
    
  
  //These are just the definitions of the interface TYPES, not the interfaces themselves. These must be called in a set of loops
  //below, and invoke these surface definitions.
  if( !fConstructed ){
    fSiNbInterface = new G4CMPSurfaceProperty("SiNbInterface",
					      1.0, 0.0, 0.0, 0.0,
					      0.1, 1.0, 0.0, 0.0);
    fSiCopperInterface = new G4CMPSurfaceProperty("SiCopperInterface",
						  1.0, 0.0, 0.0, 0.0,
						  1.0, 0.0, 0.0, 0.0 );
    fSiVacuumInterface = new G4CMPSurfaceProperty("SiVacuumInterface",
						  0.0, 1.0, 0.0, 0.0,
						  0.0, 1.0, 0.0, 0.0 );
    

    fSiNbInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					    diffCoeffs, specCoeffs, GHz, GHz, GHz);  
    fSiCopperInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
						diffCoeffs, specCoeffs, GHz, GHz, GHz);  
    fSiVacuumInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
						diffCoeffs, specCoeffs, GHz, GHz, GHz);

    //Add a phonon sensor to the interface properties here.
    AttachPhononSensor(fSiNbInterface);
  }
















  //---------------------------------------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------------------------------------
  // Now we start constructing the various components and their interfaces  
  //     
  // World
  //
  G4VSolid* solid_world = new G4Box("World",55.*cm,55.*cm,55.*cm);
  G4LogicalVolume* log_world = new G4LogicalVolume(solid_world,fLiquidHelium,"World");
  //  worldLogical->SetUserLimits(new G4UserLimits(10*mm, DBL_MAX, DBL_MAX, 0, 0));
  log_world->SetVisAttributes(G4VisAttributes::Invisible);
  fWorldPhys = new G4PVPlacement(0,
				 G4ThreeVector(),
				 log_world,
				 "World",
				 0,
                                 false,
				 0);
  
  
  bool checkOverlaps = true;

  



  //-------------------------------------------------------------------------------------------------------------------
  //First, set up the qubit chip substrate. By default, assume that we're using this. Otherwise, it's hard to establish
  //a sensitivity object for this.
  G4Box * solid_siliconChip = new G4Box("QubitChip_solid",
					0.5*dp_siliconChipDimX,
					0.5*dp_siliconChipDimY,
					0.5*dp_siliconChipDimZ);
  
  //Now attribute a physical material to the chip
  G4LogicalVolume * log_siliconChip = new G4LogicalVolume(solid_siliconChip,
							  fSilicon,
							  "SiliconChip_log");
    
  //Now, create a physical volume and G4PVPlacement for storing as the final output
  G4ThreeVector siliconChipTranslate(0,0,0.5*(dp_housingDimZ - dp_siliconChipDimZ) + dp_eps); 
  G4VPhysicalVolume * phys_siliconChip = new G4PVPlacement(0,
							   siliconChipTranslate,
							   log_siliconChip,
							   "SiliconChip", 
							   log_world,
							   false,
							   0,
							   checkOverlaps);

  G4VisAttributes* siliconChipVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.5));
  siliconChipVisAtt->SetVisibility(true);
  log_siliconChip->SetVisAttributes(siliconChipVisAtt);



  //Set up the G4CMP silicon lattice information using the G4LatticeManager
  // G4LatticeManager gives physics processes access to lattices by volume
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  LM->SetVerboseLevel(3);
  G4LatticeLogical* log_siliconLattice = LM->LoadLattice(fSilicon, "Si");
    
  // G4LatticePhysical assigns G4LatticeLogical a physical orientation
  G4LatticePhysical* phys_siliconLattice = new G4LatticePhysical(log_siliconLattice);
  phys_siliconLattice->SetMillerOrientation(1,0,0); 
  LM->RegisterLattice(phys_siliconChip,phys_siliconLattice);

  //Set up border surfaces
  G4CMPLogicalBorderSurface * border_siliconChip_world = new G4CMPLogicalBorderSurface("border_siliconChip_world", phys_siliconChip, fWorldPhys, fSiVacuumInterface);

    



  //-------------------------------------------------------------------------------------------------------------------
  //If desired, set up the copper qubit housing
  if( dp_useQubitHousing ){
      
      
    RISQTutorialQubitHousing * qubitHousing = new RISQTutorialQubitHousing(0,
									   G4ThreeVector(0,0,0),
									   "QubitHousing",
									   log_world,
									   false,
									   0,
									   checkOverlaps);
    G4LogicalVolume * log_qubitHousing = qubitHousing->GetLogicalVolume();
    G4VPhysicalVolume * phys_qubitHousing = qubitHousing->GetPhysicalVolume();
    
    //Set up the logical border surface
    G4CMPLogicalBorderSurface * border_siliconChip_qubitHousing = new G4CMPLogicalBorderSurface("border_siliconChip_qubitHousing", phys_siliconChip, phys_qubitHousing, fSiCopperInterface);
  }
    


    

  //-------------------------------------------------------------------------------------------------------------------
  //Now set up the ground plane, in which the transmission line, resonators, and qubits will be located.
  if( dp_useGroundPlane ){
    
    
    G4Box * solid_groundPlane = new G4Box("GroundPlane_solid",
					  0.5*dp_groundPlaneDimX,
					  0.5*dp_groundPlaneDimY,
					  0.5*dp_groundPlaneDimZ);
    
    
    //Now attribute a physical material to the chip
    G4LogicalVolume * log_groundPlane = new G4LogicalVolume(solid_groundPlane,
							    fNiobium,
							    "GroundPlane_log");
    
    
    //Now, create a physical volume and G4PVPlacement for storing as the final output
    G4ThreeVector groundPlaneTranslate(0,0,0.5*(dp_housingDimZ) + dp_eps + dp_groundPlaneDimZ*0.5);
    G4VPhysicalVolume * phys_groundPlane = new G4PVPlacement(0,
							     groundPlaneTranslate,
							     log_groundPlane,
							     "GroundPlane", 
							     log_world,
							     false,
							     0,
							     checkOverlaps);
    
    G4VisAttributes* groundPlaneVisAtt= new G4VisAttributes(G4Colour(0.0,1.0,1.0,0.5));
    groundPlaneVisAtt->SetVisibility(true);
    log_groundPlane->SetVisAttributes(groundPlaneVisAtt);
    
    
    //Set up the logical border surface
    G4CMPLogicalBorderSurface * border_siliconChip_groundPlane = new G4CMPLogicalBorderSurface("border_siliconChip_groundPlane", phys_siliconChip, phys_groundPlane, fSiNbInterface);


    

    //-------------------------------------------------------------------------------------------------------------------
    //Now set up the transmission line
    if( dp_useTransmissionLine ){
	
      G4ThreeVector transmissionLineTranslate(0,0,0.0);//Since it's within the ground plane exactly; 0.5*(dp_housingDimZ) + dp_eps + dp_groundPlaneDimZ*0.5 ); 
      RISQTutorialTransmissionLine * tLine = new RISQTutorialTransmissionLine(0,
									      transmissionLineTranslate,
									      "TransmissionLine",
									      log_groundPlane,
									      false,
									      0,
									      checkOverlaps);
      G4LogicalVolume * log_tLine = tLine->GetLogicalVolume();
      G4VPhysicalVolume * phys_tLine = tLine->GetPhysicalVolume();



      //Now, if we're using the chip and ground plane AND the transmission line
      //This gets a bit hairy, since the transmission line is composite of both Nb and vacuum.
      //So we'll access the list of physical objects present in it and link those one-by-one to the
      //silicon chip.
      for( int iSubVol = 0; iSubVol < tLine->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol){
	std::cout << "TLine sub volume names (to be used for boundaries): " << std::get<1>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << " with material " << std::get<0>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << std::endl;

	std::string tempName = "border_siliconChip_" + std::get<1>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]);	  
	if( std::get<0>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos ){
	  G4CMPLogicalBorderSurface * border_siliconChip_transmissionLineEmpty = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiVacuumInterface);
	}
	if( std::get<0>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Niobium") != std::string::npos ){
	  G4CMPLogicalBorderSurface * border_siliconChip_transmissionLineConductor = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(tLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiNbInterface);
	}
      }
    }


    //-------------------------------------------------------------------------------------------------------------------
    //Now set up a set of 6 resonator assemblies
    if( dp_useResonatorAssembly ){
      int nR = 6;
      for( int iR = 0; iR < nR; ++iR ){
      
	//First, get the translation vector for the resonator assembly
	//For the top three, don't do a rotation. For the bottom three, do
	G4ThreeVector resonatorAssemblyTranslate(0,0,0);
	G4RotationMatrix * rotAssembly = 0;
	if( iR <= 2 ){
	  resonatorAssemblyTranslate = G4ThreeVector(dp_resonatorLateralSpacing*(iR-1)+dp_centralResonatorOffsetX,
						     0.5 * dp_resonatorAssemblyBaseNbDimY + 0.5 * dp_transmissionLineCavityFullWidth,
						     0.0);
	  rotAssembly = 0;
	}
	else{
	  resonatorAssemblyTranslate = G4ThreeVector(dp_resonatorLateralSpacing*(iR-4)-dp_centralResonatorOffsetX, //Negative offset because qubit is mirrored on underside
						     -1*(0.5 * dp_resonatorAssemblyBaseNbDimY + 0.5 * dp_transmissionLineCavityFullWidth),
						     0.0);
	  rotAssembly = new G4RotationMatrix();
	  rotAssembly->rotateZ(180*deg);
	}
	
	char name[400];
	sprintf(name,"ResonatorAssembly_%d",iR);
	G4String resonatorAssemblyName(name);
	RISQTutorialResonatorAssembly * resonatorAssembly = new RISQTutorialResonatorAssembly(rotAssembly,
											      resonatorAssemblyTranslate,
											      resonatorAssemblyName,
											      log_groundPlane,
											      false,
											      0,
											      checkOverlaps);
	G4LogicalVolume * log_resonatorAssembly = resonatorAssembly->GetLogicalVolume();
	G4VPhysicalVolume * phys_resonatorAssembly = resonatorAssembly->GetPhysicalVolume();
	
	
	//Do the logical border creation now
	for( int iSubVol = 0; iSubVol < resonatorAssembly->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol){
	  std::cout << "TLine sub volume names (to be used for boundaries): " << std::get<1>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]) << " with material " << std::get<0>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]) << std::endl;
	  
	  std::string tempName = "border_siliconChip_" + std::get<1>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]);	  
	  if( std::get<0>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos ){
	    G4CMPLogicalBorderSurface * border_siliconChip_resonatorAssemblyEmpty = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiVacuumInterface);
	  }
	  if( std::get<0>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Niobium") != std::string::npos ){
	    G4CMPLogicalBorderSurface * border_siliconChip_resonatorAssemblyConductor = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(resonatorAssembly->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiNbInterface);
	  }
	}
      }
    }
    
    
    
    //-------------------------------------------------------------------------------------------------------------------
    // Flux lines
    if( dp_useFluxLines ){
      
      
      //--------------------
      G4ThreeVector topStraightFluxLineTranslate(dp_topCenterFluxLineOffsetX,dp_topCenterFluxLineOffsetY,0);
      RISQTutorialStraightFluxLine * topStraightFLine = new RISQTutorialStraightFluxLine(0,
											 topStraightFluxLineTranslate,
											 "TopStraightFluxLine",
											 log_groundPlane,
											 false,
											 0,
											 checkOverlaps);
      G4LogicalVolume * log_topStraightFline = topStraightFLine->GetLogicalVolume();
      G4VPhysicalVolume * phys_topStraightFline = topStraightFLine->GetPhysicalVolume();
      
      //Do the logical border creation now
      for( int iSubVol = 0; iSubVol < topStraightFLine->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol){
	std::cout << "TLine sub volume names (to be used for boundaries): " << std::get<1>(topStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << " with material " << std::get<0>(topStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << std::endl;
	
	std::string tempName = "border_siliconChip_" + std::get<1>(topStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]);	  
	if( std::get<0>(topStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos ){
	  G4CMPLogicalBorderSurface * border_siliconChip_topStraightFluxLineEmpty = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(topStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiVacuumInterface);
	}
	if( std::get<0>(topStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Niobium") != std::string::npos ){
	  G4CMPLogicalBorderSurface * border_siliconChip_topStraightFluxLineConductor = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(topStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiNbInterface);
	}
      }




	


      //--------------------
      G4ThreeVector bottomStraightFluxLineTranslate(dp_topCenterFluxLineOffsetX,-1*dp_topCenterFluxLineOffsetY,0);
      G4RotationMatrix * rotBottomCenter = new G4RotationMatrix();
      rotBottomCenter->rotateZ(180.*deg);
      RISQTutorialStraightFluxLine * bottomStraightFLine = new RISQTutorialStraightFluxLine(rotBottomCenter,
											    bottomStraightFluxLineTranslate,
											    "BottomStraightFluxLine",
											    log_groundPlane,
											    false,
											    0,
											    checkOverlaps);
      G4LogicalVolume * log_bottomStraightFline = bottomStraightFLine->GetLogicalVolume();
      G4VPhysicalVolume * phys_bottomStraightFline = bottomStraightFLine->GetPhysicalVolume();

      //Do the logical border creation now
      for( int iSubVol = 0; iSubVol < bottomStraightFLine->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol){
	std::cout << "TLine sub volume names (to be used for boundaries): " << std::get<1>(bottomStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << " with material " << std::get<0>(bottomStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << std::endl;

	std::string tempName = "border_siliconChip_" + std::get<1>(bottomStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]);	  
	if( std::get<0>(bottomStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos ){
	  G4CMPLogicalBorderSurface * border_siliconChip_bottomStraightFluxLineEmpty = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(bottomStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiVacuumInterface);
	}
	if( std::get<0>(bottomStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Niobium") != std::string::npos ){
	  G4CMPLogicalBorderSurface * border_siliconChip_bottomStraightFluxLineConductor = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(bottomStraightFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiNbInterface);
	}

      }

	

       

      //--------------------
      //Corner flux line 1
      G4ThreeVector topLeftCornerFluxLineTranslate(dp_topLeftFluxLineOffsetX,dp_topLeftFluxLineOffsetY,0);
      G4RotationMatrix * rotTopLeftCenter = new G4RotationMatrix();
      rotTopLeftCenter->rotateZ(0.*deg);
      RISQTutorialCornerFluxLine * topLeftCornerFLine = new RISQTutorialCornerFluxLine(rotTopLeftCenter,
										       topLeftCornerFluxLineTranslate,
										       "GroundPlane_TopLeftCornerFluxLine",
										       log_groundPlane,
										       false,
										       0,
										       checkOverlaps);


	
      G4LogicalVolume * log_topLeftCornerFline = topLeftCornerFLine->GetLogicalVolume();
      G4VPhysicalVolume * phys_topLeftCornerFline = topLeftCornerFLine->GetPhysicalVolume();

	
      //Do the logical border creation now
      for( int iSubVol = 0; iSubVol < topLeftCornerFLine->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol){
	std::cout << "TLine sub volume names (to be used for boundaries): " << std::get<1>(topLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << " with material " << std::get<0>(topLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << std::endl;

	std::string tempName = "border_siliconChip_" + std::get<1>(topLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]);	  
	if( std::get<0>(topLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos ){
	  G4CMPLogicalBorderSurface * border_siliconChip_topLeftCornerFluxLineEmpty = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(topLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiVacuumInterface);
	}
	if( std::get<0>(topLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Niobium") != std::string::npos ){
	  G4CMPLogicalBorderSurface * border_siliconChip_topLeftCornerFluxLineConductor = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(topLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiNbInterface);
	}	  
      }






	

      //--------------------
      //Corner flux line 2
      G4ThreeVector topRightCornerFluxLineTranslate(-1*dp_topLeftFluxLineOffsetX,dp_topLeftFluxLineOffsetY,0);
      G4RotationMatrix * rotTopRightCenter = new G4RotationMatrix();
      rotTopRightCenter->rotateY(180.*deg);
      RISQTutorialCornerFluxLine * topRightCornerFLine = new RISQTutorialCornerFluxLine(rotTopRightCenter,
											topRightCornerFluxLineTranslate,
											"GroundPlane_TopRightCornerFluxLine",
											log_groundPlane,
											false,
											0,
											checkOverlaps);
      G4LogicalVolume * log_topRightCornerFline = topRightCornerFLine->GetLogicalVolume();
      G4VPhysicalVolume * phys_topRightCornerFline = topRightCornerFLine->GetPhysicalVolume();

      //Do the logical border creation now
      for( int iSubVol = 0; iSubVol < topRightCornerFLine->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol){
	std::cout << "TLine sub volume names (to be used for boundaries): " << std::get<1>(topRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << " with material " << std::get<0>(topRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << std::endl;

	std::string tempName = "border_siliconChip_" + std::get<1>(topRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]);	  
	if( std::get<0>(topRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos ){
	  G4CMPLogicalBorderSurface * border_siliconChip_topRightCornerFluxLineEmpty = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(topRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiVacuumInterface);
	}
	if( std::get<0>(topRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Niobium") != std::string::npos ){
	  G4CMPLogicalBorderSurface * border_siliconChip_topRightCornerFluxLineConductor = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(topRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiNbInterface);
	}

      }

    


      //--------------------
      //Corner flux line 3
      G4ThreeVector bottomLeftCornerFluxLineTranslate(dp_topLeftFluxLineOffsetX,-1*dp_topLeftFluxLineOffsetY,0);
      G4RotationMatrix * rotBottomLeftCenter = new G4RotationMatrix();
      rotBottomLeftCenter->rotateX(180.*deg);
      RISQTutorialCornerFluxLine * bottomLeftCornerFLine = new RISQTutorialCornerFluxLine(rotBottomLeftCenter,
											  bottomLeftCornerFluxLineTranslate,
											  "BottomLeftCornerFluxLine",
											  log_groundPlane,
											  false,
											  0,
											  checkOverlaps);
      G4LogicalVolume * log_bottomLeftCornerFline = bottomLeftCornerFLine->GetLogicalVolume();
      G4VPhysicalVolume * phys_bottomLeftCornerFline = bottomLeftCornerFLine->GetPhysicalVolume();

      //Do the logical border creation now
      for( int iSubVol = 0; iSubVol < bottomLeftCornerFLine->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol){
	std::cout << "TLine sub volume names (to be used for boundaries): " << std::get<1>(bottomLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << " with material " << std::get<0>(bottomLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << std::endl;

	  
	std::string tempName = "border_siliconChip_" + std::get<1>(bottomLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]);	  
	if( std::get<0>(bottomLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos ){
	  G4CMPLogicalBorderSurface * border_siliconChip_bottomLeftCornerFluxLineEmpty = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(bottomLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiVacuumInterface);
	}
	if( std::get<0>(bottomLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Niobium") != std::string::npos ){
	  G4CMPLogicalBorderSurface * border_siliconChip_bottomLeftCornerFluxLineConductor = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(bottomLeftCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiNbInterface);
	}

	  
      }

    
    
      //--------------------
      //Corner flux line 4
      G4ThreeVector bottomRightCornerFluxLineTranslate(-1*dp_topLeftFluxLineOffsetX,-1*dp_topLeftFluxLineOffsetY,0);
      G4RotationMatrix * rotBottomRightCenter = new G4RotationMatrix();
      rotBottomRightCenter->rotateX(180.*deg);
      rotBottomRightCenter->rotateY(180.*deg);
      RISQTutorialCornerFluxLine * bottomRightCornerFLine = new RISQTutorialCornerFluxLine(rotBottomRightCenter,
											   bottomRightCornerFluxLineTranslate,
											   "BottomRightCornerFluxLine",
											   log_groundPlane,
											   false,
											   0,
											   checkOverlaps);
      G4LogicalVolume * log_bottomRightCornerFline = bottomRightCornerFLine->GetLogicalVolume();
      G4VPhysicalVolume * phys_bottomRightCornerFline = bottomRightCornerFLine->GetPhysicalVolume();

      //Do the logical border creation now
      for( int iSubVol = 0; iSubVol < bottomRightCornerFLine->GetListOfAllFundamentalSubVolumes().size(); ++iSubVol){
	std::cout << "TLine sub volume names (to be used for boundaries): " << std::get<1>(bottomRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << " with material " << std::get<0>(bottomRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]) << std::endl;

	std::string tempName = "border_siliconChip_" + std::get<1>(bottomRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]);	  
	if( std::get<0>(bottomRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Vacuum") != std::string::npos ){
	  G4CMPLogicalBorderSurface * border_siliconChip_bottomRightCornerFluxLineEmpty = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(bottomRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiVacuumInterface);
	}
	if( std::get<0>(bottomRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]).find("Niobium") != std::string::npos ){
	  G4CMPLogicalBorderSurface * border_siliconChip_bottomRightCornerFluxLineConductor = new G4CMPLogicalBorderSurface(tempName, phys_siliconChip, std::get<2>(bottomRightCornerFLine->GetListOfAllFundamentalSubVolumes()[iSubVol]), fSiNbInterface);
	}
      }
    }
  }




  //---------------------------------------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------------------------------------
  // Now we establish a sensitivity object
  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  if (!fSuperconductorSensitivity)
    fSuperconductorSensitivity = new RISQTutorialSensitivity("PhononElectrode");
  SDman->AddNewDetector(fSuperconductorSensitivity);
  log_siliconChip->SetSensitiveDetector(fSuperconductorSensitivity);



}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Set up a phonon sensor for this surface property object. I'm pretty sure that this
// phonon sensor doesn't get stapled to individual geometrical objects, but rather gets
// stapled to a surface property, but I'm not sure... have to ask mKelsey
void RISQTutorialDetectorConstruction::AttachPhononSensor(G4CMPSurfaceProperty * surfProp)
{
  //If no surface, don't do anything
  if(!surfProp) return;

  //Specify properties of the niobium sensors
  auto sensorProp = surfProp->GetPhononMaterialPropertiesTablePointer();
  sensorProp->AddConstProperty("filmAbsorption",0.0);              //NOT WELL MOTIVATED - probably parametrize and put on slider?
  sensorProp->AddConstProperty("filmThickness",90.*CLHEP::nm);     //Accurate for our thin film.
  sensorProp->AddConstProperty("gapEnergy",1.6e-3*CLHEP::eV);       //Reasonably motivated. Actually, looks like Novotny and Meincke are quoting 2Delta, and this is delta. Nuss and Goossen mention that Nb has a delta value closer to this.
  sensorProp->AddConstProperty("lowQPLimit",3.);                   //NOT WELL MOTIVATED YET -- Dunno how to inform this...
  sensorProp->AddConstProperty("phononLifetime",4.17*CLHEP::ps);   //Kaplan paper says 242ps for Al, same table says 4.17ps for characteristic time for Nb.
  sensorProp->AddConstProperty("phononLifetimeSlope",0.29);        //Based on guessing from Kaplan paper, I think this is material-agnostic?
  sensorProp->AddConstProperty("vSound",3.480*CLHEP::km/CLHEP::s); //True for room temperature, probably good to 10%ish - should follow up
  sensorProp->AddConstProperty("subgapAbsorption",0.0);            //Assuming that since we're mostly sensitive to quasiparticle density, phonon "heat" here isn't something that we're sensitive to? Unsure how to select this.

  //  sensorProp->AddConstProperty("gapEnergy",3.0e-3*CLHEP::eV);      //Reasonably motivated. Novotny and Meincke, 1975 (2.8-3.14 meV)
  //  sensorProp->AddConstProperty("phononLifetime",242.*ps);      //Kaplan paper says 242ps for Al, same table says 4.17ps for characteristic time for Nb.
  
  surfProp->SetPhononElectrode(new G4CMPPhononElectrode);
  
}
