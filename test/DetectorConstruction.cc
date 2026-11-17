#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"


DetectorConstruction::DetectorConstruction() {}

G4VPhysicalVolume* DetectorConstruction::Construct() {

  std::cout << "DetectorConstruction::Construct" << std::endl;

  //////////////////////////////
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  //////////////////////////////
  // Envelope parameters
  G4double env_sizeXY = 1000*cm, env_sizeZ = 1000*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_Al");

  //////////////////////////////
  // World
  G4double world_sizeXY = env_sizeXY;
  G4double world_sizeZ  = env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");

  auto solidWorld = new G4Box("World",
                              0.5 * world_sizeXY,
                              0.5 * world_sizeXY,
                              0.5 * world_sizeZ);

  auto logicWorld = new G4LogicalVolume(solidWorld,
                                        world_mat,
                                        "World");

  auto physWorld = new G4PVPlacement(nullptr,
                                     G4ThreeVector(),
                                     logicWorld,
                                     "World",
                                     nullptr,
                                     false,
                                     0,
                                     false);


  return physWorld;
}
