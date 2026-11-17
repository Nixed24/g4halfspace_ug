#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4HalfIntersectionSolid.hh"
#include "G4HalfUnionSolid.hh"
#include "G4HalfPlane.hh"


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

    G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
    const G4double d = 0.5 * m;

    G4ThreeVector shapePosition;

    auto makeCube = [d](G4HalfIntersectionSolid& s) {
        s.AddIntersectionSurface<G4HalfPlane>(G4ThreeVector(1. ,0., 0.), d);
        s.AddIntersectionSurface<G4HalfPlane>(G4ThreeVector(0., 1., 0.), d);
        s.AddIntersectionSurface<G4HalfPlane>(G4ThreeVector(0., 0., 1.), d);
        s.AddIntersectionSurface<G4HalfPlane>(G4ThreeVector(-1., 0., 0.), d);
        s.AddIntersectionSurface<G4HalfPlane>(G4ThreeVector(0., -1., 0.), d);
        s.AddIntersectionSurface<G4HalfPlane>(G4ThreeVector(0., 0., -1.), d);
    };

    auto* solidShape = new G4HalfUnionSolid("Shape");

    auto& outerCube = solidShape->AddIntersectionSolid<G4HalfIntersectionSolid>("OuterCube");
    makeCube(outerCube);

    auto* logicShape = new G4LogicalVolume(solidShape,
                                            shape2_mat,
                                            "Shape");
    new G4PVPlacement(nullptr,  // no rotation
                    shapePosition,  // at position
                    logicShape,  // its logical volume
                    "Shape",  // its name
                    logicWorld,  // its mother  volume
                    false,  // no boolean operation
                    0,  // copy number
                    false);
  return physWorld;
}
