#include "DetectorConstruction.hh"

#include "G4CADLoader.hh"
#include "G4CADTopology.hh"

#include "G4Box.hh"
#include "G4ios.hh"
#include "G4NistManager.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

DetectorConstruction::DetectorConstruction(const G4String& cadFilename)
    : fCADFilename(cadFilename)
{}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  const TopoDS_Shape cadShape = G4CADLoader().Load(fCADFilename);
  std::size_t topologyNodeCount = 0;
  G4CADTopology(cadShape).Traverse(
      [&topologyNodeCount](const TopoDS_Shape&, std::size_t) {
        ++topologyNodeCount;
      });
  G4cout << "Loaded CAD topology with " << topologyNodeCount << " nodes"
         << G4endl;

  G4NistManager* nist = G4NistManager::Instance();
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");

  auto* worldSolid = new G4Box("World", 5.0 * m, 5.0 * m, 5.0 * m);
  auto* worldLogical = new G4LogicalVolume(worldSolid, world_mat, "World");
  auto* worldPhysical = new G4PVPlacement(
      nullptr, G4ThreeVector(), worldLogical, "World", nullptr, false, 0);

  // CAD-derived solids will be placed here once the geometry conversion
  // classes are available.
  auto* placeholderSolid =
      new G4Box("CADPlaceholder", 0.5 * m, 0.5 * m, 0.5 * m);
  auto* placeholderLogical = new G4LogicalVolume(
      placeholderSolid, nist->FindOrBuildMaterial("G4_Al"), "CADPlaceholder");
  new G4PVPlacement(nullptr, G4ThreeVector(), placeholderLogical,
                    "CADPlaceholder", worldLogical, false, 0);

  return worldPhysical;
}

