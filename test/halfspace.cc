#include "G4HalfQuadric.hh"

#include "G4BooleanSolid.hh"
#include "G4GDMLParser.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4RunManagerFactory.hh"
#include "G4TransportationManager.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4ios.hh"
#include "FTFP_BERT.hh"
#include "G4SystemOfUnits.hh"

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

int main(int argc, char** argv)
{

  auto* runManager = G4RunManagerFactory::CreateRunManager();

    auto cyllinder = G4HalfQuadric();
    cyllinder.set_R(-3);
    G4RotationMatrix rotation;
    rotation.rotateZ(30*deg);

    G4AffineTransform rotate(rotation);
    sphere.Transform(rotate);
    G4ThreeVector eigenvals;

    G4cout << "SPHERE TEST" << G4endl;
    sphere.GetRotationMatrix(eigenvals);

  runManager->SetUserInitialization(new DetectorConstruction);
  runManager->SetUserInitialization(new ActionInitialization);
  runManager->SetUserInitialization(new FTFP_BERT);

  runManager->Initialize();

  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  runManager->BeamOn(0);

  if (argc == 5)  // batch mode
  {
      G4String command = "/control/execute ";
      G4String fileName = argv[4];
      UImanager->ApplyCommand(command + fileName);
  }
  else  // interactive mode
  {
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
      UImanager->ApplyCommand("/control/execute vis.mac");
      ui->SessionStart();
      delete ui;
  }

  delete visManager;
  delete runManager;

  return 0;
}