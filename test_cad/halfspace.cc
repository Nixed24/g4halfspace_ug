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

  const G4String cadFilename = argc > 1 ? argv[1] : "cad_files/cube_of_cube_skew.step";
  runManager->SetUserInitialization(new DetectorConstruction(cadFilename));
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
