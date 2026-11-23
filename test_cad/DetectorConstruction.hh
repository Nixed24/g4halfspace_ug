#pragma once

#include "G4String.hh"
#include "G4VUserDetectorConstruction.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    explicit DetectorConstruction(const G4String& cadFilename = "model.step");

    G4VPhysicalVolume* Construct() override;

private:
    G4String fCADFilename;
};
