#pragma once

#include "G4VUserDetectorConstruction.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

private:
    G4VPhysicalVolume* fWorld;
};