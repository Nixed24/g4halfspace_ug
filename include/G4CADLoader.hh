#ifndef GEANT4_G4CADLOADER_HH
#define GEANT4_G4CADLOADER_HH

#include "G4String.hh"
#include "TopoDS_Shape.hxx"

class G4CADLoader
{
public:
  TopoDS_Shape Load(const G4String& filename) const;
};

#endif
