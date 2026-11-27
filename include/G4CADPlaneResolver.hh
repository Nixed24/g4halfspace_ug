//
// Created by deck on 11/26/26.
//

#ifndef GEANT4HALFSPACE_G4CADPLANERESOLVER_HH
#define GEANT4HALFSPACE_G4CADPLANERESOLVER_HH

#include "G4HalfPlane.hh"
#include "TopoDS_Shape.hxx"
class G4CADPlaneResolver
{
    struct G4CADPlaneDeck;
    G4HalfPlane MakePlane(const TopoDS_Shape& shape);
};
#endif //GEANT4HALFSPACE_G4CADPLANERESOLVER_HH
