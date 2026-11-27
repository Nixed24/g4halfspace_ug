//
// Created by deck on 11/26/26.
//

#include "G4CADPlaneResolver.hh"

struct G4CADPlaneResolver::G4CADPlaneDeck
{
    std::vector<std::pair<G4HalfPlane, G4ThreeVector>> planes_and_points;
    void Tidy()
    {
        if (planes_and_points.size() < 3) return;
        //G4ThreeVector current_point = planes_and_points.begin();
        for (auto& pair : planes_and_points)
        {
            if (pair.first.Inside(current_point) == kOutside) // implicit function with an epsilon would also work
            {
                pair.first.set_normal(-pair.first.get_normal());
            }
        }
    }
};