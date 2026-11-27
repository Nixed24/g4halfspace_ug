#ifndef G4HALFPLANE_HH
#define G4HALFPLANE_HH

#include "G4HalfSurface.hh"
#include "G4Types.hh"

class G4HalfPlane : public G4HalfSurface
{
    G4ThreeVector normal;
    G4double distance;

public:

    void set_normal(const G4ThreeVector& n);
    void set_distance(G4double d);
    G4ThreeVector get_normal() const;
    G4double get_distance() const;

    G4HalfPlane(const G4ThreeVector& normal, const G4double& distance);

    G4double ImplicitFunction(const G4ThreeVector& p) const override;

    std::vector<G4ThreeVector> RayIntersection(const G4ThreeVector& p, const G4ThreeVector& v) const override;

    G4double DistanceToSurface(const G4ThreeVector& p) const override;
    G4double DistanceToSurface(const G4ThreeVector& p, const G4ThreeVector& v) const override;

    G4ThreeVector ClosestSurfacePoint(const G4ThreeVector& p) const override;
    G4ThreeVector SurfaceNormal(const G4ThreeVector& /*p*/) const override {return normal;}

    G4HalfSurface* Transform(const G4AffineTransform& t) override;

    void GetExtent(G4ThreeVector& maxBounds, G4ThreeVector& minBounds, const G4ThreeVector& boundLimit) const override;

    G4SurfaceMeshCGAL* GetSurfaceMesh() override;

};
#endif