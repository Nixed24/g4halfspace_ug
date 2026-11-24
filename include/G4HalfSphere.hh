#ifndef G4_HALF_ELLIPSOID
#define G4_HALF_ELLIPSOID

#include "G4HalfQuadric.hh"
#include "G4HalfSurface.hh"

class G4HalfSphere : public G4HalfSurface
{
    using G4HalfSurface::Transform; // Doesn't allow users to end up using the parent Transform

    G4double radius_;
    G4ThreeVector centre_;

    G4HalfQuadric GetQuadric() const;

public:
    G4HalfSphere(const G4double& semiMajor, const G4double& semiMinor, const G4double& focus, const G4ThreeVector& centre);
    ~G4HalfSphere() override = default;

    G4double ImplicitFunction(const G4ThreeVector &p) const override;

    std::vector<G4ThreeVector> RayIntersection(const G4ThreeVector &p, const G4ThreeVector &v) const override;

    G4double DistanceToSurface(const G4ThreeVector &p) const override;
    G4double DistanceToSurface(const G4ThreeVector &p, const G4ThreeVector &v) const override;

    void GetExtent(G4ThreeVector& maxBounds, G4ThreeVector& minBounds, const G4ThreeVector& boundLimit) const override;

    G4ThreeVector SurfaceNormal(const G4ThreeVector &p) const override;
    G4ThreeVector ClosestSurfacePoint(const G4ThreeVector &p) const override;

    G4HalfSurface* Transform(const G4AffineTransform& t) override;
    std::ostream& StreamInfo(std::ostream&) const override;

    G4SurfaceMeshCGAL* GetSurfaceMesh() override;
};

#endif
