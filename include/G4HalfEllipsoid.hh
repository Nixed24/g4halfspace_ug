#ifndef G4_HALF_ELLIPSOID
#define G4_HALF_ELLIPSOID

#include "G4HalfQuadric.hh"
#include "G4HalfSurface.hh"

class G4HalfEllipsoid : public G4HalfSurface
{
    using G4HalfSurface::Transform; // Doesn't allow users to end up using the parent Transform

    G4double semiMajor_;
    G4double semiMinor_;
    G4double focus_;
    G4ThreeVector centre_;

    std::vector<G4double> GetLambdas(const G4ThreeVector& p) const;

    G4HalfQuadric GetQuadric() const;

public:
    G4HalfEllipsoid(const G4double& semiMajor, const G4double& semiMinor, const G4double& focus, const G4ThreeVector& centre);
    ~G4HalfEllipsoid() override = default;

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