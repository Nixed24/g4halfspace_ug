#ifndef G4_HALF_HYPERBOLIC_PARABOLOID
#define G4_HALF_HYPERBOLIC_PARABOLOID

#include "G4HalfQuadric.hh"
#include "G4HalfSurface.hh"

class G4HalfHyperbolicParaboloid : public G4HalfSurface
{
    using G4HalfSurface::Transform; // Doesn't allow users to end up using the parent Transform

    G4double a_;
    G4double b_;
    G4double c_; // Sets scale for z

    std::vector<G4double> GetLambdas(const G4ThreeVector& p) const;
    std::vector<G4ThreeVector> GetLMIntersections(const G4ThreeVector& p) const;
    G4HalfQuadric GetQuadric() const;

public:
    G4HalfHyperbolicParaboloid(const G4double& a, const G4double& b, const G4double& c);
    ~G4HalfHyperbolicParaboloid() override = default;

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