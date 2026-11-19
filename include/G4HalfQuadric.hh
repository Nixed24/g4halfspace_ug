#ifndef G4HALF_QUADRIC_HH
#define G4HALF_QUADRIC_HH

#include <vector>
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/Vector.h>

#include "G4HalfSurface.hh"
#include "G4Types.hh"
#include "G4VSolid.hh"

class G4HalfQuadric : public G4HalfSurface {
  private:
    static G4RotationMatrix _transpose(const G4RotationMatrix & M);

    std::vector<G4ThreeVector> GetLMContourIntersections(const G4ThreeVector& p) const;

    CLHEP::HepMatrix Q {CLHEP::HepMatrix(3,3,1)};
    CLHEP::HepVector P {CLHEP::HepVector(3,0)};
    G4double R {-1}; // subsequently R = -1 gives a unit sphere at the origin

    enum QuadricType
    { // Excluding imaginary types
      kEllipsoid,
      kHyperboloidOneSheet,
      kHyperboloidTwoSheets,
      kQuadricCone,
      kEllipticParaboloid,
      kHyperbolicParaboloid,
      kEllipticCylinder,
      kIntersectingPlanes,
      kParabolicCylinder,
      kCoincidentPlanes
    };

  public:

    CLHEP::HepMatrix get_Q() const;
    CLHEP::HepVector get_P() const;
    G4double get_R() const;
    void set_Q(const CLHEP::HepMatrix& new_Q);
    void set_P(const CLHEP::HepVector& new_P);
    void set_R(G4double new_R);

    G4HalfQuadric() = default;
    G4HalfQuadric(G4double R);
    G4HalfQuadric(const CLHEP::HepMatrix& Q, const CLHEP::HepVector& P, G4double R);
    G4HalfQuadric(const G4HalfQuadric& other) = default;
    ~G4HalfQuadric() override = default;

    G4double ImplicitFunction(const G4ThreeVector& p) const override;

    std::vector<G4ThreeVector> RayIntersection(const G4ThreeVector& p, const G4ThreeVector& v) const override;

    G4double DistanceToSurface(const G4ThreeVector& p, const G4ThreeVector& v) const override;
    G4double DistanceToSurface(const G4ThreeVector&p) const override;

    G4RotationMatrix GetRotationMatrix(G4ThreeVector &QEigenvalues) const;
    QuadricType IdentifyQuadricType() const;
    void GetQuadricBound(G4ThreeVector& maxBounds, G4ThreeVector& minBounds, const G4ThreeVector& boundLimit) const;

    // This will involve taking a gradient of the implicit function and evaluating it at p
    G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const override;
    G4ThreeVector ClosestSurfacePoint(const G4ThreeVector& p) const override;

    void Transform(const G4AffineTransform& t) override;

    std::ostream& StreamInfo(std::ostream&) const override;

    G4SurfaceMeshCGAL* GetSurfaceMesh() override;

};

#endif