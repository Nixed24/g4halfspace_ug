#include "G4HalfSphere.hh"
#include "G4Orb.hh"

G4HalfSphere::G4HalfSphere(const G4double& radius, const G4ThreeVector& centre) : radius_(radius), centre_(centre)
{
    if (radius <= 0.0)
        throw std::invalid_argument("G4HalfSphere: radius must be positive");
}

G4double G4HalfSphere::ImplicitFunction(const G4ThreeVector &p) const
{
    const G4ThreeVector q = p - centre_;
    return q[0]*q[0] + q[1]*q[1] + q[2]*q[2] - radius_*radius_;
}

std::vector<G4ThreeVector> G4HalfSphere::RayIntersection(const G4ThreeVector &p, const G4ThreeVector &v) const
{
    if (v.mag() == 0)
    {
        throw std::invalid_argument("G4HalfSphere::RayIntersection: v must be non-zero.");
    }

    const G4ThreeVector q = p - centre_;
    const G4ThreeVector vUnit = v.unit();

    const G4double a = vUnit.dot(vUnit);
    const G4double b = 2*q.dot(vUnit);
    const G4double c = q.dot(q) - radius_*radius_;

    const G4double determinant = b*b - 4*a*c;

    const G4double minIntersectionDistance = (-b - sqrt(determinant))/(2*a);
    const G4double maxIntersectionDistance = (-b + sqrt(determinant))/(2*a);

    return {p + minIntersectionDistance*vUnit, p + maxIntersectionDistance*vUnit};
}


G4double G4HalfSphere::DistanceToSurface(const G4ThreeVector &p, const G4ThreeVector &v) const
{
    if (v.mag() == 0)
        throw std::logic_error("G4HalfEllipsoid::DistanceToSurface: Cannot pass 0 length direction vector.");

    return (p - RayIntersection(p, v)[0]).mag();
}

G4double G4HalfSphere::DistanceToSurface(const G4ThreeVector &p) const
{
    const G4ThreeVector q = p - centre_;

    return sqrt(q.dot(q)) - radius_;
}

void G4HalfSphere::GetExtent(G4ThreeVector &maxBounds, G4ThreeVector &minBounds, const G4ThreeVector &boundLimit) const
{
    maxBounds =
    {
        std::min(maxBounds[0], radius_),
        std::min(maxBounds[1], radius_),
        std::min(maxBounds[2], radius_)
    };
    minBounds =
    {
        std::max(-maxBounds[0], -radius_),
        std::max(-maxBounds[1], -radius_),
        std::min(-maxBounds[2], -radius_)
    };
}

G4ThreeVector G4HalfSphere::SurfaceNormal(const G4ThreeVector &p) const
{
    return G4ThreeVector(2*p[0], 2*p[1], 2*p[2]).unit();
}

G4ThreeVector G4HalfSphere::ClosestSurfacePoint(const G4ThreeVector &p) const
{
    const G4ThreeVector v = p - centre_;
    return RayIntersection(p, v)[0];
}

G4HalfSurface* G4HalfSphere::Transform(const G4AffineTransform& t)
{
    const G4ThreeVector translate = t.NetTranslation();

    const auto sphere = new G4HalfSphere(this->radius_, this->centre_ + translate);

    return sphere;
}

std::ostream &G4HalfSphere::StreamInfo(std::ostream &os) const
{
    os << ">>> G4HalfSphere <<<" << G4endl;
    os << "Radius = " << this->radius_ << G4endl;
    os << "Centre = " << this->centre_ << G4endl;

    return os;
}

G4SurfaceMeshCGAL *G4HalfSphere::GetSurfaceMesh()
{
    const auto o = G4Orb("rendering sphere", radius_);
    G4Polyhedron *g4poly = o.GetPolyhedron();
    auto *sm = new G4SurfaceMeshCGAL();
    sm->Fill(g4poly);
    sm->Translate(centre_);
    return sm;
}
