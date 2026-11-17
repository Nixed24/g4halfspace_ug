#include "G4HalfUnionSolid.hh"

#include "G4VGraphicsScene.hh"

G4HalfUnionSolid::G4HalfUnionSolid(const G4String& name) : G4HalfSolid(name)
{
}

G4HalfUnionSolid::G4HalfUnionSolid(const std::vector<G4HalfSolid*>& unionSolids,
                                   const G4String& name)
                                      :  G4HalfSolid(name),
                                         unionSolids_(unionSolids)

{
}

G4HalfUnionSolid::~G4HalfUnionSolid()
{
    for (const auto s : unionSolids_)
        if (s) delete s;
}

void G4HalfUnionSolid::AddUnionSolid(G4HalfSolid *unionSolid)
{
    unionSolids_.push_back(unionSolid);
}

G4double G4HalfUnionSolid::ImplicitFunction(const G4ThreeVector& p) const
{
    G4double implicitFunction = kInfinity;

    for (const auto &s : unionSolids_)
    {
        implicitFunction = std::min(implicitFunction, s->ImplicitFunction(p));
    }

    return implicitFunction;
}

void G4HalfUnionSolid::GetPotentialRayIntersections(const G4ThreeVector& p, const G4ThreeVector& v,
    std::vector<RayHit>& hits, const G4bool keepOnlySolidHits) const
{
    if (v.mag() == 0) throw std::logic_error("Cannot pass zero magnitude unit vector");

    const G4ThreeVector vUnit = v.unit();

    for (const auto& solid : unionSolids_) {
        solid->GetPotentialRayIntersections(p, vUnit, hits, keepOnlySolidHits);
    }


    std::sort(hits.begin(), hits.end(),
            [](const RayHit& a, const RayHit& b)
            {
                return a.distance < b.distance;
            });
}

void G4HalfUnionSolid::GetOnSurfaceNormal(const G4ThreeVector& p, G4ThreeVector& norm) const
{
    const G4bool onSolid = this->Inside(p) == kSurface;

    if (not onSolid)
        throw std::logic_error("G4HalfIntersectionSolid::GetOnSurfaceNormal: Point passed in is not on the surface.");

    for (const auto &solid : unionSolids_)
    {
        if (solid->Inside(p) == kSurface)
            norm = solid->SurfaceNormal(p);
    }
}

G4double G4HalfUnionSolid::SignedDistance(const G4ThreeVector& p, G4ThreeVector* n, const G4bool calcNorm) const
{
    G4double signedDistance = kInfinity;

    for (const auto &solid : unionSolids_)
    {
        G4ThreeVector candidateNorm;
        G4double candidateSignedDistance = solid->SignedDistance(p, &candidateNorm, calcNorm);

        if (candidateSignedDistance < signedDistance)
        {
            signedDistance = candidateSignedDistance;
            if (calcNorm && n != nullptr)
            {
                n = &candidateNorm;
            }
        }
    }

    return signedDistance;
}

G4double G4HalfUnionSolid::DistanceToIn(const G4ThreeVector& p, const G4bool calcNorm, G4ThreeVector* n) const
{
    return std::max(0.0, this->SignedDistance(p, n, calcNorm));
}


G4double G4HalfUnionSolid::DistanceToIn(const G4ThreeVector& p) const
{
    return DistanceToIn(p, false, nullptr);
}


G4double G4HalfUnionSolid::DistanceToOut(const G4ThreeVector& p, const G4bool calcNorm, G4ThreeVector* n) const
{
    return std::max(0.0, -this->SignedDistance(p, n, calcNorm));
}


G4double G4HalfUnionSolid::DistanceToOut(const G4ThreeVector& p) const
{
    return DistanceToOut(p, false, nullptr);
}

void G4HalfUnionSolid::Transform(const G4AffineTransform& t) const
{
    for (auto &solid : unionSolids_) solid->Transform(t);
}


G4GeometryType G4HalfUnionSolid::GetEntityType() const
{
    return "G4HalfUnionSolid";
}

std::ostream& G4HalfUnionSolid::StreamInfo(std::ostream& os) const
{
    os << "--- Half Union Solid: " << GetName() << " ---" << std::endl;

    return os;
}

G4SurfaceMeshCGAL* G4HalfUnionSolid::GetSurfaceMesh() {

    G4cout << "G4HalfUnionSolid::GetSurfaceMesh" << G4endl;

    G4SurfaceMeshCGAL *sm1 = new G4SurfaceMeshCGAL();

    for(const auto &s : unionSolids_)
    {
        auto sm2 = s->GetSurfaceMesh();
        G4bool valid;
        sm1 = sm1->Union(sm2, valid);
    }

    return sm1;

}

void G4HalfUnionSolid::DescribeYourselfTo(G4VGraphicsScene& scene) const {
    G4cout << "G4HalfUnionSolid::DescribeYourselfTo" << std::endl;

    auto sm = const_cast<G4HalfUnionSolid*>(this)->GetSurfaceMesh();
    auto ph = sm->GetG4Polyhedron();
    scene.AddPrimitive(*ph);
}
