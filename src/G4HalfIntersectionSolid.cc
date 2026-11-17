#include "G4HalfIntersectionSolid.hh"
#include <algorithm>

#include "G4HalfPlane.hh"
#include "G4VGraphicsScene.hh"
#include "G4VoxelLimits.hh"

const G4ThreeVector ORIGIN = {0,0,0};

G4HalfIntersectionSolid::G4HalfIntersectionSolid(const G4String& name) : G4HalfSolid(name)
{
}

G4HalfIntersectionSolid::G4HalfIntersectionSolid(const std::vector<G4HalfSurface*>& intersectionSurfaces,
                                                 const std::vector<G4HalfSurface*>& subtractionSurfaces,
                                                 const G4String& name)
                                                     :  G4HalfSolid(name),
                                                        intersectionSurfaces_(intersectionSurfaces),
                                                        subtractionSurfaces_(subtractionSurfaces)

{
}

G4HalfIntersectionSolid::~G4HalfIntersectionSolid()
{
    for (const auto s : intersectionSurfaces_)
        if (s) delete s;


    for (const auto s : subtractionSurfaces_)
        if (s) delete s;

}

void G4HalfIntersectionSolid::AddIntersectionSurface(G4HalfSurface *surface)
{
    intersectionSurfaces_.push_back(surface);
}

void G4HalfIntersectionSolid::AddSubtractionSurface(G4HalfSurface *surface)
{
    subtractionSurfaces_.push_back(surface);
}

G4double G4HalfIntersectionSolid::ImplicitFunction(const G4ThreeVector &p) const
{
    G4double implicitFunction = -kInfinity;

    for (const auto &s : intersectionSurfaces_)
    {
        implicitFunction = std::max(implicitFunction, s->ImplicitFunction(p));
    }

    for (const auto &s : subtractionSurfaces_)
    {
        implicitFunction = std::max(implicitFunction, -s->ImplicitFunction(p));
    }

    return implicitFunction;
}


void G4HalfIntersectionSolid::GetPotentialRayIntersections(const G4ThreeVector& p, const G4ThreeVector& v,
    std::vector<RayHit>& hits, const G4bool keepOnlySolidHits) const
{
    if (v.mag() == 0) throw std::logic_error("Cannot pass zero magnitude directional vector.");

    const G4ThreeVector vUnit = v.unit();

    auto getHits = [&](const  G4HalfSurface* surface, const G4bool subtraction)
    {
        for (const G4ThreeVector& intersection : surface->RayIntersection(p,vUnit))
        {
            if (keepOnlySolidHits && this->Inside(intersection) != kSurface) continue;

            RayHit hit;
            hit.intersection = intersection;
            hit.distance = (intersection - p).dot(vUnit);
            hit.targetSurface = surface;
            hit.isSubtraction = subtraction;

            hits.push_back(hit);
        }
    };

    for (const auto &s : intersectionSurfaces_)
        getHits(s, false);


    for (const auto &s : subtractionSurfaces_)
        getHits(s, true);

    std::sort(hits.begin(), hits.end(),
        [](const RayHit& a, const RayHit& b)
        {
            return a.distance < b.distance;
        });
}



G4double G4HalfIntersectionSolid::SignedDistance(const G4ThreeVector& p, G4ThreeVector *n, const G4bool calcNorm) const
{
    G4double signedDistance = -kInfinity;
    const G4HalfSurface* candidateSurface = nullptr;
    G4int normSign;

    auto setDistance = [&](const G4HalfSurface* surface, const G4int sign)
    {
        G4double candidateDistance = surface->DistanceToSurface(p);
        if (sign*candidateDistance > signedDistance)
        {
            signedDistance = candidateDistance;
            candidateSurface = surface;
            normSign = sign;
        }
    };

    for (const auto &surface : intersectionSurfaces_)
    {
        setDistance(surface, +1);
    }

    for (const auto &surface : subtractionSurfaces_)
    {
        setDistance(surface, -1);
    }

    if (calcNorm &&  n != nullptr)
    {
        const G4ThreeVector surfacePoint = candidateSurface->ClosestSurfacePoint(p);
        *n = normSign*candidateSurface->SurfaceNormal(surfacePoint);
    }

    return signedDistance;
}

G4double G4HalfIntersectionSolid::DistanceToIn(const G4ThreeVector& p, const G4bool calcNorm, G4ThreeVector* n) const
{
    return std::max(0.0, SignedDistance(p, n, calcNorm));
}

G4double G4HalfIntersectionSolid::DistanceToIn(const G4ThreeVector& p) const
{
    return DistanceToIn(p, false, nullptr);
}


G4double G4HalfIntersectionSolid::DistanceToOut(const G4ThreeVector& p, const G4bool calcNorm, G4ThreeVector* n) const
{
    return std::max(0.0, -SignedDistance(p, n, calcNorm));
}


G4double G4HalfIntersectionSolid::DistanceToOut(const G4ThreeVector& p) const
{
    return DistanceToOut(p, false, nullptr);
}


void G4HalfIntersectionSolid::GetOnSurfaceNormal(const G4ThreeVector& p, G4ThreeVector& norm) const
{
    const G4bool onSolid = this->Inside(p) == kSurface;
    if (not onSolid)
        throw std::logic_error("G4HalfIntersectionSolid::GetOnSurfaceNormal: Point passed in is not on the surface.");

    for (const auto &s : intersectionSurfaces_)
    {
        if (s->Inside(p) == kSurface)
        {
            norm = (onSolid) ? s->SurfaceNormal(p) : s->SurfaceNormal(s->ClosestSurfacePoint(p));
        }
    }
    for (const auto &s : subtractionSurfaces_)
    {
        if (s->Inside(p) == kSurface)
        {
            norm = (onSolid) ? -s->SurfaceNormal(p) : -s->SurfaceNormal(s->ClosestSurfacePoint(p));
        }
    }
}


void G4HalfIntersectionSolid::Transform(const G4AffineTransform& t) const
{
    for (auto &s : intersectionSurfaces_) s->Transform(t);

    for (auto &s : subtractionSurfaces_) s->Transform(t);
}


G4GeometryType G4HalfIntersectionSolid::GetEntityType() const
{
    return "G4HalfIntersectionSolid";
}

G4bool G4HalfIntersectionSolid::CalculateExtent(const EAxis pAxis, const G4VoxelLimits& pVoxelLimit, const G4AffineTransform& pTransform, G4double& pMin, G4double& pMax) const
{
    return false;
}

std::ostream& G4HalfIntersectionSolid::StreamInfo(std::ostream& os) const
{
    os << "--- Half Intersection Solid: " << GetName() << " ---" << std::endl;
    os << "Surface Intersections: " << std::endl;
    for (const auto &s : intersectionSurfaces_)
    {
        s->StreamInfo(os);
    }
    os << "Surface Subtractions: " << std::endl;
    for (const auto &s : subtractionSurfaces_)
    {
        s->StreamInfo(os);
    }
    return os;
}

G4SurfaceMeshCGAL* G4HalfIntersectionSolid::GetSurfaceMesh() {

    G4cout << "G4HalfIntersectionSolid::GetSurfaceMesh" << G4endl;

    G4SurfaceMeshCGAL *sm1 = new G4SurfaceMeshCGAL();

    for(const auto &s : intersectionSurfaces_) {
        auto sm2 = s->GetSurfaceMesh();
        G4bool valid;
        sm1 = sm1->Intersection(sm2, valid);
    }
    for(const auto &s : subtractionSurfaces_) {
        auto sm2 = s->GetSurfaceMesh();
        G4bool valid;
        sm1 = sm1->Subtraction(sm2, valid);
    }

    return sm1;

}

void G4HalfIntersectionSolid::DescribeYourselfTo(G4VGraphicsScene& scene) const {
    G4cout << "G4HalfIntersectionSolid::DescribeYourselfTo" << std::endl;

    auto sm = const_cast<G4HalfIntersectionSolid*>(this)->GetSurfaceMesh();
    auto ph = sm->GetG4Polyhedron();
    scene.AddPrimitive(*ph);
}
