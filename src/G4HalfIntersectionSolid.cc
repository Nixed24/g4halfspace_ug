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
    G4double intersectionImplicitFunction = -kInfinity;
    G4double subtractionImplicitFunction = -kInfinity;

    for (const auto &s : intersectionSurfaces_)
    {
        intersectionImplicitFunction = std::max(intersectionImplicitFunction, s->ImplicitFunction(p));
    }

    for (const auto &s : subtractionSurfaces_)
    {
        subtractionImplicitFunction = std::max(subtractionImplicitFunction, s->ImplicitFunction(p));
    }

    if (subtractionSurfaces_.empty()) subtractionImplicitFunction = kInfinity;

    return std::max(intersectionImplicitFunction, -subtractionImplicitFunction);
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

            const G4double distance = (intersection - p).dot(vUnit);
            if (distance < -kCarTolerance) continue; // Keep only forward-ray intersectionsinside

            RayHit hit;
            hit.intersection = intersection;
            hit.distance = std::max(0.0, distance);
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

void G4HalfIntersectionSolid::SetDistance(G4double &distance, G4int &normSign, const G4ThreeVector &p,
    const G4HalfSurface &surface, G4HalfSurface &candidateSurface, const G4int &sign)
{
    G4double candidateDistance = surface.DistanceToSurface(p);

    if (sign*candidateDistance < distance)
    {
        distance = candidateDistance;
        candidateSurface = surface;
        normSign = sign;
    }
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
        //SetDistance(signedDistance, normSign, p, *surface, *candidateSurface,  +1);
        setDistance(surface, +1);
    }

    for (const auto &surface : subtractionSurfaces_)
    {
        //SetDistance(signedDistance, normSign, p, *surface, *candidateSurface,  -1);
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


G4HalfSolid* G4HalfIntersectionSolid::Transform(const G4AffineTransform& t) const
{
    std::vector<G4HalfSurface*> transformedIntersectionSurfaces;
    std::vector<G4HalfSurface*> transformedSubtractionSurfaces;

    transformedIntersectionSurfaces.reserve(intersectionSurfaces_.size());
    transformedSubtractionSurfaces.reserve(subtractionSurfaces_.size());

    for (const auto &s : intersectionSurfaces_)
    {
        transformedIntersectionSurfaces.push_back(s->Transform(t));
    }

    for (const auto &s : subtractionSurfaces_)
    {
        transformedSubtractionSurfaces.push_back(s->Transform(t));
    }

    return new G4HalfIntersectionSolid(transformedIntersectionSurfaces, transformedSubtractionSurfaces, this->GetName());
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
    G4bool first = true;
    auto intersectionMesh = new G4SurfaceMeshCGAL();
    auto subtractionMesh = new G4SurfaceMeshCGAL();
    auto sm = new G4SurfaceMeshCGAL();
    G4bool valid;
    for(const auto &s : intersectionSurfaces_) {
        if (first)
        {
            intersectionMesh = s->GetSurfaceMesh();
            first = false;
            continue;
        }
        auto intersectionMesh2 = s->GetSurfaceMesh();
        intersectionMesh = intersectionMesh->Intersection(intersectionMesh2, valid);
        //G4cout << "G4HalfIntersectionSolid::GetSurfaceMesh post-intersection facets = " << sm1->NumberOfFaces() << G4endl;
        //G4cout << "valid variable output : " << valid << G4endl;
    }
    first = true;
    for(const auto &s : subtractionSurfaces_) {
        if (first)
        {
            subtractionMesh = s->GetSurfaceMesh();
            first = false;
            continue;
        }
        auto subtractionMesh2 = s->GetSurfaceMesh();
        subtractionMesh = subtractionMesh->Intersection(subtractionMesh2, valid);
        G4cout << "G4HalfIntersectionSolid::GetSurfaceMesh post-subtraction facets = " << subtractionMesh->NumberOfFaces() << G4endl;
        G4cout << "valid variable output : " << valid << G4endl;
    }

    if (subtractionSurfaces_.empty())
        return intersectionMesh;

    sm = intersectionMesh->Subtraction(subtractionMesh, valid);

    return sm;

}

void G4HalfIntersectionSolid::DescribeYourselfTo(G4VGraphicsScene& scene) const {
    G4cout << "G4HalfIntersectionSolid::DescribeYourselfTo" << std::endl;

    auto sm = const_cast<G4HalfIntersectionSolid*>(this)->GetSurfaceMesh();
    auto ph = sm->GetG4Polyhedron();
    scene.AddPrimitive(*ph);
}
