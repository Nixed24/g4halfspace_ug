#include "G4HalfPlane.hh"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Extended_cartesian.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include "G4ImplicitSurfaceCGAL.hh"

G4HalfPlane::G4HalfPlane(const G4ThreeVector& normalIn, const G4double& distanceIn)
            :normal(normalIn.unit()),
             distance(distanceIn) {}


G4double G4HalfPlane::ImplicitFunction(const G4ThreeVector& p) const
{
    return normal.dot(p) - distance;
}


std::vector<G4ThreeVector> G4HalfPlane::RayIntersection(const G4ThreeVector& p, const G4ThreeVector& v) const
{
    const G4ThreeVector vUnit = v.unit();

    const G4double distanceToSurface = this->DistanceToSurface(p,v);

    if (distanceToSurface == kInfinity && this->Inside(p) == kSurface) return {p};
    if (distanceToSurface == kInfinity) return {};

    return {vUnit*distanceToSurface + p};
}


// Returns the signed distance to surface
G4double G4HalfPlane::DistanceToSurface(const G4ThreeVector& p) const
{
    // For a plane, the implicit function is the same as the signed distance
    return this->ImplicitFunction(p);
}


G4double G4HalfPlane::DistanceToSurface(const G4ThreeVector& p, const G4ThreeVector& v) const
{
    const G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

    const G4ThreeVector vUnit = v.unit();

    G4double projection = normal.dot(vUnit);
    G4double implicitFunction = this->ImplicitFunction(p);

    if (std::abs(implicitFunction) <= kCarTolerance) return 0;

    if (projection == 0) return kInfinity;
    if (projection < 0 && implicitFunction < -kCarTolerance) return kInfinity;
    if (projection > 0 && implicitFunction > kCarTolerance) return kInfinity;

    return -implicitFunction / projection;
    if (projection == 0) return kInfinity; // Parallel to plane

    return -implicitFunction/projection;
}


G4ThreeVector G4HalfPlane::ClosestSurfacePoint(const G4ThreeVector& p) const
{
    return p - this->ImplicitFunction(p)*this->normal;
}


void G4HalfPlane::Transform(const G4AffineTransform& t)
{
    const G4ThreeVector& normalPoint = this->distance*this->normal;
    const G4ThreeVector& transformedPoint = t.TransformPoint(normalPoint);

    this->normal = t.TransformAxis(normal).unit();
    this->distance = normal.dot(transformedPoint);
}


G4SurfaceMeshCGAL* G4HalfPlane::GetSurfaceMesh() {

    Nef_polyhedron_3_ECER nef = Nef_polyhedron_3_ECER(Nef_polyhedron_3_ECER::COMPLETE);
    nef *= Nef_polyhedron_3_ECER(Plane_3_ECER(Point_3_ECER(1000000000,0,0),
                                         Vector_3_ECER(1, 0, 0)));
    nef *= Nef_polyhedron_3_ECER(Plane_3_ECER(Point_3_ECER(-1000000000,0,0),
                                         Vector_3_ECER(-1, 0, 0)));
    nef *= Nef_polyhedron_3_ECER(Plane_3_ECER(Point_3_ECER(0,1000000000,0),
                                         Vector_3_ECER(0, 1, 0)));
    nef *= Nef_polyhedron_3_ECER(Plane_3_ECER(Point_3_ECER(0,-1000000000,0),
                                         Vector_3_ECER(0, -1, 0)));
    nef *= Nef_polyhedron_3_ECER(Plane_3_ECER(Point_3_ECER(0,0,1000000000),
                                         Vector_3_ECER(0, 0, 1)));
    nef *= Nef_polyhedron_3_ECER(Plane_3_ECER(Point_3_ECER(0,0,-1000000000),
                                         Vector_3_ECER(0, 0, -1)));

    nef *= Nef_polyhedron_3_ECER(Plane_3_ECER(Point_3_ECER(normal.x() * distance, normal.y() * distance, normal.z() * distance),
                                         Vector_3_ECER(normal.x(), normal.y(), normal.z())));

    return (new G4SurfaceMeshCGAL(nef));
}
