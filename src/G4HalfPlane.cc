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


G4HalfSurface* G4HalfPlane::Transform(const G4AffineTransform& t)
{
    const G4ThreeVector& normalPoint = this->distance*this->normal;
    const G4ThreeVector& transformedPoint = t.TransformPoint(normalPoint);

    G4ThreeVector transformedNormal = t.TransformAxis(normal).unit();
    G4double transformedDistance = normal.dot(transformedPoint);

    return new G4HalfPlane(transformedNormal, transformedDistance);
}

void G4HalfPlane::GetExtent(G4ThreeVector &maxBounds, G4ThreeVector &minBounds, const G4ThreeVector &boundLimit) const
{
    const G4double maxRawBoundX = (distance + normal[1]*boundLimit[1] + normal[2]*boundLimit[2]) / boundLimit[0];
    const G4double maxRawBoundY = (distance + normal[0]*boundLimit[0] + normal[2]*boundLimit[2]) / boundLimit[1];
    const G4double maxRawBoundZ = (distance + normal[0]*boundLimit[0] + normal[1]*boundLimit[1]) / boundLimit[2];
    const G4double minRawBoundX = (distance - normal[1]*boundLimit[1] - normal[2]*boundLimit[2]) / boundLimit[0];
    const G4double minRawBoundY = (distance - normal[0]*boundLimit[0] - normal[2]*boundLimit[2]) / boundLimit[1];
    const G4double minRawBoundZ = (distance - normal[0]*boundLimit[0] - normal[1]*boundLimit[1]) / boundLimit[2];

    maxBounds[0] = std::min(boundLimit[0], maxRawBoundX);
    maxBounds[1] = std::min(boundLimit[1], maxRawBoundY);
    maxBounds[2] = std::min(boundLimit[2], maxRawBoundZ);
    minBounds[0] = std::max(-boundLimit[0], minRawBoundX);
    minBounds[1] = std::max(-boundLimit[1], minRawBoundY);
    minBounds[2] = std::max(-boundLimit[2], minRawBoundZ);
}

G4SurfaceMeshCGAL* G4HalfPlane::GetSurfaceMesh() {

    Nef_polyhedron_3_ECER nef = Nef_polyhedron_3_ECER(Nef_polyhedron_3_ECER::COMPLETE);
    nef *= Nef_polyhedron_3_ECER(Plane_3_ECER(Point_3_ECER(1000000000,0,0),
                                         Vector_3_ECER(1, 0, 0)));
    //G4cout << "G4HalfPlane::GetSurfaceMesh, bound 1, faces on nef = " << nef.number_of_facets() << G4endl;
    nef *= Nef_polyhedron_3_ECER(Plane_3_ECER(Point_3_ECER(-1000000000,0,0),
                                         Vector_3_ECER(-1, 0, 0)));
    //G4cout << "G4HalfPlane::GetSurfaceMesh, bound 2, faces on nef = " << nef.number_of_facets() << G4endl;
    nef *= Nef_polyhedron_3_ECER(Plane_3_ECER(Point_3_ECER(0,1000000000,0),
                                         Vector_3_ECER(0, 1, 0)));
    //G4cout << "G4HalfPlane::GetSurfaceMesh, bound 3, faces on nef = " << nef.number_of_facets() << G4endl;
    nef *= Nef_polyhedron_3_ECER(Plane_3_ECER(Point_3_ECER(0,-1000000000,0),
                                         Vector_3_ECER(0, -1, 0)));
    //G4cout << "G4HalfPlane::GetSurfaceMesh, bound 4, faces on nef = " << nef.number_of_facets() << G4endl;
    nef *= Nef_polyhedron_3_ECER(Plane_3_ECER(Point_3_ECER(0,0,1000000000),
                                         Vector_3_ECER(0, 0, 1)));
    //G4cout << "G4HalfPlane::GetSurfaceMesh, bound 5, faces on nef = " << nef.number_of_facets() << G4endl;
    nef *= Nef_polyhedron_3_ECER(Plane_3_ECER(Point_3_ECER(0,0,-1000000000),
                                         Vector_3_ECER(0, 0, -1)));
    //G4cout << "G4HalfPlane::GetSurfaceMesh, bound 6, faces on nef = " << nef.number_of_facets() << G4endl;

    nef *= Nef_polyhedron_3_ECER(Plane_3_ECER(Point_3_ECER(normal.x() * distance, normal.y() * distance, normal.z() * distance),
                                         Vector_3_ECER(normal.x(), normal.y(), normal.z())));
    auto DEBUG_MESH = G4SurfaceMeshCGAL(nef);
    /*
    G4cout << "G4HalfPlane::GetSurfaceMesh, add plane, faces on nef = " << nef.number_of_facets() <<
        ", faces on G4SurfaceMeshCGAL = " << DEBUG_MESH.NumberOfFaces() << G4endl;
    DEBUG_MESH.DebugOutput();
    G4cout << "Normal: (" << normal.x() << ", " << normal.y() << ", " << normal.z() << ")" << G4endl;
    for (int i = 0; i < DEBUG_MESH.NumberOfVertices(); ++i)
    {
        std::vector<G4double> vertex = DEBUG_MESH.GetVertex(i);
        G4cout << "Vertex " << i << ": (" << vertex[0] << ", " << vertex[1] << ", " << vertex[2] << ")" << G4endl;
    }
    */
    DEBUG_MESH.WriteMesh("final_mesh.off");
    return (new G4SurfaceMeshCGAL(DEBUG_MESH));
}