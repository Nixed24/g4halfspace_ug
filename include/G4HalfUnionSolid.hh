#ifndef G4HALF_UNION_SOLID_H
#define G4HALF_UNION_SOLID_H

#include "G4VSolid.hh"
#include "G4HalfSolid.hh"
#include "G4HalfIntersectionSolid.hh"
#include "G4SolidStore.hh"
#include <vector>
#include <utility>
#include <memory>

class G4HalfUnionSolid : public G4HalfSolid
{
    std::vector<G4HalfSolid*> unionSolids_;

protected:
    void GetOnSurfaceNormal(const G4ThreeVector& p, G4ThreeVector& norm) const override;

public:
    using G4HalfSolid::DistanceToIn;
    using G4HalfSolid::DistanceToOut;

    explicit G4HalfUnionSolid(const G4String& name = "G4HalfIntersectionSolid");
    G4HalfUnionSolid(const std::vector<G4HalfSolid*>& unionSolids,
                    const G4String& name = "G4HalfIntersectionSolid");

    ~G4HalfUnionSolid() override;

    void AddUnionSolid(G4HalfSolid* unionSolid);

    // --- Inherited functions ---
    /*
     * Calculates the minimum and maximum length of the solid
     * along a given axis, clipped to pVoxelLimit.
     *
     * Used for rendering.
     *
     * Leaving this virtual for now until stewart gives me his renderer
     */
    G4bool CalculateExtent(const EAxis pAxis, const G4VoxelLimits& pVoxelLimit, const G4AffineTransform& pTransform,
                           G4double& pMin, G4double& pMax) const override {return false;}

    G4double ImplicitFunction(const G4ThreeVector &p) const override;

    void GetPotentialRayIntersections(const G4ThreeVector&, const G4ThreeVector&,
        std::vector<RayHit>&, const G4bool) const override;

    G4double SignedDistance(const G4ThreeVector& p, G4ThreeVector* n, const G4bool calcNorm) const override;

    /*
     * Calculates shortest distance to the inside of the solid.
     */
    G4double DistanceToIn(const G4ThreeVector& p) const override;

    /*
     * Calculates shortest distance to the inside of the solid,
     * passing the closest surface normal by reference
     */
    G4double DistanceToIn(const G4ThreeVector& p, const G4bool calcNorm, G4ThreeVector* n) const override;

    G4double DistanceToOut(const G4ThreeVector& p) const override;
    G4double DistanceToOut(const G4ThreeVector& p, const G4bool calcNorm, G4ThreeVector* n) const override;

    G4HalfSolid* Transform(const G4AffineTransform& t) const override;

    G4GeometryType GetEntityType() const override;
    void DescribeYourselfTo(G4VGraphicsScene& scene) const override;
    std::ostream& StreamInfo(std::ostream& os) const override;
    G4SurfaceMeshCGAL* GetSurfaceMesh() override;
};

#endif
