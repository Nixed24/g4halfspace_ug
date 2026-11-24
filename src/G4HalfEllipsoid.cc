#include "G4HalfEllipsoid.hh"

#include "G4PolynomialSolve.hh"
#include "G4ImplicitSurfaceCGAL.hh"

#include "G4Ellipsoid.hh"

G4HalfEllipsoid::G4HalfEllipsoid(const G4double &semiMajor,
                                 const G4double &semiMinor,
                                 const G4double &focus,
                                 const G4ThreeVector &centre):
semiMajor_(semiMajor), semiMinor_(semiMinor), focus_(focus), centre_(centre)
{}

G4HalfQuadric G4HalfEllipsoid::GetQuadric() const
{
    const G4double a2= semiMajor_*semiMajor_;
    const G4double b2 = semiMinor_*semiMinor_;
    const G4double c2 = focus_*focus_;

    CLHEP::HepMatrix Q(3,3,0);
    Q[0][0] = 1/a2;
    Q[1][1] = 1/b2;
    Q[2][2] = 1/c2;

    CLHEP::HepVector P(3,0);
    P[0] = -2*centre_[0]/a2;
    P[1] = -2*centre_[1]/b2;
    P[2] = -2*centre_[2]/c2;

    const G4double x02 = centre_[0]*centre_[0];
    const G4double y02 = centre_[1]*centre_[1];
    const G4double z02 = centre_[2]*centre_[2];

    const G4double R = x02/a2 + y02/b2 + z02/c2 - 1;

    const G4HalfQuadric ellipsoid(Q, P, R);

    return ellipsoid;
}

std::vector<G4double> G4HalfEllipsoid::GetLambdas(const G4ThreeVector& p) const
{
    const G4double rx = (p - centre_)[0];
    const G4double ry = (p - centre_)[1];
    const G4double rz = (p - centre_)[2];

    const G4double a = semiMajor_;
    const G4double b = semiMinor_;
    const G4double c = focus_;

    const G4double p6 = -1;
    const G4double p5 = 2*std::pow(a, 2) + 2*std::pow(b, 2) + 2*std::pow(c, 2);
    const G4double p4 = -std::pow(a, 4) - 4*std::pow(a, 2)*std::pow(b, 2) - 4*std::pow(a, 2)*std::pow(c, 2)
                        + std::pow(a, 2)*std::pow(rx, 2) - std::pow(b, 4) - 4*std::pow(b, 2)*std::pow(c, 2)
                        + std::pow(b, 2)*std::pow(ry, 2) - std::pow(c, 4) + std::pow(c, 2)*std::pow(rz, 2);
    const G4double p3 = 2*std::pow(a, 4)*std::pow(b, 2) + 2*std::pow(a, 4)*std::pow(c, 2)
                        + 2*std::pow(a, 2)*std::pow(b, 4) + 8*std::pow(a, 2)*std::pow(b, 2)*std::pow(c, 2)
                        + 2*std::pow(a, 2)*std::pow(c, 4) + 2*std::pow(b, 4)*std::pow(c, 2)
                        + 2*std::pow(b, 2)*std::pow(c, 4) + std::pow(rx, 2)*(-2*std::pow(a, 2)*std::pow(b, 2)
                        - 2*std::pow(a, 2)*std::pow(c, 2)) + std::pow(ry, 2)*(-2*std::pow(a, 2)*std::pow(b, 2)
                        - 2*std::pow(b, 2)*std::pow(c, 2)) + std::pow(rz, 2)*(-2*std::pow(a, 2)*std::pow(c, 2)
                        - 2*std::pow(b, 2)*std::pow(c, 2));
    const G4double p2 = -std::pow(a, 4)*std::pow(b, 4) - 4*std::pow(a, 4)*std::pow(b, 2)*std::pow(c, 2)
                        - std::pow(a, 4)*std::pow(c, 4) - 4*std::pow(a, 2)*std::pow(b, 4)*std::pow(c, 2)
                        - 4*std::pow(a, 2)*std::pow(b, 2)*std::pow(c, 4) - std::pow(b, 4)*std::pow(c, 4)
                        + std::pow(rx, 2)*(std::pow(a, 2)*std::pow(b, 4) + 4*std::pow(a, 2)*std::pow(b, 2)*std::pow(c, 2)
                        + std::pow(a, 2)*std::pow(c, 4)) + std::pow(ry, 2)*(std::pow(a, 4)*std::pow(b, 2)
                        + 4*std::pow(a, 2)*std::pow(b, 2)*std::pow(c, 2) + std::pow(b, 2)*std::pow(c, 4))
                        + std::pow(rz, 2)*(std::pow(a, 4)*std::pow(c, 2) + 4*std::pow(a, 2)*std::pow(b, 2)*std::pow(c, 2)
                        + std::pow(b, 4)*std::pow(c, 2));
    const G4double p1 = 2*std::pow(a, 4)*std::pow(b, 4)*std::pow(c, 2) + 2*std::pow(a, 4)*std::pow(b, 2)*std::pow(c, 4)
                        + 2*std::pow(a, 2)*std::pow(b, 4)*std::pow(c, 4) + std::pow(rx, 2)*(-2*std::pow(a, 2)*std::pow(b, 4)*std::pow(c, 2)
                        - 2*std::pow(a, 2)*std::pow(b, 2)*std::pow(c, 4)) + std::pow(ry, 2)*(-2*std::pow(a, 4)*std::pow(b, 2)*std::pow(c, 2)
                        - 2*std::pow(a, 2)*std::pow(b, 2)*std::pow(c, 4)) + std::pow(rz, 2)*(-2*std::pow(a, 4)*std::pow(b, 2)*std::pow(c, 2)
                        - 2*std::pow(a, 2)*std::pow(b, 4)*std::pow(c, 2));
    const G4double p0 = -std::pow(a, 4)*std::pow(b, 4)*std::pow(c, 4)
                        + std::pow(a, 4)*std::pow(b, 4)*std::pow(c, 2)*std::pow(rz, 2)
                        + std::pow(a, 4)*std::pow(b, 2)*std::pow(c, 4)*std::pow(ry, 2)
                        + std::pow(a, 2)*std::pow(b, 4)*std::pow(c, 4)*std::pow(rx, 2);

    const std::vector coefficients{p6, p5, p4, p3, p2, p1, p0};
    const std::vector<G4double> lambdas = G4PolynomialSolve::SolvePolynomial(coefficients);

    if (lambdas.empty())
    {
        throw std::logic_error("G4PolynomialSolve::SolvePolynomial: Lambda solve failed!");
    }

    return lambdas;
}

G4double G4HalfEllipsoid::ImplicitFunction(const G4ThreeVector &p) const
{
    const G4ThreeVector q = p - centre_;

    const G4double x = q[0]*q[0] / (semiMajor_*semiMajor_);
    const G4double y = q[1]*q[1] / (semiMinor_*semiMinor_);
    const G4double z = q[2]*q[2] / (focus_*focus_);

    return x + y + z - 1;
}

std::vector<G4ThreeVector> G4HalfEllipsoid::RayIntersection(const G4ThreeVector &p, const G4ThreeVector &v) const
{
    if (v.mag() == 0)
        throw std::logic_error("G4HalfEllipsoid::RayIntersection: Cannot pass 0 length direction vector.");

    const G4ThreeVector q = p - centre_;

    G4ThreeVector vUnit = v.unit();

    const G4double a2 = semiMajor_*semiMajor_;
    const G4double b2 = semiMinor_*semiMinor_;
    const G4double c2 = focus_*focus_;
    const G4double vx2 = vUnit[0]*vUnit[0];
    const G4double vy2 = vUnit[1]*vUnit[1];
    const G4double vz2 = vUnit[2]*vUnit[2];
    const G4double qx2 = q[0]*q[0];
    const G4double qy2 = q[1]*q[1];
    const G4double qz2 = q[2]*q[2];

    const G4double A  = vx2*b2*c2 + vy2*a2*c2 + vz2*a2*b2;
    const G4double B = 2*(q[0]*vUnit[0]*b2*c2 + q[1]*vUnit[1]*a2*c2 + q[2]*vUnit[2]*a2*b2);
    const G4double C = qx2*b2*c2 + qy2*a2*c2 + qz2*a2*b2 - a2*b2*c2;

    const G4double determinant = B*B - 4*A*C;
    if (determinant == 0) return {p + (-B/(2*A))*vUnit};
    if (determinant < 0) return {};

    const G4double minDistance = (-B - sqrt(determinant))/(2*A);
    const G4double maxDistance = (-B + sqrt(determinant))/(2*A);

    return {p + minDistance*vUnit, p + maxDistance*vUnit};
}

G4double G4HalfEllipsoid::DistanceToSurface(const G4ThreeVector &p, const G4ThreeVector &v) const
{
    if (v.mag() == 0)
        throw std::logic_error("G4HalfEllipsoid::DistanceToSurface: Cannot pass 0 length direction vector.");

    return (p - RayIntersection(p, v)[0]).mag();
}

G4double G4HalfEllipsoid::DistanceToSurface(const G4ThreeVector &p) const
{
    const G4ThreeVector q = p - centre_;
    const std::vector<G4double> lambdas = GetLambdas(p);
    G4double signedDistance = kInfinity;

    for (const auto &lambda : lambdas)
    {
        G4double candidateDistance = lambda * std::sqrt(
              (q[0]*q[0] / ((semiMajor_*semiMajor_ - lambda)*(semiMajor_*semiMajor_ - lambda)))
            + (q[1]*q[1] / ((semiMinor_*semiMinor_ - lambda)*(semiMinor_*semiMinor_ - lambda)))
            + (q[2]*q[2] / ((focus_*focus_ - lambda)*(focus_*focus_ - lambda))));

        if (candidateDistance < signedDistance)
            signedDistance = candidateDistance;
    }

    return signedDistance;
}

void G4HalfEllipsoid::GetExtent(G4ThreeVector &maxBounds, G4ThreeVector &minBounds, const G4ThreeVector &boundLimit) const
{
    maxBounds =
    {
        std::min(boundLimit[0], centre_[0] + semiMajor_),
        std::min(boundLimit[1], centre_[1] + semiMinor_),
        std::min(boundLimit[2], centre_[2] + focus_)
    };

    minBounds =
    {
        std::max(-boundLimit[0], centre_[0] - semiMajor_),
        std::max(-boundLimit[1], centre_[1] - semiMinor_),
        std::max(-boundLimit[2], centre_[2] - focus_)
    };
}

G4ThreeVector G4HalfEllipsoid::ClosestSurfacePoint(const G4ThreeVector &p) const
{
    const G4ThreeVector q = p - centre_;
    std::vector<G4double> lambdas = GetLambdas(p);
    const G4double lambda = *std::min_element(std::begin(lambdas), std::end(lambdas));

    return G4ThreeVector
    {
        q[0]*(semiMajor_*semiMajor_/(semiMajor_*semiMajor_ - lambda)),
        q[1]*(semiMinor_*semiMinor_/(semiMinor_*semiMinor_ - lambda)),
        q[2]*(focus_*focus_/(focus_*focus_ - lambda))
    };
}

G4ThreeVector G4HalfEllipsoid::SurfaceNormal(const G4ThreeVector &p) const
{
    G4ThreeVector surfacePoint = p - centre_;
    if (this->Inside(p) != kSurface)
    {
        surfacePoint = ClosestSurfacePoint(p);
    }

    return G4ThreeVector
    {
        2*surfacePoint[0]/(semiMajor_*semiMajor_),
        2*surfacePoint[1]/(semiMinor_*semiMinor_),
        2*surfacePoint[2]/(focus_*focus_)
    };
}

G4HalfSurface* G4HalfEllipsoid::Transform(const G4AffineTransform& t)
{
    G4HalfQuadric ellipsoid = this->GetQuadric();

    G4HalfSurface* transformedEllipsoid = ellipsoid.Transform(t);

    return transformedEllipsoid;
}

G4SurfaceMeshCGAL* G4HalfEllipsoid::GetSurfaceMesh()
{
    const auto t = G4Ellipsoid("test", semiMajor_,semiMinor_, focus_);
    G4Polyhedron *g4poly = t.GetPolyhedron();
    auto *sm = new G4SurfaceMeshCGAL();
    sm->Fill(g4poly);

    G4ThreeVector axis;
    G4double angle;

    sm->Translate(centre_);

    return sm;
}

std::ostream &G4HalfEllipsoid::StreamInfo(std::ostream &os) const
{
    os << ">>> G4HalfEllipsoid <<<" << G4endl;
    os << "Semi Major = " << semiMajor_ << G4endl;
    os << "Semi Minor = " << semiMinor_ << G4endl;
    os << "Focus = " << focus_ << G4endl;
    os << "Centre = " << centre_ << G4endl;

    return os;
}
