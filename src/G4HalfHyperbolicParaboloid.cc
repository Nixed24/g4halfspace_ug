#include "G4HalfHyperbolicParaboloid.hh"
#include "G4SystemOfUnits.hh"

#include "G4PolynomialSolve.hh"
#include "G4ImplicitSurfaceCGAL.hh"

G4HalfHyperbolicParaboloid::G4HalfHyperbolicParaboloid(const G4double& a, const G4double& b, const G4double& c) : a_(a), b_(b), c_(c)
{
    if (a == 0 || b == 0 || c == 0)
        throw std::invalid_argument("G4HalfHyperbolicParaboloid: Constructor only takes non-zero values of a and b.");
}

G4HalfQuadric G4HalfHyperbolicParaboloid::GetQuadric() const
{
    CLHEP::HepMatrix Q(3,3,0);
    Q[0][0] = 1/(a_*a_);
    Q[1][1] = -1/(b_*b_);

    CLHEP::HepVector P(3,0);
    P[2] = -1/c_;

    G4HalfQuadric HyperbolicParaboloid(Q, P, 0);

    return HyperbolicParaboloid;
}

std::vector<G4double> G4HalfHyperbolicParaboloid::GetLambdas(const G4ThreeVector& p) const
{
    const G4double p5 = 1;
    const G4double p4 = -2*std::pow(a_, 2) + 2*std::pow(b_, 2) - 2*c_*p[2];
    const G4double p3 = std::pow(a_, 4) - 4*std::pow(a_, 2)*std::pow(b_, 2) + std::pow(b_, 4) + p[2]*(4*std::pow(a_, 2)*c_ - 4*std::pow(b_, 2)*c_);
    const G4double p2 = 2*std::pow(a_, 4)*std::pow(b_, 2) - 2*std::pow(a_, 2)*std::pow(b_, 4) + 2*std::pow(a_, 2)*c_*std::pow(p[0], 2) - 2*std::pow(b_, 2)*c_*std::pow(p[1], 2) + p[2]*(-2*std::pow(a_, 4)*c_ + 8*std::pow(a_, 2)*std::pow(b_, 2)*c_ - 2*std::pow(b_, 4)*c_);
    const G4double p1 = std::pow(a_, 4)*std::pow(b_, 4) + 4*std::pow(a_, 2)*std::pow(b_, 2)*c_*std::pow(p[0], 2) + 4*std::pow(a_, 2)*std::pow(b_, 2)*c_*std::pow(p[1], 2) + p[2]*(-4*std::pow(a_, 4)*std::pow(b_, 2)*c_ + 4*std::pow(a_, 2)*std::pow(b_, 4)*c_);
    const G4double p0 = -2*std::pow(a_, 4)*std::pow(b_, 4)*c_*p[2] - 2*std::pow(a_, 4)*std::pow(b_, 2)*c_*std::pow(p[1], 2) + 2*std::pow(a_, 2)*std::pow(b_, 4)*c_*std::pow(p[0], 2);

    const std::vector coefficients{p5, p4, p3, p2, p1, p0};
    const std::vector<G4double> lambdas = G4PolynomialSolve::SolvePolynomial(coefficients);
    if (lambdas.empty())
        throw std::logic_error("G4PolynomialSolve::SolvePolynomial: Lambda solve failed!");
    return lambdas;
}

std::vector<G4ThreeVector> G4HalfHyperbolicParaboloid::GetLMIntersections(const G4ThreeVector& p) const
{
    std::vector<G4double> lambdas = GetLambdas(p);
    std::vector<G4ThreeVector> points;
    for (const auto& lambda : lambdas)
    {
        const G4double determinant = (1 - lambda/(a_*a_))*(1 + lambda/(b_*b_));
        if (determinant == 0) continue;
        points.push_back({
            p[0]*a_*a_/(a_*a_ - lambda),
            p[1]*b_*b_/(b_*b_ + lambda),
            p[2] - lambda/(2*c_)
        });
    }
    return points;
}

G4double G4HalfHyperbolicParaboloid::ImplicitFunction(const G4ThreeVector& p) const
{
    return p[0]*p[0]/(a_*a_) - p[1]*p[1]/(b_*b_) - p[2]/c_;
}

std::vector<G4ThreeVector> G4HalfHyperbolicParaboloid::RayIntersection(
    const G4ThreeVector& p, const G4ThreeVector& v) const
{
    if (v.mag() == 0)
        throw std::logic_error("G4HalfHyperbolicParaboloid::RayIntersection: v magnitude cannot be 0.");
    const G4ThreeVector vUnit = v.unit();
    const G4double A = vUnit[0]*vUnit[0]/(a_*a_) - vUnit[1]*vUnit[1]/(b_*b_);
    const G4double B = 2*(vUnit[0]*p[0]/(a_*a_) - vUnit[1]*p[1]/(b_*b_)) - vUnit[2]/c_;
    const G4double C = ImplicitFunction(p);

    // Need this while trailing 0s are not removed at the start of the polynomial solver
    if (A == 0)
    {
        if (B == 0) return {};
        return {p - (C/B)*vUnit};
    }
    const std::vector coefficients{A, B, C};
    std::vector<G4double> distances = G4PolynomialSolve::SolvePolynomial(coefficients);
    std::sort(distances.begin(), distances.end());
    std::vector<G4ThreeVector> intersections;
    for (const auto distance : distances)
    {
        intersections.push_back(p + vUnit*distance);
    }

    return intersections;
}

G4double G4HalfHyperbolicParaboloid::DistanceToSurface(
    const G4ThreeVector& p, const G4ThreeVector& v) const
{
    const auto intersections = RayIntersection(p, v);
    if (intersections.empty()) return kInfinity;
    return (p - intersections[0]).mag();
}

G4double G4HalfHyperbolicParaboloid::DistanceToSurface(const G4ThreeVector& p) const
{
    const G4double sign = Inside(p) == kInside ? -1 : 1;
    return sign*(p - ClosestSurfacePoint(p)).mag();
}

void G4HalfHyperbolicParaboloid::GetExtent(
    G4ThreeVector& maxBounds, G4ThreeVector& minBounds, const G4ThreeVector& boundLimit) const
{
    maxBounds = boundLimit;
    minBounds = -boundLimit;
}

G4ThreeVector G4HalfHyperbolicParaboloid::SurfaceNormal(const G4ThreeVector& p) const
{
    return {2*p[0]/(a_*a_), -2*p[1]/(b_*b_), -1/c_};
}

G4ThreeVector G4HalfHyperbolicParaboloid::ClosestSurfacePoint(const G4ThreeVector& p) const
{
    if (Inside(p) == kSurface) return p;
    G4ThreeVector closestPoint;
    const auto intersections = GetLMIntersections(p);
    G4double distance = kInfinity;
    for (const auto& intersection : intersections)
    {
        const G4double candidateDistance = (p - intersection).mag();
        if (candidateDistance < distance)
        {
            distance = candidateDistance;
            closestPoint = intersection;
        }
    }
    return closestPoint;
}
G4HalfSurface *G4HalfHyperbolicParaboloid::Transform(const G4AffineTransform &t)
{
    G4HalfQuadric hyperbolicParaboloid = this->GetQuadric();

    G4HalfSurface* transformedPara = hyperbolicParaboloid.Transform(t);

    return transformedPara;
}

std::ostream& G4HalfHyperbolicParaboloid::StreamInfo(std::ostream& os) const
{
    os << ">>> G4HalfHyperbolicParaboloid <<<" << G4endl;
    os << "a = " << a_ << G4endl;
    os << "b = " << b_ << G4endl;
    return os;
}

G4SurfaceMeshCGAL* G4HalfHyperbolicParaboloid::GetSurfaceMesh()
{
    G4HalfQuadric hp = this->GetQuadric();
    return make_mesh(hp, 2.5*2.5*m*m);
}