#include "G4HalfSurface.hh"


EInside G4HalfSurface::Inside(const G4ThreeVector& p) const
{
    const G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

    if (this->ImplicitFunction(p) < -kCarTolerance) return kInside;
    if (this->ImplicitFunction(p) > kCarTolerance) return kOutside;
    return kSurface;
}

std::ostream& G4HalfSurface::StreamInfo(std::ostream& os) const
{
    os << "--- Half Surface ---" << std::endl;
    return os;
}