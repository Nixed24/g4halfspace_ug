#ifndef GEANT4_G4CADTOPOLOGY_HH
#define GEANT4_G4CADTOPOLOGY_HH

#include "TopoDS_Shape.hxx"

#include <cstddef>
#include <functional>

class G4CADTopology
{
public:
  using Visitor = std::function<void(const TopoDS_Shape&, std::size_t)>;

  explicit G4CADTopology(const TopoDS_Shape& shape);

  void Traverse(const Visitor& visitor) const;

private:
  static void TraverseShape(const TopoDS_Shape& shape,
                            std::size_t depth,
                            const Visitor& visitor);

  TopoDS_Shape shape_;
};

#endif
