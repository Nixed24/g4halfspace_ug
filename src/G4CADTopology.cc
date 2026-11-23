#include "G4CADTopology.hh"

#include <TopoDS_Iterator.hxx>

#include <stdexcept>

G4CADTopology::G4CADTopology(const TopoDS_Shape& shape) : shape_(shape)
{
  if (shape_.IsNull()) {
    throw std::invalid_argument("G4CADTopology requires a non-null shape");
  }
}

void G4CADTopology::Traverse(const Visitor& visitor) const
{
  if (!visitor) {
    throw std::invalid_argument("G4CADTopology requires a traversal visitor");
  }

  TraverseShape(shape_, 0, visitor);
}

void G4CADTopology::TraverseShape(const TopoDS_Shape& shape,
                                  const std::size_t depth,
                                  const Visitor& visitor)
{
  visitor(shape, depth);

  for (TopoDS_Iterator iterator(shape); iterator.More(); iterator.Next()) {
    TraverseShape(iterator.Value(), depth + 1, visitor);
  }
}
