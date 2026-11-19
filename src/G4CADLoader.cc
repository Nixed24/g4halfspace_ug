#include "G4CADLoader.hh"

#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <IGESControl_Reader.hxx>
#include <STEPControl_Reader.hxx>

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>

namespace
{
std::string FileExtension(const G4String& filename)
{
  const std::string path = filename;
  const std::size_t dot = path.find_last_of('.');
  if (dot == std::string::npos) {
    return {};
  }

  std::string extension = path.substr(dot + 1);
  std::transform(extension.begin(), extension.end(), extension.begin(),
                 [](unsigned char character) {
                   return static_cast<char>(std::tolower(character));
                 });
  return extension;
}
}  // namespace

TopoDS_Shape G4CADLoader::Load(const G4String& filename) const
{
  if (filename.empty()) {
    throw std::invalid_argument("G4CADLoader requires a CAD filename");
  }

  const std::string extension = FileExtension(filename);
  if (extension == "step" || extension == "stp") {
    STEPControl_Reader reader;
    if (reader.ReadFile(filename.c_str()) != IFSelect_RetDone ||
        reader.TransferRoots() == 0) {
      throw std::runtime_error("Could not read STEP file: " +
                               std::string(filename));
    }
    return reader.OneShape();
  }

  if (extension == "iges" || extension == "igs") {
    IGESControl_Reader reader;
    if (reader.ReadFile(filename.c_str()) != IFSelect_RetDone ||
        reader.TransferRoots() == 0) {
      throw std::runtime_error("Could not read IGES file: " +
                               std::string(filename));
    }
    return reader.OneShape();
  }

  if (extension == "brep" || extension == "brp") {
    BRep_Builder builder;
    TopoDS_Shape shape;
    if (!BRepTools::Read(shape, filename.c_str(), builder)) {
      throw std::runtime_error("Could not read BREP file: " +
                               std::string(filename));
    }
    return shape;
  }

  throw std::invalid_argument("Unsupported CAD file extension: ." +
                              extension);
}
