/*!----------------------------------------------------------------------
 \file lib_headers.hpp
 \brief header file handling for input_generator

 \maintainer Jonas Koban

 \level 3

 */
/*-----------------------------------------------------------*/

#ifndef LIB_HEADERS_H_
#define LIB_HEADERS_H_

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include "../../libs/blocks/src/blocks.hpp"

using namespace Blocks;

namespace Headers {
enum class HeaderSection {
  ProblemTyp = 0,
  Discretisation,
  IO,
  StructuralDynamic,
  BinningStrategy,
  BrownianDynamics,
  BeamInteraction,
  BeamInteractionSphereBeamLink,
  BeamInteractionCrosslinking,
  BeamContact,
  BeamPotential,
  Solver,
  StructNox
};

static std::string HeaderSectionNames[] = {"PROBLEM_TYP",
                                           "DISCRETISATION",
                                           "IO",
                                           "STRUCTURAL_DYNAMIC",
                                           "BINNING_STRATEGY",
                                           "BROWNIAN_DYNAMICS",
                                           "BEAM_INTERACTION",
                                           "BEAM_INTERACTION_SPHERE_BEAM_LINK",
                                           "BEAM_INTERACTION_CROSSLINKING",
                                           "BEAM_CONTACT",
                                           "BEAM_POTENTIAL",
                                           "SOLVER_1",
                                           "STRUCT_NOX"};

class Header : public BlockKeyValueSectionCollection {
public:
  Header();
  ~Header();

  std::string ToString() const override;
  void ParseFromString(std::string &data) override;

  void StreamInto(std::stringstream &target);

  // overrides the key value pairs in all sections with the ones specified
  // (returns true), or adds them, if not present (returns false)
  void Override(Header header);

  template<typename T> static T FromFile(std::string filePath);
  template<typename T> static T FromHeaderFile(std::string headerName);
  template<typename T> static T FromString(std::string data);
  static std::string GetPath(std::string headerName);
  static std::string GetSectionName(HeaderSection section);
};

class ArgumentList : public BlockKeyValueCollection {
public:
  ArgumentList();
  ~ArgumentList();
  std::vector<std::string> ToStringVector() const;
  static ArgumentList FromFile(std::string filePath);
};

} // namespace Headers

// template implementation added for compiler
#include "lib_headers_templates.cpp"

#endif /* LIB_PROMPT_H_ */