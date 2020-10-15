/*!----------------------------------------------------------------------
 \file lib_headers.cpp
 \brief header file handling for input_generator

 \maintainer Jonas Koban

 \level 3

 */
/*-----------------------------------------------------------*/

#include "lib_headers.hpp"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

namespace Headers {
Header::Header() : BlockKeyValueSectionCollection() {}
Header::~Header() {}

void Header::ParseFromString(std::string &data) {
  BlockKeyValueSectionCollection::ParseFromString(data);
}

std::string Header::ToString() const {
  return BlockKeyValueSectionCollection::ToString();
}

void Header::StreamInto(std::stringstream &target) {
  target << this->ToString();
}

void Header::Override(Header header) {
  for (size_t sectionId = 0; sectionId < header.children_.size(); ++sectionId)
    this->OverrideSection(header.children_[sectionId]);
}

std::string Header::GetPath(std::string headerName) {
  return "./Input/" + headerName + ".txt";
}

std::string Header::GetSectionName(HeaderSection section) {
  return HeaderSectionNames[(int)section];
}

ArgumentList::ArgumentList(){};
ArgumentList::~ArgumentList(){};

ArgumentList ArgumentList::FromFile(std::string filePath) {
  std::string errorMessage = "File:" + filePath + " not found.";
  std::ifstream inFile;
  inFile.open(filePath);
  ArgumentList result;

  if (!inFile) {
    result = ArgumentList();
    result.ParseFromString(errorMessage);
  } else {
    std::stringstream buffer;
    buffer << inFile.rdbuf();
    std::string data = buffer.str();
    result.ParseFromString(data);
  }

  inFile.close();
  return result;
}

std::vector<std::string> ArgumentList::ToStringVector() const {
  size_t childCount = this->children_.size();
  std::vector<std::string> result(childCount);
  for (size_t childId = 0; childId < childCount; ++childId)
    result[childId] = this->children_[childId].GetKey();
  // remove comments
  for (size_t childId = 0; childId < childCount; ++childId)
    if (result[childId].length() > 1 && result[childId].substr(0, 2) == "//") {
      result.erase(result.begin() + childId);
      --childId;
      --childCount;
    }

  return result;
}

} // namespace Headers
