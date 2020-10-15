/*!----------------------------------------------------------------------
 \file lib_headers_templates.cpp
 \brief header file handling for input_generator

 \maintainer Jonas Koban

 \level 3

 */
/*-----------------------------------------------------------*/

#ifndef LIB_HEADERS_TEMPLATES_CPP_
#define LIB_HEADERS_TEMPLATES_CPP_

#include "lib_headers.hpp"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

namespace Headers {

template<typename T> T Header::FromFile(std::string filePath) {
  std::string errorMessage = "-----File " + filePath + " not found.";
  std::ifstream inFile;
  inFile.open(filePath);
  T result;

  if (!inFile) {
    result = T();
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

template<typename T> T Header::FromHeaderFile(std::string headerName) {
  std::string path = Header::GetPath(headerName);
  return Header::FromFile<T>(path);
}

template<typename T> T Header::FromString(std::string data) {
  T result = T();
  result.ParseFromString(data);
  return result;
}

} // namespace Headers

#endif //LIB_HEADERS_TEMPLATES_CPP