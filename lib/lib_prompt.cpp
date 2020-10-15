/*!----------------------------------------------------------------------
 \file lib_prompt.cpp
 \brief prompt functions for input_generator

 \maintainer Jonas Koban

 \level 3

 */
/*-----------------------------------------------------------*/
#include "lib_prompt.hpp"
#include "filanetwork_utils.hpp"
#include <string>
#include <vector>
#include <cmath>
#include <string>

bool PromptYesNo(std::string message)
{
  std::cout << message + "\n";
  std::string readString;
  while (true)
  {
    std::cout << "                             1/y(es)   0/n(o): ";
    getline(std::cin, readString);
    if (readString == "y" || readString == "1")
      return true;
    else if (readString == "n" || readString == "0")
      return false;
    else
      std::cout
          << "Input \'" + readString + "\' invalid. " + "Try again!\n";
  }
} // PromptYesNo()
bool PromptYesNoCMD(std::string message, std::vector<std::string> &arguments)
{

  if (arguments.size() == 0)
    return PromptYesNo(message);
  else
  {
    std::cout << message + "\n";
    std::cout << "                             1/y(es)   0/n(o): ";
    std::cout << arguments[0] << "\n";
    std::string arg = arguments[0];
    arguments.erase(arguments.begin());
    if (arg == "y" || arg == "1")
      return true;
    else if (arg == "n" || arg == "0")
      return false;
    else
    {
      std::cout << "Argument '" << arg << "' invalid.\n";
      return PromptYesNo(message);
    }
  }
}
unsigned int Prompt10(std::string message)
{
  std::cout << message + "\n";
  std::string readString;
  while (true)
  {
    std::cout << "                             1/y(es)   0/n(o): ";
    getline(std::cin, readString);
    if (readString == "y" || readString == "1")
      return 1;
    else if (readString == "n" || readString == "0")
      return 0;
    else
      std::cout
          << "Input \'" + readString + "\' invalid. " + "Try again!\n";
  }
} // PromptYesNo()
unsigned int Prompt10CMD(std::string message, std::vector<std::string> &arguments)
{

  if (arguments.size() == 0)
    return Prompt10(message);
  else
  {
    std::cout << message + "\n";
    std::cout << "                             1/y(es)   0/n(o): ";
    std::cout << arguments[0] << "\n";
    std::string arg = arguments.at(0);
    arguments.erase(arguments.begin());
    if (arg == "y" || arg == "1")
      return 1;
    else if (arg == "n" || arg == "0")
      return 0;
    else
    {
      std::cout << "Argument '" << message << " = " << arg << "' invalid.\n";
      return Prompt10(message);
    }
  }
}
int PromptInt(std::string message, std::string name, bool includeZero)
{
  std::cout << message + "\n";
  int param;
  int min = 1;
  if (includeZero)
    min = 0;

  std::string readString;
  while (true)
  {
    std::cout << "                                          " << min << "-...: ";
    getline(std::cin, readString);
    // This code converts from string to number safely.
    std::stringstream stream(readString);
    if (stream >> param)
    {
      if (param < min)
      {
        std::cout << "Invalid input, " << name << "=" << param << "<"
                  << min << "! Please reenter value: ";
        continue;
      }
      else
        break;
    }
    else
      std::cout << "Invalid input, " << name << "=" << readString
                << "! Please reenter value: ";
  }
  return param;
} // PromptInt()
int PromptIntCMD(std::string message, std::string name, std::vector<std::string> &arguments,
                 bool includeZero)
{
  if (arguments.size() == 0)
    return PromptInt(message, name, includeZero);
  else
  {
    std::string arg = arguments.at(0);
    int param;
    std::stringstream stream(arg);
    std::cout << message + "\n";
    std::cout << arguments[0] << "\n";
    arguments.erase(arguments.begin());
    if (stream >> param)
    {
      return param;
    }
    else
    {
      std::cout << "Argument '" << message << " = " << arg << "' invalid.\n";
      return PromptInt(message, name, includeZero);
    }
  }
}
unsigned int PromptUInt(std::string message, std::string name, bool includeZero)
{
  std::cout << message + "\n";
  int param;

  std::string readString;
  while (true)
  {
    getline(std::cin, readString);
    // This code converts from string to number safely.
    std::stringstream stream(readString);
    if (stream >> param)
      if (includeZero || param != 0)
        break;
      else
        std::cout << "Input '" << message << " = " << param
                  << "' out of bounds. Minimum: "
                  << 1 << ". Please reenter value: \n";
    else
      std::cout << "Invalid input, " << name << "=" << readString
                << "! Please reenter value: ";
  }
  return param;
} // PromptInt()
unsigned int PromptUIntCMD(std::string message, std::string name, std::vector<std::string> &arguments,
                           bool includeZero)
{
  if (arguments.size() == 0)
    return PromptInt(message, name, includeZero);
  else
  {
    std::string arg = arguments.at(0);
    int param;
    int min = 1;
    if (includeZero)
      min = 0;
    std::stringstream stream(arg);
    std::cout << message + "\n";
    std::cout << "                                        " << min << "-...: ";
    std::cout << arguments[0] << "\n";
    arguments.erase(arguments.begin());
    if (stream >> param)
    {
      if (param < min)
      {
        std::cout << "Argument '" << message << " = " << arg
                  << "' out of bounds. Minimum: "
                  << min << "\n";
        return PromptInt(message, name, includeZero);
      }
      else
        return param;
    }
    else
    {
      std::cout << "Argument '" << message << " = " << arg << "' invalid.\n";
      return PromptInt(message, name, includeZero);
    }
  }
}
double PromptDouble(std::string message, std::string name, double minimum, double maximum)
{
  std::cout << message + "\n";
  double param;
  std::string stringRead;
  while (true)
  {
    std::cout << "                                   " << minimum << "-" << maximum << ": ";
    getline(std::cin, stringRead);
    // This code converts from string to number safely.
    std::stringstream stream(stringRead);
    if (stream >> param)
    {
      if (param < minimum || param > maximum)
      {
        std::cout << "Invalid input, " << name << "=" << param << "<"
                  << minimum << "! Please reenter value: ";
        continue;
      }
      break;
    }
    else
      std::cout << "Invalid input, " << name << "=" << stringRead
                << "! Please reenter value: ";
  }
  return param;
}
double PromptDoubleCMD(std::string message, std::string name, std::vector<std::string> &arguments,
                       double minimum, double maximum)
{
  if (arguments.size() == 0)
    return PromptDouble(message, name, minimum, maximum);
  else
  {
    std::cout << message + "\n";
    std::cout << "                                   " << minimum << "-" << maximum << ": ";
    std::cout << arguments[0] << "\n";
    std::string arg = arguments.at(0);
    arguments.erase(arguments.begin());
    double param;

    std::stringstream stream(arg);
    if (stream >> param)
    {
      if (param < minimum || param > maximum)
      {
        std::cout << "Argument '" << message << " = " << arg
                  << "' out of bounds. Minimum: "
                  << minimum << "Maximum: " << maximum << "\n";
        return PromptDouble(message, name, minimum, maximum);
      }
      else
        return param;
    }
    else
    {
      std::cout << "Argument '" << message << " = " << arg << "' invalid.\n";
      return PromptDouble(message, name, minimum);
    }
  }
}

void ConvertCharPointerArrayToStringVector(char *c[], int arrayLength, std::vector<std::string> &s)
{
  if (arrayLength == 0)
    return;
  for (int i = 0; i < arrayLength; ++i)
  {
    std::string str = std::string(c[i]);
    s.insert(s.end(), str);
  }
}