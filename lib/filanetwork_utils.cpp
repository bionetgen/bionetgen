/*!----------------------------------------------------------------------
 \file filanetwork_utils.cpp
 \brief contains some helper methods for filament network generation

 \maintainer Jonas Eichinger

 \level 3

 */
/*-----------------------------------------------------------*/

#include "filanetwork_utils.hpp"
#include <vector>
#include <cmath>
/*----------------------------------------------------------------------*
 | Read double parameter                         (protected) mueller 01/12|
 *----------------------------------------------------------------------*/
void ReadDoubleParam(double& param, std::string& name, double minimum)
{
  std::string paramstring;
  while (true)
  {
    getline(std::cin, paramstring);
    // This code converts from string to number safely.
    std::stringstream myStream(paramstring);
    if(myStream >> param)
    {
      if(param <= minimum)
      {
        std::cout << "Invalid input, " << name << "=" << param << "<"
            << minimum << "! Please reenter value: ";
        continue;
      }
      break;
    }
    else
      std::cout << "Invalid input, " << name << "=" << paramstring
          << "! Please reenter value: ";
  }
  return;
} //ReadDoubleParam

/*----------------------------------------------------------------------*
 | Read integer parameter                  (public,static) mueller 01/12|
 *----------------------------------------------------------------------*/
void ReadIntParam(int& param, std::string& name, bool includezero)
{
  int min = 0;
  if(includezero)
    min = -1;

  std::string paramstring;
  while (true)
  {
    getline(std::cin, paramstring);
    // This code converts from string to number safely.
    std::stringstream myStream(paramstring);
    if(myStream >> param)
    {
      if(param <= min)
      {
        std::cout << "Invalid input, " << name << "=" << param << "<"
            << min << "! Please reenter value: ";
        continue;
      }
      break;
    }
    else
      std::cout << "Invalid input, " << name << "=" << paramstring
          << "! Please reenter value: ";
  }
  return;
} // ReadIntParam()
bool AskYesNo(std::string message)
{
  std::cout << message + "\n";
  std::string readString;
  while (true)
  {
    std::cout << "                             1/y(es)   0/n(o): ";
    getline(std::cin, readString);
    if(readString == "y" || readString == "1")
      return true;
    else if(readString == "n" || readString == "0")
      return false;
    else
      std::cout
      << "Input \'" + readString + "\' invalid. "
          + "Try again!\n";
  }
} // AskYesNo()
int Ask10(std::string message)
{
  std::cout << message + "\n";
  std::string readString;
  while (true)
  {
    std::cout << "                             1/y(es)   0/n(o): ";
    getline(std::cin, readString);
    if(readString == "y" || readString == "1")
      return 1;
    else if(readString == "n" || readString == "0")
      return 0;
    else
      std::cout
      << "Input \'" + readString + "\' invalid. "
          + "Try again!\n";
  }
} // AskYesNo()
int AskInt(std::string message, std::string& name, bool includezero)
{
  std::cout << message + "\n";
  int param;
  int min = 1;
  if(includezero)
    min = 0;

  std::string readString;
  while (true)
  {
    std::cout << "                                         " << min << "-...: ";
    getline(std::cin, readString);
    // This code converts from string to number safely.
    std::stringstream stream(readString);
    if(stream >> param)
    {
      if(param < min)
      {
        std::cout << "Invalid input, " << name << "=" << param << "<"
            << min << "! Please reenter value: ";
        continue;
      }
      break;
    }
    else
      std::cout << "Invalid input, " << name << "=" << readString
          << "! Please reenter value: ";
  }
  return param;
} // AskInt()
double AskDouble(std::string message, std::string& name, double minimum)
{
  std::cout << message + "\n";
  double param;
  std::string stringRead;
    while (true)
    {
      std::cout << "                                   " << minimum << "-...: ";
      getline(std::cin, stringRead);
      // This code converts from string to number safely.
      std::stringstream stream(stringRead);
      if(stream >> param)
      {
        if(param <= minimum)
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
