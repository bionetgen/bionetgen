/*!----------------------------------------------------------------------
\file filanetwork_utils.hpp
\brief contains some helper methods for filament network generation

\maintainer Jonas Eichinger

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FILANETWORK_UTILS_H_
#define FILANETWORK_UTILS_H_

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <sstream>

// read double parameter
void ReadDoubleParam(double& param, std::string& name, double minimum = 1e-16);
// read integer parameter
void ReadIntParam(int& param, std::string& name, bool includezero=false);

bool AskYesNo(std::string message);
int Ask10(std::string message);
int AskInt(std::string message, std::string& name, bool includezero);
double AskDouble(std::string message, std::string& name, double minimum = 1e-16);

#endif /* FILANETWORK_UTILS_H_ */
