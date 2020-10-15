/*!----------------------------------------------------------------------
 \file lib_prompt.hpp
 \brief Prompt functions for input_generator

 \maintainer Jonas Koban

 \level 3

 */
/*-----------------------------------------------------------*/

#ifndef LIB_PROMPT_H_
#define LIB_PROMPT_H_

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <sstream>

bool PromptYesNo(std::string message);
bool PromptYesNoCMD(std::string message, std::vector<std::string> &arguments);
unsigned int Prompt10(std::string message);
unsigned int Prompt10CMD(std::string message, std::vector<std::string> &arguments);
int PromptInt(std::string message, std::string name, bool includeZero);
int PromptIntCMD(std::string message, std::string name, std::vector<std::string> &arguments, bool includeZero);
unsigned int PromptUInt(std::string message, std::string name, bool includeZero);
unsigned int PromptUIntCMD(std::string message, std::string name, std::vector<std::string> &arguments, bool includeZero);
double PromptDouble(std::string message, std::string name, double minimum = 1e-16, double maximum = 1e16);
double PromptDoubleCMD(std::string message, std::string name, std::vector<std::string> &arguments, double minimum = 1e-16, double maximum = 1e16);

void ConvertCharPointerArrayToStringVector(char *c[], int arrayLength, std::vector<std::string> &s);

#endif /* LIB_PROMPT_H_ */
