/*!----------------------------------------------------------------------
 \file lib_elements.cpp
 \brief element related data container

 \maintainer Jonas Koban

 \level 3
 */
/*-----------------------------------------------------------*/
#include "lib_elements.hpp"
#include <cstdio>
#include <vector>
#include <string>
namespace FilamentNetwork
{
namespace Elements
{

RigidSphereElement::RigidSphereElement(std::vector<double> position, double radius, double density) : Element()
{
  this->SetPosition(position);
  this->SetRadius(radius);
  this->SetDensity(density);
}
void RigidSphereElement::SetPosition(std::vector<double> position)
{
  this->_position = position;
}
std::vector<double> RigidSphereElement::GetPosition()
{
  return this->_position;
}
void RigidSphereElement::SetRadius(double radius)
{
  this->_radius = radius;
}
double RigidSphereElement::GetRadius()
{
  return this->_radius;
}
void RigidSphereElement::SetDensity(double density)
{
  this->_density = density;
}
double RigidSphereElement::GetDensity()
{
  return this->_density;
}
void RigidSphereElement::SetFixedByDBC(bool fixedByDBC)
{
  this->fixedByDBC_ = fixedByDBC;
}
bool RigidSphereElement::GetFixedByDBC()
{
  return this->fixedByDBC_;
}
}
}
