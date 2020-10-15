/*!----------------------------------------------------------------------
 \file lib_elements.hpp
 \brief element related data container

 \maintainer Jonas Koban

 \level 3

 */
/*-----------------------------------------------------------*/

#ifndef LIB_ELEMENTS_HPP_
#define LIB_ELEMENTS_HPP_

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <sstream>

namespace FilamentNetwork
{

namespace Elements
{
enum ElementType
{
  none,
  rigid_sphere
};
struct Element
{
public:
  Element(){};
  virtual ~Element() = 0;

protected:
  //old
  //ElementType _type;

public:
  virtual ElementType GetType() { return ElementType::none; };
};
inline Element::~Element(){};

struct RigidSphereElement : Element
{
public:
  RigidSphereElement() : Element(){};
  RigidSphereElement(std::vector<double> position, double radius, double density);
  virtual ~RigidSphereElement(){};

private:
  //member hiding! this is WRONG!
  //ElementType _type
  double _radius;
  std::vector<double> _position;
  double _density;
  bool fixedByDBC_ = false;

public:
  ElementType GetType() override { return ElementType::rigid_sphere; };

  void SetPosition(std::vector<double> position);
  std::vector<double> GetPosition();
  void SetRadius(double radius);
  double GetRadius();
  void SetDensity(double density);
  double GetDensity();
  void SetFixedByDBC(bool fixedByDBC = true);
  bool GetFixedByDBC();
};
}
}
#endif /* LIB_ELEMENTS_HPP_ */
