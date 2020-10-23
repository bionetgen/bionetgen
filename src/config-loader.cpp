#include "./voronoi.cpp"

#include <boost/bind/bind.hpp>

class ConfigLoader
{
    boost::property_tree::ptree _root;

public:
    void load(std::string path)
    {
        // Load the json file in this ptree
        boost::property_tree::read_json(path, this->_root);
    }
    void configure(Voronoi& voronoi)
    {
        voronoi.configure(this->_root);
    }
};
