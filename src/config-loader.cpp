#include "./voronoi.cpp"

#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/property_tree/json_parser.hpp>
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
