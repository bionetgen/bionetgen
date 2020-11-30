#include "./voronoi.cpp"

#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/property_tree/json_parser.hpp>
class ConfigLoader
{
    boost::property_tree::ptree _root;
    boost::filesystem::path _config_path;

public:
    void load(std::string path)
    {
        this->_config_path = boost::filesystem::path(path).parent_path();
        boost::property_tree::read_json(path, this->_root);
    }
    void configure(Voronoi& voronoi)
    {
        voronoi.configure(this->_config_path, this->_root);
    }
};
