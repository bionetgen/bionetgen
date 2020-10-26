#include "./config-loader.cpp"

int main(int argc, char const *argv[])
{
    Voronoi voronoi = Voronoi();
    FilamentNetworkProblem *filanetprob = new FilamentNetworkProblem();
    ConfigLoader loader = ConfigLoader();
    std::string configFilePath = argv[1];
    loader.load(configFilePath);
    loader.configure(voronoi);
    voronoi.run(filanetprob);
    return 0;
}
