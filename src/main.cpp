#include "./voronoi.cpp"

int main(int argc, char const *argv[])
{
    Voronoi voronoi = Voronoi();
    FilamentNetworkProblem *filanetprob = new FilamentNetworkProblem();
    filanetprob->boxsize_ = std::vector<double>{1, 1, 1};
    filanetprob->boxZeroPoint_ = std::vector<double>{0, 0, 0};
    voronoi.Initialize();
    voronoi.ComputeVoronoi(filanetprob);
    return 0;
}
