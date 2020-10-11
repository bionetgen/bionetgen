#include "./voronoi.cpp"

int main(int argc, char const *argv[])
{
    Voronoi voronoi = Voronoi();
    FilamentNetworkProblem *filanetprob = new FilamentNetworkProblem();
    voronoi.ComputeVoronoi(filanetprob);
    return 0;
}
