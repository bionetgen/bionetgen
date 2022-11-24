/* _________________________________________________________________________________
 *
 * BSD 3-Clause License
 *
 * Copyright (c) 2021, bionetgen
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * _________________________________________________________________________________
 */


#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <random>
#include <chrono>
#include <map>
#include <set>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>

#include "../voro++-0.4.6/src/voro++.hh"
// #include "../lib/lib_vec.hpp"

const double one_third = 1.0 / 3.0;
using uint = unsigned;

class Voronoi
{
  bool generate_;
  bool simulate_;
  uint mode_;
  std::filesystem::path outputPrefix_;
  std::filesystem::path inputPrefix_;
  std::vector<double> boxSize_;
  std::vector<double> boxOrigin_;
  double seed_;
  uint voronoiParticleCount_;
  uint currnumfils_;
  // position of center of particles
  std::vector<std::vector<double>> particlePositions_ = std::vector<std::vector<double>>();
  // contains vertex id with its position
  std::vector<std::vector<double>> vertices_ = std::vector<std::vector<double>>();
  // contains edges and the respective nodes it is attached to
  std::vector<std::vector<unsigned int>> edges_ = std::vector<std::vector<unsigned int>>();
  // contains order of a vertices
  std::vector<unsigned int> vertexEdgeCount_ = std::vector<unsigned int>();
  // nodes to edges (contains dead nodes)
  std::vector<std::vector<unsigned int>> node_to_edges_;
  // the following map does not contain dead nodes, i.e. nodes with z == 0
  // so in here only really existing nodes
  std::map<unsigned int, std::vector<double>> vertices_map_;
  // all really existing nodes with their respective edges
  std::map<unsigned int, std::vector<unsigned int>> edge_map_;
  // a vector containing all vertex ids of really existing verteces. Note
  // that index of vector is not equal to to node id
  std::vector<unsigned int> vertices_for_random_draw_;
  // shifted vertex positions
  std::vector<std::vector<double>> vtxs_shifted_;
  // contains all finite element node ids that belong to vertex
  std::vector<std::vector<unsigned int>> vertexNodeIds_ = std::vector<std::vector<unsigned int>>();
  // Input parameters for Simulated Annealing
  unsigned int max_iter;
  unsigned int max_subiter;
  double weight_line;
  double weight_cosine;
  double tolerance;
  double temperature_inital;
  double decay_rate_temperature;
  double max_movement;
  unsigned int screen_output_every;
  // for binning
  unsigned int p_num_bins_lengths;
  unsigned int p_num_bins_cosines;

public:
  void configure(std::filesystem::path config_path, boost::property_tree::ptree config)
  {
    this->outputPrefix_ = std::filesystem::path(config_path) / std::filesystem::path(config.get<std::string>("output-prefix"));
    if (this->outputPrefix_.has_parent_path())
      if (!std::filesystem::exists(this->outputPrefix_.parent_path()))
        std::filesystem::create_directory(this->outputPrefix_.parent_path());
    this->inputPrefix_ = std::filesystem::path(config_path) / std::filesystem::path(config.get<std::string>("input-prefix"));
    if (this->inputPrefix_.has_parent_path())
      if (!std::filesystem::exists(this->inputPrefix_.parent_path()))
        std::filesystem::create_directory(this->inputPrefix_.parent_path());

    this->seed_ = config.get<double>("seed");
    this->voronoiParticleCount_ = config.get<uint>("particles");

    this->generate_ = config.get<bool>("generate");
    this->simulate_ = config.get<bool>("simulate");

    for (auto child : config.get_child("box-size"))
      this->boxSize_.push_back(child.second.get_value<double>());
    for (auto child : config.get_child("box-origin"))
      this->boxOrigin_.push_back(child.second.get_value<double>());

    // Simulated Annealing
    auto config_sa = config.get_child("simulated-annealing");
    std::string mode = config_sa.get<std::string>("mode");
    if (mode == "1")
      this->mode_ = 1;
    else if (mode == "2")
      this->mode_ = 2;
    else if (mode == "both")
      this->mode_ = 3;
    else
      throw "Invalid mode. Can only be '1', '2' or 'both'.";
    this->max_iter = config_sa.get<uint>("max-iter");
    this->max_subiter = config_sa.get<uint>("max-subiter");
    this->weight_line = config_sa.get<double>("weight-line");
    this->weight_cosine = config_sa.get<double>("weight-cosine");
    this->tolerance = config_sa.get<double>("tolerance");
    this->temperature_inital = config_sa.get<double>("temperature-initial");
    this->decay_rate_temperature = config_sa.get<double>("temperature-decay-rate");
    const double max_movementFrac = config_sa.get<double>("max-movement-frac");
    this->max_movement = max_movementFrac * this->boxSize_[0];
    this->screen_output_every = config_sa.get<uint>("screen-output-every");

    // for binning
    this->p_num_bins_lengths = config_sa.get<uint>("num-bins-length");
    this->p_num_bins_cosines = config_sa.get<uint>("num-bins-cosine");
  }

  void run()
  {
    std::random_device rd;
    std::mt19937 gen(rd());
    gen.seed(this->seed_);
    // random number between 0 and 1
    std::uniform_real_distribution<> dis_uni(0, 1);
    if (this->generate_)
      this->ComputeVoronoi(gen, dis_uni);
    else
      this->ReadGeometry();

    gen.seed(this->seed_);
    if (this->simulate_)
      this->SimulatedAnnealing(this->mode_, gen, dis_uni);

    this->OutputGeometry();
  }

  void ComputeVoronoi(std::mt19937 &gen, std::uniform_real_distribution<> &dis_uni)
  {
    std::cout << "\n\nNetwork Generation started." << std::endl;
    std::cout << "------------------------------------------------------\n"
              << std::flush;
    std::cout << "1) Computing Voronoi.\n"
              << std::flush;
    auto start_voro = std::chrono::high_resolution_clock::now();

    // reset variables in case of unsuccessful computation
    this->particlePositions_.clear();
    this->vertices_.clear();
    this->edges_.clear();
    this->vertexEdgeCount_.clear();

    double compareTolerance = 1e-13;
    std::string nullString = "";

    double x_min;
    double x_max;
    double y_min;
    double y_max;
    double z_min;
    double z_max;

    x_min = this->boxOrigin_[0] - this->boxSize_[0] / 2;
    x_max = this->boxOrigin_[0] + this->boxSize_[0] / 2;
    y_min = this->boxOrigin_[1] - this->boxSize_[1] / 2;
    y_max = this->boxOrigin_[1] + this->boxSize_[1] / 2;
    z_min = this->boxOrigin_[2] - this->boxSize_[2] / 2;
    z_max = this->boxOrigin_[2] + this->boxSize_[2] / 2;
    int n_x = 1, n_y = 1, n_z = 1;

    unsigned int particleCount = this->voronoiParticleCount_;
    // allocate
    // this->cellVertexOrders_ = std::vector<std::vector<uint32_t>>(particleCount);
    // this->cellRads_ = std::vector<std::vector<double>>(particleCount);
    // this->cellRadAvg_ = std::vector<double>(particleCount);
    // this->cellRadMax_ = std::vector<double>(particleCount);
    // this->cellRadMin_ = std::vector<double>(particleCount);
    // this->cellRad_in_ = std::vector<double>(particleCount);

    double x, y, z;

    // Create a container with the geometry given above. Allocate space for
    // eight particles within each computational block
    voro::container con(x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z,
                        true,
                        true,
                        true, 8);

    this->particlePositions_.reserve(particleCount);

    for (unsigned int i = 0; i < particleCount; ++i)
    {
      x = x_min + dis_uni(gen) * (x_max - x_min);
      y = y_min + dis_uni(gen) * (y_max - y_min);
      z = z_min + dis_uni(gen) * (z_max - z_min);
      con.put(i, x, y, z);
      this->particlePositions_.emplace_back(std::vector<double>{x, y, z});
    }

    voro::c_loop_all loop = voro::c_loop_all(con);

    // cell index = particle index
    unsigned int cellIndex = 0;
    unsigned int cellCount = particleCount;
    std::vector<double> vertices = std::vector<double>();

    int debug_lastUnique = -1;
    unsigned int count = 1;
    if (loop.start())
      do
      {
        voro::voronoicell_neighbor cell;
        if (con.compute_cell(cell, loop))
        {
          // unused
          int cellId;
          double x, y, z, radius;
          loop.pos(cellId, x, y, z, radius);

          cell.vertices(x, y, z, vertices);

          std::vector<std::vector<double>> cellVertices = std::vector<std::vector<double>>(vertices.size() / 3);

          for (unsigned int vertexIndex = 0; vertexIndex < vertices.size(); vertexIndex += 3)
          {
            cellVertices[vertexIndex / 3] = std::vector<double>{vertices[vertexIndex],
                                                                vertices[vertexIndex + 1],
                                                                vertices[vertexIndex + 2]};
          }

          unsigned int vertexCount = cellVertices.size();

          // statistical data for network characterization
          std::vector<double> cellRads = std::vector<double>();
          double maxRad = 0;
          double minRad = 0;
          double avgRad = 0;

          std::vector<int> vertexEdgeCount;
          cell.vertex_orders(vertexEdgeCount);

          std::vector<std::vector<int>> vertexPartners = std::vector<std::vector<int>>(vertexCount);
          // std::vector<std::vector<uint32_t>> cellEdges = std::vector<std::vector<uint32_t>>();

          unsigned int partner1Index;
          unsigned int partner2Index;

          for (int vertexIndex = 0; vertexIndex < vertexCount; ++vertexIndex)
          {
            // edges
            vertexPartners[vertexIndex] = std::vector<int>(vertexEdgeCount[vertexIndex]);

            for (int partnerIndex = 0; partnerIndex < vertexEdgeCount[vertexIndex]; ++partnerIndex)
              vertexPartners[vertexIndex][partnerIndex] = cell.ed[vertexIndex][partnerIndex];

            // statistical
            double radius =
                std::sqrt((cellVertices[vertexIndex][0] - this->particlePositions_[cellIndex][0]) *
                              (cellVertices[vertexIndex][0] - this->particlePositions_[cellIndex][0]) +
                          (cellVertices[vertexIndex][1] - this->particlePositions_[cellIndex][1]) *
                              (cellVertices[vertexIndex][1] - this->particlePositions_[cellIndex][1]) +
                          (cellVertices[vertexIndex][2] - this->particlePositions_[cellIndex][2]) *
                              (cellVertices[vertexIndex][2] - this->particlePositions_[cellIndex][2]));
            cellRads.push_back(radius);
            avgRad += radius;
            if (radius > maxRad)
              maxRad = radius;
            if (radius < minRad | minRad == 0)
              minRad = radius;

            // actual network creation from output of voro++
            std::vector<double> vertexPositionCurrent = cellVertices[vertexIndex];
            bool vertexIsUnique = true;

            for (int vertexIndex = 0; vertexIndex < this->vertices_.size(); ++vertexIndex)
            {
              if (std::abs(this->vertices_[vertexIndex][0] - vertexPositionCurrent[0]) < compareTolerance &&
                  std::abs(this->vertices_[vertexIndex][1] - vertexPositionCurrent[1]) < compareTolerance &&
                  std::abs(this->vertices_[vertexIndex][2] - vertexPositionCurrent[2]) < compareTolerance)
              {
                vertexIsUnique = false;
                partner1Index = vertexIndex;
                // std::cout << "1: " <<
                // std::to_string(partner1Index) << "\n";
              }
            }

            if (vertexIsUnique)
            {
              partner1Index = this->vertices_.size();
              this->vertices_.push_back(vertexPositionCurrent);

              if (debug_lastUnique > 0 && partner1Index < debug_lastUnique)
              {
                std::cout << "1U: " << std::to_string(partner1Index)
                          << "\n";
                std::cout << "^ ERROR! Unique vertex id is not higher than last! "
                             "Press any key to continue...\n";
                getline(std::cin, nullString);
              }
              debug_lastUnique++;
            }

            std::vector<bool> onPlanesPartner1 = std::vector<bool>{
                this->VertexIsOnHighPlane(partner1Index, 0) || this->VertexIsOnLowPlane(partner1Index, 0),
                this->VertexIsOnHighPlane(partner1Index, 1) || this->VertexIsOnLowPlane(partner1Index, 1),
                this->VertexIsOnHighPlane(partner1Index, 2) || this->VertexIsOnLowPlane(partner1Index, 2)};

            unsigned int edgesCreated = 0;
            for (int partnerIndex = 0; partnerIndex < vertexEdgeCount[vertexIndex]; ++partnerIndex)
            {
              int vertexPartnerIndexCurrent = vertexPartners[vertexIndex][partnerIndex];
              std::vector<double> vertexPartnerPositionCurrent = cellVertices[vertexPartnerIndexCurrent];

              bool vertexPartnerIsUnique = true;
              for (int vertexIndex = 0; vertexIndex < this->vertices_.size(); ++vertexIndex)
              {
                if (std::abs(this->vertices_[vertexIndex][0] - vertexPartnerPositionCurrent[0]) < compareTolerance &&
                    std::abs(this->vertices_[vertexIndex][1] - vertexPartnerPositionCurrent[1]) < compareTolerance &&
                    std::abs(this->vertices_[vertexIndex][2] - vertexPartnerPositionCurrent[2]) < compareTolerance)
                {
                  vertexPartnerIsUnique = false;
                  partner2Index = vertexIndex;
                  // std::cout << "2: " <<
                  // std::to_string(partner2Index) << "\n";
                }
              }

              if (vertexPartnerIsUnique)
              {
                partner2Index = this->vertices_.size();
                this->vertices_.push_back(vertexPartnerPositionCurrent);

                if (debug_lastUnique > 0 && partner2Index < debug_lastUnique)
                {
                  std::cout << "2U: "
                            << std::to_string(partner2Index)
                            << "\n";
                  std::cout << "^ ERROR! Unique vertex id is not higher than "
                               "last! Press any key to continue...\n";
                  getline(std::cin, nullString);
                }
                debug_lastUnique++;
              }

              if (partner1Index == partner2Index)
              {
                std::cout << "Error: VertexUIds equal. This should not happen!\n";
                continue;
              }

              bool edgeIsUnique = true;
              for (int edgeIndex = 0; edgeIndex < this->edges_.size(); ++edgeIndex)
                if ((this->edges_[edgeIndex][0] == partner1Index &&
                     this->edges_[edgeIndex][1] == partner2Index) ||
                    (this->edges_[edgeIndex][0] == partner2Index &&
                     this->edges_[edgeIndex][1] == partner1Index))
                {
                  edgeIsUnique = false;
                }

              if (edgeIsUnique)
              {
                //! OLD
                std::vector<bool> onPlanesPartner2 = std::vector<bool>{
                    this->VertexIsOnHighPlane(partner2Index, 0) || this->VertexIsOnLowPlane(partner2Index, 0),
                    this->VertexIsOnHighPlane(partner2Index, 1) || this->VertexIsOnLowPlane(partner2Index, 1),
                    this->VertexIsOnHighPlane(partner2Index, 2) || this->VertexIsOnLowPlane(partner2Index, 2)};

                // 1 1 1 1
                if ((!true && onPlanesPartner1[0] && onPlanesPartner2[0]) ||
                    (!true && onPlanesPartner1[1] && onPlanesPartner2[1]) ||
                    (!true && onPlanesPartner1[2] && onPlanesPartner2[2])) // if(onPlanesPartner1 >= (1) &&
                                                                           // vertexIsUnique && edgesCreated ==
                                                                           // 0)
                {
                  if (vertexPartnerIsUnique)
                  {
                    // remove 2
                    this->vertices_.erase(this->vertices_.begin() + partner2Index);
                    debug_lastUnique--;

                    // WARNING
                    if (partner2Index < partner1Index)
                      --partner1Index;
                  }
                  continue;
                }
                else
                {
                  this->edges_.push_back(std::vector<unsigned int>{partner1Index,
                                                                   partner2Index});
                  ++edgesCreated;
                }
              }

              // cellEdges.push_back(std::vector<unsigned int>{
              //     partner1Index, partner2Index});
            }

            if (((!true && onPlanesPartner1[0]) ||
                 (!true && onPlanesPartner1[1]) ||
                 (!true && onPlanesPartner1[2])) &&
                vertexIsUnique && edgesCreated == 0) // no partner2 depends on partner1?
            {
              // remove 1
              this->vertices_.erase(this->vertices_.begin() + partner1Index);
              --debug_lastUnique;
            }
          }

          // finalize statistical data for network characterization
          // for (uint32_t vtxOrder : vertexEdgeCount_)
          //   this->cellVertexOrders_[cellIndex].push_back(vtxOrder);
          // double cellRad_in = 0;
          // size_t proj_num = 0;
          // for (std::vector<uint32_t> edge : edges)
          // {
          //   vec3 vtxP1_pos(this->uniqueVertices_[edge[0]]);
          //   vec3 vtxP2_pos(this->uniqueVertices_[edge[1]]);
          //   vec3 dir = vtxP2_pos - vtxP1_pos;
          //   double dir_abs = dir.length();
          //   //* projection *
          //   //           part
          //   //     loc  . `
          //   //       . `  | dist
          //   //     . ` cos |
          //   // p1 ---------------------------- p2
          //   //            x   dir
          //   vec3 part_pos(this->particlePositions_[cellIndex]);
          //   vec3 loc = part_pos - vtxP1_pos;

          //   double cos_abs = loc.dot(dir);
          //   vec3 cos = dir * (cos_abs / dir_abs);
          //   vec3 x = vtxP1_pos + cos;
          //   vec3 dist = part_pos - x;
          //   double dist_abs = dist.length();
          //   if ((cos_abs > 0) & (cos_abs < dir_abs)) // projection is inside
          //   {
          //     cellRad_in += dist_abs;
          //     ++proj_num;
          //   }
          // }
          // if (proj_num)
          //   cellRad_in /= proj_num;
          // else
          //   std::cout << " ! | Projections rejected. Outside of all "
          //                "edges.\nParticle pos: "
          //             << this->particlePositions_[cellIndex][0] << " "
          //             << this->particlePositions_[cellIndex][1] << " "
          //             << this->particlePositions_[cellIndex][2] << "\n";
          // this->cellRad_in_[cellIndex] = cellRad_in;
          // avgRad /= vertexCount;
          // this->cellRads_[cellIndex] = cellRads;
          // this->cellRadAvg_[cellIndex] = avgRad;
          // this->cellRadMax_[cellIndex] = maxRad;
          // this->cellRadMin_[cellIndex] = minRad;

          if (cellIndex >= 0.05 * count * cellCount - 1)
          {
            auto stop_voro = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed_voro = stop_voro - start_voro;
            std::cout << "   " << elapsed_voro.count() / 60 << " minutes needed for " << std::flush;
            std::cout << "   " << static_cast<int>(0.05 * count * 100.0) << " %\r" << std::flush;
            ++count;
          }

          ++cellIndex;
        }
      } while (loop.inc());

    std::cout << std::endl;

    //! NEW shifting to add all filaments that are not present on both periodic
    //! boundaries
    // this->ShiftEdges(this->uniqueVertexEdgePartners_,
    // this->uniqueVertices_);
    //! this breaks calc of fiber vol frac

    auto stop_voro = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_voro = stop_voro - start_voro;
    std::cout << "\n   Computation of Voronoi is done now. It took " << elapsed_voro.count() / 60 << " minutes. \n";

    // remove double edges in periodic BC dimension
    std::cout << "\n2) Removing double edges. Current edge count: " << this->edges_.size() << "\n";
    auto start_removing_doubles = std::chrono::high_resolution_clock::now();
    // now clean up
    node_to_edges_ = std::vector<std::vector<unsigned int>>(vertices_.size(), std::vector<unsigned int>());
    for (unsigned int i_edge = 0; i_edge < edges_.size(); ++i_edge)
    {
      node_to_edges_[this->edges_[i_edge][0]].push_back(i_edge);
      node_to_edges_[this->edges_[i_edge][1]].push_back(i_edge);
    }

    //! Careful: creates dead nodes! (vertices still exist!)
    // if (this->applyPeriodicBCsPerDim_[0] or this->applyPeriodicBCsPerDim_[1] or this->applyPeriodicBCsPerDim_[2])
    this->RemoveDoubles();

    auto stop_removing_doubles = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_removing_doubles = stop_removing_doubles - start_removing_doubles;
    std::cout << "   Double edges are removed now. It took " << elapsed_removing_doubles.count() / 60 << " minutes. \n";

    std::cout << "\n3) Cleaning up ... " << std::endl;
    ;
    auto start_cleaning = std::chrono::high_resolution_clock::now();

    // now clean up
    node_to_edges_ = std::vector<std::vector<unsigned int>>(vertices_.size(), std::vector<unsigned int>());
    for (unsigned int i_edge = 0; i_edge < edges_.size(); ++i_edge)
    {
      node_to_edges_[this->edges_[i_edge][0]].push_back(i_edge);
      node_to_edges_[this->edges_[i_edge][1]].push_back(i_edge);
    }

    vertices_map_.clear();
    this->vtxs_shifted_ = this->vertices_;
    this->ShiftVertices(this->vtxs_shifted_);

    std::vector<unsigned int> vertices_with_order_1;
    std::vector<unsigned int> vertices_with_order_3;
    for (unsigned int i_node = 0; i_node < node_to_edges_.size(); ++i_node)
    {
      // if (node_to_edges_[i_node].size() == 4)
      //   vertices_map_[i_node] = vtxs_shifted_[i_node];

      if (node_to_edges_[i_node].size() == 1)
        vertices_with_order_1.push_back(i_node);

      if (node_to_edges_[i_node].size() == 3)
        vertices_with_order_3.push_back(i_node);
    }

    // find each partner
    for (unsigned int i = 0; i < vertices_with_order_3.size(); ++i)
    {
      for (unsigned int j = 0; j < vertices_with_order_1.size(); ++j)
      {
        if (VerticesMatch(vtxs_shifted_[vertices_with_order_3[i]],
                          vtxs_shifted_[vertices_with_order_1[j]]))
        {
          if (edges_[node_to_edges_[vertices_with_order_1[j]][0]][0] == vertices_with_order_1[j])
            edges_[node_to_edges_[vertices_with_order_1[j]][0]][0] = vertices_with_order_3[i];
          else
            edges_[node_to_edges_[vertices_with_order_1[j]][0]][1] = vertices_with_order_3[i];

          // vertices_map_[vertices_with_order_3[i]] = vtxs_shifted_[vertices_with_order_3[i]];
        }
      }
    }

    //! remove 'dead' vertices
    std::vector<bool> removeVertices = std::vector<bool>(vertices_.size(), true);
    auto edgeId = 0;
    for (auto &edge : this->edges_)
    {
      removeVertices[edge[0]] = false;
      removeVertices[edge[1]] = false;
    }

    auto finalSize = std::count_if(removeVertices.begin(), removeVertices.end(), [](bool remove) { return !remove; });
    std::vector<uint> newVertexIds = std::vector<uint>(vertices_.size());
    auto newVertices = std::vector<std::vector<double>>(finalSize);
    auto newVertexId = 0;
    for (uint vertexId = 0; vertexId < vertices_.size(); vertexId++)
    {
      if (!removeVertices[vertexId])
      {
        newVertexIds[vertexId] = newVertexId;
        newVertices[newVertexId] = vertices_[vertexId];
        vertices_map_[newVertexId] = vtxs_shifted_[vertexId];
        newVertexId++;
      }
    }
    vertices_ = newVertices;

    for (auto &edge : this->edges_)
    {
      edge[0] = newVertexIds[edge[0]];
      edge[1] = newVertexIds[edge[1]];
    }

    node_to_edges_ = std::vector<std::vector<unsigned int>>(vertices_.size(), std::vector<unsigned int>());
    for (unsigned int i_edge = 0; i_edge < edges_.size(); ++i_edge)
    {
      node_to_edges_[this->edges_[i_edge][0]].push_back(i_edge);
      node_to_edges_[this->edges_[i_edge][1]].push_back(i_edge);
    }

    auto stop_cleaning = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapses_cleaning = stop_cleaning - start_cleaning;
    std::cout << "   Cleaning is done now. It took " << elapses_cleaning.count() / 60 << " minutes. \n";

    std::cout << "   Number of lines: " << this->edges_.size() << "\n"
              << std::flush;
    std::cout << "   Number of nodes: " << this->vertices_for_random_draw_.size() << "\n"
              << std::flush;

    bool adapt_connectivity = true;
    if (adapt_connectivity)
      this->AdaptConnectivity(gen, dis_uni);

    // Compute Vertex order
    this->vertexEdgeCount_.clear();
    for (unsigned int i_node = 0; i_node < node_to_edges_.size(); ++i_node)
    {
      this->vertexEdgeCount_.push_back(node_to_edges_[i_node].size());
    }
    // std::cout << "   Average vertex order is " << CalcAvgVtxOrder() << std::endl;

    // statistical network data
    // this->networkCellRadAvg_ = 0;
    // this->networkCellRadAvg_in_ = 0;
    // this->networkCellRadMax_ = 0;
    // this->networkCellRadMin_ = 0;
    // for (unsigned int cellIndex = 0; cellIndex < cellCount; ++cellIndex)
    // {
    //   this->networkCellRadAvg_ += this->cellRadAvg_[cellIndex];
    //   this->networkCellRadAvg_in_ += this->cellRad_in_[cellIndex];
    //   this->networkCellRadMax_ += this->cellRadMax_[cellIndex];
    //   this->networkCellRadMin_ += this->cellRadMin_[cellIndex];
    // }

    // this->networkCellRadAvg_ /= cellCount;
    // this->networkCellRadAvg_in_ /= cellCount;
    // this->networkCellRadMax_ /= cellCount;
    // this->networkCellRadMin_ /= cellCount;

    // this->networkCellRadStdDev_ = 0;
    // double radiusDeltaSquared = 0;
    // for (unsigned int cellIndex = 0; cellIndex < cellCount; ++cellIndex)
    // {
    //   radiusDeltaSquared +=
    //       (this->cellRadAvg_[cellIndex] - this->networkCellRadAvg_) *
    //       (this->cellRadAvg_[cellIndex] - this->networkCellRadAvg_);
    // }

    // this->networkCellRadStdDev_ = std::sqrt(radiusDeltaSquared / (cellCount - 1));
    // this->networkCellRadStdDevNorm_ = this->networkCellRadStdDev_ / this->networkCellRadAvg_;

    // allocate memory for dnodes (finite element nodes of actual discretization)
    this->vertexNodeIds_ = std::vector<std::vector<unsigned int>>(this->vertices_.size());

    int numberOfFilaments = this->edges_.size();
    this->currnumfils_ = numberOfFilaments;
  }

  double GetFilamentLength(
      unsigned int filamentIndex) const
  {
    return this->GetEdgeLength(filamentIndex);
  }

  double GetEdgeLength(unsigned int edgeUId) const
  {
    std::vector<unsigned int> partners = this->edges_[edgeUId];
    std::vector<double> partner1Position = this->vertices_[partners[0]];
    std::vector<double> partner2Position = this->vertices_[partners[1]];

    UnShift3D(partner1Position, partner2Position);

    double deltaX = partner2Position[0] - partner1Position[0];
    double deltaY = partner2Position[1] - partner1Position[1];
    double deltaZ = partner2Position[2] - partner1Position[2];
    double squareX = deltaX * deltaX;
    double squareY = deltaY * deltaY;
    double squareZ = deltaZ * deltaZ;
    return std::sqrt(squareX + squareY + squareZ);
  }

  bool VertexIsOnHighPlane(
      unsigned int vertexUId,
      unsigned int dimension) const
  {
    return this->PointIsOnHighPlane(this->vertices_[vertexUId],
                                    dimension);
  }

  bool VertexIsOnLowPlane(
      unsigned int vertexUId,
      unsigned int dimension) const
  {
    return this->PointIsOnLowPlane(this->vertices_[vertexUId],
                                   dimension);
  }

  bool PointIsOverHighPlane(
      std::vector<double> const &point,
      unsigned int dimension) const
  {
    return (this->boxOrigin_[dimension] +
            this->boxSize_[dimension] / 2) <= point[dimension];
  }

  bool PointIsOverLowPlane(

      std::vector<double> const &point,
      unsigned int dimension) const
  {
    return (this->boxOrigin_[dimension] -
            this->boxSize_[dimension] / 2) >= point[dimension];
  }

  bool PointIsOnHighPlane(
      std::vector<double> const &point,
      unsigned int dimension) const
  {
    double compareTolerance = 1e-13;
    return std::abs(this->boxOrigin_[dimension] +
                    this->boxSize_[dimension] / 2 - point[dimension]) <
           compareTolerance;
  }

  bool PointIsOnLowPlane(
      std::vector<double> const &point,
      unsigned int dimension) const
  {
    double compareTolerance = 1e-13;
    return std::abs(this->boxOrigin_[dimension] -
                    this->boxSize_[dimension] / 2 - point[dimension]) <
           compareTolerance;
  }

  bool VerticesMatch(
      std::vector<double> const &vertex1,
      std::vector<double> const &vertex2) const
  {
    double compareTolerance = 1e-7;
    if (std::abs(vertex2[0] - vertex1[0]) > compareTolerance)
      return false;
    if (std::abs(vertex2[1] - vertex1[1]) > compareTolerance)
      return false;
    if (std::abs(vertex2[2] - vertex1[2]) > compareTolerance)
      return false;

    return true;
  }

  std::vector<unsigned int> ShiftVertices(

      std::vector<std::vector<double>> &vertices) const
  {
    std::vector<unsigned int> shifted_lines;
    shifted_lines.reserve(static_cast<int>(vertices.size() * 0.2));
    for (unsigned int i = 0; i < vertices.size(); ++i)
    {
      for (auto dim : {0, 1, 2})
      {
        if (!true)
          continue;

        if (this->PointIsOverHighPlane(vertices[i], dim))
        {
          this->ShiftPointDown(vertices[i], this->boxSize_, dim);
          for (unsigned int k = 0; k < node_to_edges_[i].size(); ++k)
            shifted_lines.push_back(node_to_edges_[i][k]);
        }

        else if (this->PointIsOverLowPlane(vertices[i], dim))
        {
          this->ShiftPointUp(vertices[i], this->boxSize_, dim);

          for (unsigned int k = 0; k < node_to_edges_[i].size(); ++k)
            shifted_lines.push_back(node_to_edges_[i][k]);
        }
      }
    }
    return shifted_lines;
  }

  void ShiftPointDown(std::vector<double> &point,
                      std::vector<double> const &boxSize,
                      unsigned int dim) const
  {
    point[dim] -= boxSize[dim];
  }

  void ShiftPointUp(std::vector<double> &point,
                    std::vector<double> const &boxSize,
                    unsigned int dim) const
  {
    point[dim] += boxSize[dim];
  }

  bool UnShift1D(
      int dim, double &d, double const &ref, double const &X) const
  {
    bool unshifted = false;

    if (not true)
      return unshifted;

    double x = d + X;

    if (x - ref < -0.5 * boxSize_[dim])
    {
      unshifted = true;
      d += boxSize_[dim];
    }
    else if (x - ref > 0.5 * boxSize_[dim])
    {
      unshifted = true;
      d -= boxSize_[dim];
    }

    return unshifted;
  }

  void UnShift3D(
      std::vector<double> &d, std::vector<double> const &ref, std::vector<double> const X = std::vector<double>(3, 0.0)) const
  {
    for (int dim = 0; dim < 3; ++dim)
      UnShift1D(dim, d[dim], ref[dim], X[dim]);
  }

  void get_unshifted_dir_vec(
      std::vector<double> x_1, std::vector<double> const &x_2, std::vector<double> &dirvec) const
  {
    UnShift3D(x_1, x_2);

    for (int idim = 0; idim < 3; ++idim)
      dirvec[idim] = x_1[idim] - x_2[idim];
  };

  std::vector<int> Permutation(int number,
                               std::mt19937 &gen,
                               std::uniform_real_distribution<> &dis_uni) const
  {
    // auxiliary variable
    int j = 0;

    // result vector initialized with ordered numbers from 0 to N-1
    std::vector<int> randorder(number, 0);
    for (int i = 0; i < (int)randorder.size(); i++)
      randorder[i] = i;

    for (int i = 0; i < number; ++i)
    {
      // generate random number between 0 and i
      j = (int)std::floor((i + 1.0) * dis_uni(gen));

      /*exchange values at positions i and j (note: value at position i is i due to above
     *initialization and because so far only positions <=i have been changed*/
      randorder[i] = randorder[j];
      randorder[j] = i;
    }

    return randorder;
  }

  double l2_norm(std::vector<double> const &u) const
  {
    double accum = 0.;
    for (int i = 0; i < u.size(); ++i)
    {
      accum += u[i] * u[i];
    }
    return sqrt(accum);
  }
  double l2_norm_dist_two_points(
      std::vector<double> x_1, std::vector<double> const &x_2) const
  {
    UnShift3D(x_1, x_2);
    std::vector<double> dirvec(3, 0.0);

    for (int idim = 0; idim < 3; ++idim)
      dirvec[idim] = x_1[idim] - x_2[idim];

    return l2_norm(dirvec);
  }

  void ComputeCosineDistributionOfNode(

      const unsigned int i_node,
      std::vector<double> &dir_vec_1,
      std::vector<double> &dir_vec_2,
      double interval_size_cosines,
      std::vector<std::vector<double>> &node_cosine_to_bin,
      std::vector<double> &cosine_distribution)
  {
    // undo old
    for (unsigned int p = 0; p < node_cosine_to_bin[i_node].size(); ++p)
      --cosine_distribution[node_cosine_to_bin[i_node][p]];

    node_cosine_to_bin[i_node].clear();

    // loop over all edges of the respective node
    for (unsigned int i = 0; i < edge_map_[i_node].size(); ++i)
    {
      // edge 1
      unsigned int edge_1 = edge_map_[i_node][i];

      // get direction vector
      for (unsigned int idim = 0; idim < 3; ++idim)
      {
        if (this->edges_[edge_1][0] == i_node)
          get_unshifted_dir_vec(vertices_[edges_[edge_1][1]],
                                vertices_[edges_[edge_1][0]], dir_vec_1);
        else
          get_unshifted_dir_vec(vertices_[edges_[edge_1][0]],
                                vertices_[edges_[edge_1][1]], dir_vec_1);
      }

      for (unsigned int j = i + 1; j < edge_map_[i_node].size(); ++j)
      {
        // edge 2
        unsigned int edge_2 = edge_map_[i_node][j];

        // get direction vector
        for (unsigned int idim = 0; idim < 3; ++idim)
        {
          if (this->edges_[edge_2][0] == i_node)
            get_unshifted_dir_vec(vertices_[edges_[edge_2][1]],
                                  vertices_[edges_[edge_2][0]], dir_vec_2);
          else
            get_unshifted_dir_vec(vertices_[edges_[edge_2][0]],
                                  vertices_[edges_[edge_2][1]], dir_vec_2);
        }

        // compute cosine
        double curr_cosine = std::inner_product(std::begin(dir_vec_1), std::end(dir_vec_1), std::begin(dir_vec_2), 0.0);
        curr_cosine /= l2_norm(dir_vec_1) * l2_norm(dir_vec_2);

        // add cosine
        unsigned int bin = std::floor((curr_cosine + 1.0) / interval_size_cosines);
        bin = (bin >= cosine_distribution.size()) ? (cosine_distribution.size() - 1) : bin;
        ++cosine_distribution[bin];
        node_cosine_to_bin[i_node].push_back(bin);
      }
    }
  }

  /*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
  void UpdateLengthDistributionOfLine(

      const unsigned int i_edge,
      double length_norm_fac,
      std::vector<double> &dir_vec_1,
      double interval_size_lengths,
      std::vector<double> &edge_length_to_bin,
      std::vector<double> &length_distribution) const
  {
    if (edge_length_to_bin[i_edge] > -0.1)
      --length_distribution[edge_length_to_bin[i_edge]];

    unsigned int node_1 = this->edges_[i_edge][0];
    unsigned int node_2 = this->edges_[i_edge][1];

    double curr_new_length = l2_norm_dist_two_points(this->vertices_[node_1],
                                                     this->vertices_[node_2]) *
                             length_norm_fac;

    unsigned int curr_bin = std::floor(curr_new_length / interval_size_lengths);
    curr_bin = curr_bin >= (length_distribution.size()) ? (length_distribution.size() - 1) : curr_bin;
    edge_length_to_bin[i_edge] = curr_bin;
    ++length_distribution[curr_bin];
  }

  /*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
  void RevertUpdateLengthDistributionOfLine(
      std::set<unsigned int> const &edges_to_revert,
      std::vector<double> &length_distribution,
      std::vector<double> &edge_length_to_bin,
      std::vector<double> const &edge_length_to_bin_backup) const
  {
    for (auto const &i_edge : edges_to_revert)
    {
      --length_distribution[edge_length_to_bin[i_edge]];
      ++length_distribution[edge_length_to_bin_backup[i_edge]];

      edge_length_to_bin[i_edge] = edge_length_to_bin_backup[i_edge];
    }
  }

  void UpdateBackupOfCosineDistribution(
      std::set<unsigned int> const &nodes_to_revert,
      std::vector<std::vector<double>> const &node_cosine_to_bin,
      std::vector<std::vector<double>> &node_cosine_to_bin_backup) const
  {
    for (auto const &i_node : nodes_to_revert)
    {
      node_cosine_to_bin_backup[i_node] = node_cosine_to_bin[i_node];
    }
  }

  /*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
  void RevertUpdateOfNodes(
      std::set<unsigned int> const &nodes_to_revert,
      std::vector<std::vector<double>> const &nodes_backup)
  {
    for (auto const &i_node : nodes_to_revert)
    {
      this->vertices_[i_node] = nodes_backup[i_node];
    }
  }

  /*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
  void UpdateBackUpOfNodes(
      std::set<unsigned int> const &nodes_to_revert,
      std::vector<std::vector<double>> &nodes_backup)
  {
    for (auto const &i_node : nodes_to_revert)
    {
      nodes_backup[i_node] = this->vertices_[i_node];
    }
  }

  /*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
  void RevertUpdateOfEdges(
      std::set<unsigned int> const &edges_to_revert,
      std::vector<std::vector<unsigned int>> const &uniqueVertexEdgePartners_backup)
  {
    for (auto const &i_edge : edges_to_revert)
    {
      edges_[i_edge] = uniqueVertexEdgePartners_backup[i_edge];
    }
  }

  /*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
  void UpdateBackUpOfEdges(
      std::set<unsigned int> const &edges_to_revert,
      std::vector<std::vector<unsigned int>> &uniqueVertexEdgePartners_backup)
  {
    for (auto const &i_edge : edges_to_revert)
    {
      uniqueVertexEdgePartners_backup[i_edge] = edges_[i_edge];
    }
  }

  void UpdateBackupOfLineDistribution(
      std::set<unsigned int> const &edges_to_revert,
      std::vector<double> const &edge_length_to_bin,
      std::vector<double> &edge_length_to_bin_backup) const
  {
    for (auto const &i_edge : edges_to_revert)
    {
      edge_length_to_bin_backup[i_edge] = edge_length_to_bin[i_edge];
    }
  }

  void RevertComputeCosineDistributionOfNode(
      std::set<unsigned int> const &nodes_to_revert,
      std::vector<double> &cosine_distribution,
      std::vector<std::vector<double>> &node_cosine_to_bin,
      std::vector<std::vector<double>> const &node_cosine_to_bin_backup) const
  {
    for (auto const &i_node : nodes_to_revert)
    {
      for (unsigned int i = 0; i < node_cosine_to_bin[i_node].size(); ++i)
      {
        --cosine_distribution[node_cosine_to_bin[i_node][i]];
      }

      for (unsigned int i = 0; i < node_cosine_to_bin_backup[i_node].size(); ++i)
      {
        ++cosine_distribution[node_cosine_to_bin_backup[i_node][i]];
      }

      node_cosine_to_bin[i_node] = node_cosine_to_bin_backup[i_node];
    }
  }

  /*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
  void EnergyLineLindstrom(
      double &curr_energy_line,
      unsigned int num_lines,
      double interval_size_lengths,
      std::vector<double> const &length_distribution) const
  {
    curr_energy_line = 0.0;

    double static mue = -0.3;   // mean
    double static sigma = 0.68; // standard deviation

    double one_sixth = 1.0 / 6.0;
    double curr_x = 0.0;
    double F = 0.0;
    double M = 0.0;
    double S = 0.0;
    for (unsigned int p = 0; p < length_distribution.size(); ++p)
    {
      curr_x = interval_size_lengths * p + interval_size_lengths * 0.5;
      F = (0.5 + 0.5 * std::erf((std::log(curr_x) - mue) / (sigma * M_SQRT2)));
      if (p > 0)
        M += length_distribution[p - 1];
      S = M - num_lines * F - 0.5;
      curr_energy_line += length_distribution[p] * (one_sixth * (length_distribution[p] + 1.0) *
                                                        (6.0 * S + 2.0 * length_distribution[p] + 1.0) +
                                                    S * S);
    }
    curr_energy_line *= 1.0 / (num_lines * num_lines);
  }

  void EnergyCosineLindstrom(
      double &curr_energy_cosine,
      double interval_size_cosines,
      std::vector<double> const &cosine_distribution,
      unsigned int num_cosines) const
  {
    curr_energy_cosine = 0.0;

    double static b_1 = 0.646666666666667 / 2.0;
    double static b_2 = -0.126666666666667 / 4.0;
    double static b_3 = 0.0200000000000001 / 6.0;

    double one_sixth = 1.0 / 6.0;
    double power_2 = 0.0;
    double curr_x = 0.0;
    double F = 0.0;
    double M = 0.0;
    double S = 0.0;
    for (unsigned int p = 0; p < cosine_distribution.size(); ++p)
    {
      curr_x = interval_size_cosines * p + interval_size_cosines * 0.5 - 1.0;
      power_2 = (1.0 - curr_x) * (1.0 - curr_x);
      F = -1.0 * b_1 * power_2 - b_2 * power_2 * power_2 - b_3 * power_2 * power_2 * power_2 + 1.0;
      if (p > 0)
        M += cosine_distribution[p - 1];
      S = M - num_cosines * F - 0.5;
      curr_energy_cosine += cosine_distribution[p] * (one_sixth * (cosine_distribution[p] + 1.0) *
                                                          (6.0 * S + 2.0 * cosine_distribution[p] + 1.0) +
                                                      S * S);
    }

    curr_energy_cosine *= 1.0 / (num_cosines * num_cosines);
  }

  void RemoveDoubles()
  {
    auto start_remove = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<double>> vtxs_shifted = this->vertices_;

    std::vector<unsigned int> shifted_lines = this->ShiftVertices(vtxs_shifted);

    // make unique
    std::set<int> s;
    unsigned size = shifted_lines.size();
    for (unsigned int k = 0; k < size; ++k)
      s.insert(shifted_lines[k]);
    shifted_lines.assign(s.begin(), s.end());

    unsigned int count = 1;
    for (unsigned int i = 0; i < shifted_lines.size(); ++i)
    {
      if ((edges_[shifted_lines[i]][0] == INT32_MAX) || (edges_[shifted_lines[i]][1] == INT32_MAX))
      {
        continue;
      }

      for (unsigned int j = 0; j < edges_.size(); ++j)
      {
        if (shifted_lines[i] == j or
            (edges_[j][0] == INT32_MAX) or (edges_[j][1] == INT32_MAX))
        {
          continue;
        }
        if ((edges_[shifted_lines[i]][0] == INT32_MAX) || (edges_[shifted_lines[i]][1] == INT32_MAX))
        {
          continue;
        }

        std::vector<double> edge1vtx1 = vtxs_shifted[edges_[shifted_lines[i]][0]];
        std::vector<double> edge1vtx2 = vtxs_shifted[edges_[shifted_lines[i]][1]];
        std::vector<double> edge2vtx1 = vtxs_shifted[edges_[j][0]];
        std::vector<double> edge2vtx2 = vtxs_shifted[edges_[j][1]];

        if ((this->VerticesMatch(edge1vtx1, edge2vtx1) && this->VerticesMatch(edge1vtx2, edge2vtx2)) ||
            (this->VerticesMatch(edge1vtx1, edge2vtx2) && this->VerticesMatch(edge1vtx2, edge2vtx1)))
        {
          edges_[shifted_lines[i] > j ? shifted_lines[i] : j] = std::vector<unsigned int>{INT32_MAX, INT32_MAX};
        }
      }

      if (i >= 0.05 * count * shifted_lines.size() - 1)
      {
        auto stop_remove = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_remove = stop_remove - start_remove;
        std::cout << "   " << elapsed_remove.count() / 60 << " minutes needed for " << std::flush;
        std::cout << "   " << static_cast<int>(0.05 * count * 100.0) << " %\r" << std::flush;
        ++count;
      }
    }

    // size_t edge1Id = 0;
    // unsigned int count = 1;
    // for (auto &edge1 : this->uniqueVertexEdgePartners_)
    // {
    //   if ((edge1[0] == INT32_MAX) || (edge1[1] == INT32_MAX))
    //   {
    //     ++edge1Id;
    //     continue;
    //   }

    //   size_t edge2Id = 0;
    //   for (auto &edge2 : this->uniqueVertexEdgePartners_)
    //   {
    //     if (edge1Id >= edge2Id)
    //     {
    //       ++edge2Id;
    //       continue;
    //     }

    //     if ((edge2[0] == INT32_MAX) || (edge2[1] == INT32_MAX))
    //     {
    //       ++edge2Id;
    //       continue;
    //     }

    //     std::vector<double> edge1vtx1 = vtxs_shifted[edge1[0]];
    //     std::vector<double> edge1vtx2 = vtxs_shifted[edge1[1]];
    //     std::vector<double> edge2vtx1 = vtxs_shifted[edge2[0]];
    //     std::vector<double> edge2vtx2 = vtxs_shifted[edge2[1]];

    //     if ((this->VerticesMatch(edge1vtx1, edge2vtx1) && this->VerticesMatch(edge1vtx2, edge2vtx2)) ||
    //         (this->VerticesMatch(edge1vtx1, edge2vtx2) && this->VerticesMatch(edge1vtx2, edge2vtx1)))
    //     {
    //       edge2 = std::vector<unsigned int>{INT32_MAX, INT32_MAX};
    //     }
    //     ++edge2Id;
    //   }

    //   if (edge1Id >= 0.05 * count * uniqueVertexEdgePartners_.size() - 1)
    //   {
    //     auto stop_remove = std::chrono::high_resolution_clock::now();
    //     std::chrono::duration<double> elapsed_remove = stop_remove - start_remove;
    //     std::cout << "   " << elapsed_remove.count() / 60 << " minutes needed for " << std::flush;
    //     std::cout << "   " << static_cast<int>(0.05 * count * 100.0) << " %\r" << std::flush;
    //     ++count;
    //   }

    //   ++edge1Id;
    // }
    std::cout << std::endl;

    this->edges_.erase(
        std::remove(this->edges_.begin(),
                    this->edges_.end(),
                    std::vector<unsigned int>{INT32_MAX, INT32_MAX}),
        this->edges_.end());
  }

  void AdaptConnectivity(std::mt19937 &gen, std::uniform_real_distribution<> &dis_uni)
  {

    std::cout << "\n4) Adapting valency distribution " << std::endl;
    auto start_valency = std::chrono::high_resolution_clock::now();

    unsigned int num_nodes = this->vertices_map_.size();
    unsigned int num_lines = this->edges_.size();

    // put all vertex ids in vector to ease random draw of vertex
    vertices_for_random_draw_.reserve(vertices_map_.size());
    for (auto const &iter : vertices_map_)
      vertices_for_random_draw_.push_back(iter.first);

    // Adapt valency distribution to collagen network according to Nan2018
    // "Realizations of highly heterogeneous collagen networks via stochastic reconstruction
    //  for micromechanical analysis of tumor cell invasion" figure 6
    std::uniform_int_distribution<> dis_rand_line(0, 3);
    std::vector<double> dir_1(3, 0.0);
    std::vector<double> dir_2(3, 0.0);

    unsigned int num_z_3 = std::floor(0.72 * num_nodes);
    unsigned int num_z_4 = std::floor(0.2 * num_nodes);
    unsigned int num_z_5 = std::floor(0.054 * num_nodes);
    unsigned int num_z_6 = std::floor(0.011 * num_nodes);

    num_z_4 += num_nodes - (num_z_3 + num_z_4 + num_z_5 + num_z_6);

    // first take care of z = 6 (starting from all being z = 4)
    std::uniform_int_distribution<> rand_node(0, vertices_for_random_draw_.size() - 1);

    unsigned int num_curr_z_6 = 0;
    unsigned int num_curr_z_5 = 0;

    while (num_curr_z_6 < num_z_6)
    {
      unsigned int node_1 = vertices_for_random_draw_[rand_node(gen)];

      if (node_to_edges_[node_1].size() == 4)
      {
        for (int i = 0; i < 2; ++i)
        {
          bool success = false;

          unsigned int node_2 = 0;
          while (not success)
          {
            node_2 = vertices_for_random_draw_[rand_node(gen)];
            if (node_2 == node_1 or node_to_edges_[node_2].size() != 4)
              continue;

            // check distance
            for (unsigned int dim = 0; dim < 3; ++dim)
              dir_1[dim] = vertices_map_[node_1][dim] - vertices_map_[node_2][dim];

            UnShift3D(dir_1, dir_2);

            if (l2_norm(dir_1) < one_third * this->boxSize_[0])
            {
              success = true;
              break;
            }
          }

          // add line to these nodes
          std::vector<unsigned int> new_line(2, 0);
          new_line[0] = node_1;
          new_line[1] = node_2;
          node_to_edges_[node_1].push_back(edges_.size());
          node_to_edges_[node_2].push_back(edges_.size());
          edges_.push_back(new_line);
          ++num_curr_z_5;
        }
        ++num_curr_z_6;
      }
    }

    // now take care of z = 5
    while (num_curr_z_5 < num_z_5)
    {
      unsigned int node_1 = vertices_for_random_draw_[rand_node(gen)];

      bool success = false;
      if (node_to_edges_[node_1].size() == 4)
      {
        unsigned int node_2 = 0;
        while (not success)
        {
          node_2 = vertices_for_random_draw_[rand_node(gen)];
          if (node_2 == node_1 or node_to_edges_[node_2].size() != 4)
            continue;

          // check distance
          for (unsigned int dim = 0; dim < 3; ++dim)
            dir_1[dim] = vertices_map_[node_1][dim] - vertices_map_[node_2][dim];

          UnShift3D(dir_1, dir_2);

          if (l2_norm(dir_1) < one_third * this->boxSize_[0])
          {
            success = true;
            break;
          }
        }

        // add line to these nodes
        std::vector<unsigned int> new_line(2, 0);
        new_line[0] = node_1;
        new_line[1] = node_2;
        node_to_edges_[node_1].push_back(edges_.size());
        node_to_edges_[node_2].push_back(edges_.size());
        edges_.push_back(new_line);
        num_curr_z_5 += 2;
      }
    }

    // now do with z = 3
    unsigned int num_found = 0;
    unsigned int max_try = 100;
    std::vector<int> random_order = Permutation(vertices_for_random_draw_.size(), gen, dis_uni);
    // do twice for better results
    for (int s = 0; s < 2; ++s)
    {
      for (unsigned int rand_node_i = 0; rand_node_i < random_order.size(); ++rand_node_i)
      {
        int i_node = vertices_for_random_draw_[random_order[rand_node_i]];

        if (node_to_edges_[i_node].size() == 4)
        {
          unsigned int try_i = 1;
          bool success = false;
          unsigned int random_line = dis_rand_line(gen);
          unsigned int second_affected_node = -1;
          do
          {
            if (edges_[node_to_edges_[i_node][random_line]][0] == i_node)
              second_affected_node = edges_[node_to_edges_[i_node][random_line]][1];
            else
              second_affected_node = edges_[node_to_edges_[i_node][random_line]][0];

            if (node_to_edges_[second_affected_node].size() == 4)
            {
              success = true;
              break;
            }

            ++try_i;

          } while (try_i < max_try);

          if (not success)
            continue;

          //! don't erase yet! this destroys the order in edgeIds_valid!
          this->edges_[node_to_edges_[i_node][random_line]] =
              std::vector<unsigned int>{INT32_MAX, INT32_MAX};

          int index = std::distance(node_to_edges_[second_affected_node].begin(), std::find(node_to_edges_[second_affected_node].begin(),
                                                                                            node_to_edges_[second_affected_node].end(),
                                                                                            node_to_edges_[i_node][random_line]));

          node_to_edges_[second_affected_node].erase(node_to_edges_[second_affected_node].begin() + index);
          node_to_edges_[i_node].erase(node_to_edges_[i_node].begin() + random_line);

          num_found += 2;
        }

        if (num_found >= num_z_3)
          break;
      }
    }
    // now really remove edges
    // erase now
    this->edges_.erase(std::remove(this->edges_.begin(),
                                   this->edges_.end(),
                                   std::vector<unsigned int>{INT32_MAX, INT32_MAX}),
                       this->edges_.end());
    num_lines = edges_.size();

    // recompute node to edges after adaption of connectivity
    node_to_edges_.clear();
    node_to_edges_ = std::vector<std::vector<unsigned int>>(vertices_.size(), std::vector<unsigned int>());
    for (unsigned int i_edge = 0; i_edge < edges_.size(); ++i_edge)
    {
      node_to_edges_[this->edges_[i_edge][0]].push_back(i_edge);
      node_to_edges_[this->edges_[i_edge][1]].push_back(i_edge);
    }

    auto stop_valency = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_valency = stop_valency - start_valency;
    std::cout << "   Adaption of valency distribution is done now. It took " << elapsed_valency.count() / 60 << " minutes. \n";
  }

  void SimulatedAnnealing(uint mode, std::mt19937 &gen, std::uniform_real_distribution<> &dis_uni)
  {

    // put all vertex ids in vector to ease random draw of vertex
    vertices_for_random_draw_.reserve(vertices_map_.size());
    for (auto const &iter : vertices_map_)
      vertices_for_random_draw_.push_back(iter.first);

    //*************************************************************************************
    // DO SIMULATED ANNEALING
    bool do_simulated_annealing = true;
    if (not do_simulated_annealing)
      return;
    //*************************************************************************************

    std::cout << "\n5) Starting Simulated Annealing\n\n\n";
    // time measurement start
    auto start = std::chrono::high_resolution_clock::now();

    unsigned int num_nodes = vertices_map_.size();
    unsigned int num_lines = edges_.size();

    // ... some more
    // to select interchange movement
    std::uniform_int_distribution<> dis_action(1, 2);
    // to select random node
    std::uniform_int_distribution<> dis_node(0, num_nodes - 1);
    // to select random line
    std::uniform_int_distribution<> dis_line(0, num_lines - 1);
    // to select random node movement
    std::uniform_real_distribution<> dis_node_move(-1, 1);

    // ... and even more
    std::vector<double> rand_new_node_pos(3, 0.0);
    std::vector<std::vector<double>> uniqueVertices_backup(this->vertices_);
    std::vector<std::vector<unsigned int>> uniqueVertexEdgePartners_backup(this->edges_);
    unsigned int random_line_1 = 0;
    unsigned int random_line_2 = 0;
    unsigned int iter = 0;
    double curr_energy_line = 0.0;
    double last_energy_line = 0.0;
    double curr_energy_cosine = 0.0;
    double last_energy_cosine = 0.0;
    double temperature = temperature_inital;
    double delta_energy = 0.0;
    // normalize lengths according to Lindström
    double length_norm_fac = 1.0 / std::pow((num_nodes / (this->boxSize_[0] * this->boxSize_[1] * this->boxSize_[2])), -1.0 / 3.0);

    // for binning
    std::vector<unsigned int> m_j_lengths;
    std::vector<unsigned int> m_j_cosines;
    double interval_size_lengths = 5.0 / p_num_bins_lengths;
    double interval_size_cosines = 2.0 / p_num_bins_cosines;

    // build edge_map_
    for (unsigned int i_edge = 0; i_edge < edges_.size(); ++i_edge)
    {
      edge_map_[this->edges_[i_edge][0]].push_back(i_edge);
      edge_map_[this->edges_[i_edge][1]].push_back(i_edge);
    }

    // compute cosine distribution
    std::vector<double> cosine_distribution(p_num_bins_cosines, 0.0);
    std::vector<std::vector<double>> node_cosine_to_bin(vertices_.size(), std::vector<double>());
    std::vector<double> dir_vec_1(3, 0.0);
    std::vector<double> dir_vec_2(3, 0.0);
    for (auto const &i_node : edge_map_)
    {
      ComputeCosineDistributionOfNode(i_node.first, dir_vec_1, dir_vec_2,
                                      interval_size_cosines, node_cosine_to_bin, cosine_distribution);
    }
    std::vector<std::vector<double>> node_cosine_to_bin_backup(node_cosine_to_bin);

    unsigned int num_cosines = 0;
    for (unsigned int i_c = 0; i_c < cosine_distribution.size(); ++i_c)
      num_cosines += cosine_distribution[i_c];

    // compute length distribution
    std::vector<double> length_distribution(p_num_bins_lengths, 0.0);
    std::vector<double> edge_length_to_bin(edges_.size(), -1.0);
    for (unsigned int i_edge = 0; i_edge < num_lines; ++i_edge)
    {
      UpdateLengthDistributionOfLine(i_edge, length_norm_fac, dir_vec_1,
                                     interval_size_lengths, edge_length_to_bin, length_distribution);
    }
    std::vector<double> edge_length_to_bin_backup(edge_length_to_bin);

    // compute for first iteration
    EnergyLineLindstrom(last_energy_line, num_lines, interval_size_lengths, length_distribution);
    EnergyCosineLindstrom(last_energy_cosine, interval_size_cosines, cosine_distribution, num_cosines);

    // write initial distributions
    // output initial filament lengths
    std::ofstream filLen_file_initial(this->outputPrefix_.string() + "_fil_lengths_initial.txt");
    filLen_file_initial << "fil_lengths\n";
    for (unsigned int filId = 0; filId < this->edges_.size(); ++filId)
      filLen_file_initial << this->GetFilamentLength(filId) * length_norm_fac << "\n";

    // print final cosine distribution
    std::ofstream filcoshisto_initial_file(this->outputPrefix_.string() + "_cosine_histo_initial.txt");
    filcoshisto_initial_file << "cosine\n";
    for (unsigned int i_c = 0; i_c < cosine_distribution.size(); ++i_c)
      for (unsigned int j_c = 0; j_c < cosine_distribution[i_c]; ++j_c)
        filcoshisto_initial_file << interval_size_cosines * i_c + interval_size_cosines * 0.5 - 1.0 << "\n";

    // write temperature and energies to file
    std::ofstream fil_obj_function(this->outputPrefix_.string() + "_obj_function.txt");
    fil_obj_function << "step, temperature, length, cosine, total \n";

    //---------------------------
    // START SIMULATED ANNEALING
    //---------------------------
    do
    {
      unsigned int action = mode;

      // ************************************************
      // type one: move random point in random direction
      // ************************************************
      if (action == 1 || action == 3)
      {
        bool success = true;
        unsigned int subiter = 0;
        do
        {
          ++subiter;
          success = true;
          // select a random node
          unsigned int rand_node_id = vertices_for_random_draw_[dis_node(gen)];

          // update position of this vertex
          for (unsigned int idim = 0; idim < 3; ++idim)
          {
            this->vertices_[rand_node_id][idim] += dis_node_move(gen) * max_movement;
          }

          // recompute length and cosine distribution of affected nodes
          std::set<unsigned int> affected_nodes;
          std::set<unsigned int> affected_lines;
          for (auto const &iter_edges : node_to_edges_[rand_node_id])
          {
            affected_nodes.insert(this->edges_[iter_edges][0]);
            affected_nodes.insert(this->edges_[iter_edges][1]);
            affected_lines.insert(iter_edges);

            // check if line is longer than 1/3 of boxlength (we do not want this to
            // ensure that our RVE stays representative)
            if (GetEdgeLength(iter_edges) / length_norm_fac > one_third * this->boxSize_[0])
            {
              success = false;
              break;
            }
          }

          if (success == false)
          {
            RevertUpdateOfNodes(affected_nodes, uniqueVertices_backup);
            continue;
          }

          for (auto const &iter_nodes : affected_nodes)
            ComputeCosineDistributionOfNode(iter_nodes, dir_vec_1, dir_vec_2,
                                            interval_size_cosines, node_cosine_to_bin, cosine_distribution);

          for (auto const &iter_edges : affected_lines)
            UpdateLengthDistributionOfLine(iter_edges, length_norm_fac, dir_vec_1,
                                           interval_size_lengths, edge_length_to_bin, length_distribution);

          // compute energies
          // 1.) line
          EnergyLineLindstrom(curr_energy_line, num_lines, interval_size_lengths, length_distribution);

          // 2.) cosine
          EnergyCosineLindstrom(curr_energy_cosine, interval_size_cosines, cosine_distribution, num_cosines);

          // compute delta E
          delta_energy = weight_line * (curr_energy_line - last_energy_line) +
                         weight_cosine * (curr_energy_cosine - last_energy_cosine);

          if ((delta_energy < 0.0) or (dis_uni(gen) < std::exp(-delta_energy / temperature)))
          {
            last_energy_line = curr_energy_line;
            last_energy_cosine = curr_energy_cosine;
            UpdateBackUpOfNodes(affected_nodes, uniqueVertices_backup);
            UpdateBackupOfCosineDistribution(affected_nodes, node_cosine_to_bin, node_cosine_to_bin_backup);
            UpdateBackupOfLineDistribution(affected_lines, edge_length_to_bin, edge_length_to_bin_backup);
            success = true;
          }
          else
          {
            RevertUpdateOfNodes(affected_nodes, uniqueVertices_backup);
            RevertComputeCosineDistributionOfNode(affected_nodes, cosine_distribution, node_cosine_to_bin, node_cosine_to_bin_backup);
            RevertUpdateLengthDistributionOfLine(affected_lines, length_distribution, edge_length_to_bin, edge_length_to_bin_backup);
            success = false;
          }
        } while ((success == false) and (subiter < max_subiter));

        if (iter % screen_output_every == 0)
        {
          std::cout << "line energy move 1 " << curr_energy_line << std::endl;
          std::cout << "cosine energy move 1 " << curr_energy_cosine << std::endl;
          std::cout << " iter " << iter << std::endl;

          fil_obj_function << iter;
          fil_obj_function << ", " << temperature;
          fil_obj_function << ", " << curr_energy_line;
          fil_obj_function << ", " << curr_energy_cosine;
          fil_obj_function << ", " << curr_energy_line + curr_energy_cosine << "\n";
        }
      }

      // ************************************************
      // type two: change connection of two lines
      // ************************************************
      if (action == 2 || action == 3)
      {
        //      Current Conf.           Case_1             Case_2
        //      _____________       _____________      ______________
        //
        //      1 o------o 2         1 o     o 2        1 o       o 2
        //                              \   /             |       |
        //                               \ /              |       |
        //                                /               |       |
        //                               /  \             |       |
        //      3 o------o 4         3  o    o 4        3 o       o 4
        //
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        bool success = false;
        unsigned int subiter = 0;
        do
        {
          ++subiter;
          success = true;

          // select two (different) random lines
          random_line_1 = dis_line(gen);
          bool not_yet_found = true;
          std::vector<unsigned int> nodes_line_1(2, 0);
          std::vector<unsigned int> nodes_line_2(2, 0);

          nodes_line_1[0] = edges_[random_line_1][0];
          nodes_line_1[1] = edges_[random_line_1][1];

          unsigned int control_iter = 0;
          while (not_yet_found)
          {
            ++control_iter;
            // prevent code from getting stuck in this loop
            if (control_iter > 1e06)
            {
              std::cout << " Not able to find close enough second edge for option 2 " << std::endl;
              exit(0);
            }

            random_line_2 = dis_line(gen);
            nodes_line_2[0] = edges_[random_line_2][0];
            nodes_line_2[1] = edges_[random_line_2][1];

            // already check distance here
            bool to_far_away = false;
            for (unsigned int j = 0; j < 2; ++j)
            {
              if (l2_norm_dist_two_points(vertices_[nodes_line_1[j]], vertices_[nodes_line_2[0]]) > one_third * this->boxSize_[0] or
                  l2_norm_dist_two_points(vertices_[nodes_line_1[j]], vertices_[nodes_line_2[1]]) > one_third * this->boxSize_[0])
              {
                to_far_away = true;
                break;
              }
            }

            if (to_far_away)
              continue;

            // check if lines share a node, if so, choose another one
            if (not((nodes_line_1[0] == nodes_line_2[0]) or (nodes_line_1[0] == nodes_line_2[1]) or
                    (nodes_line_1[1] == nodes_line_2[0]) or (nodes_line_1[1] == nodes_line_2[1])))
              not_yet_found = false;
          }

          std::set<unsigned int> affected_nodes;
          affected_nodes.insert(edges_[random_line_1][0]);
          affected_nodes.insert(edges_[random_line_1][1]);
          affected_nodes.insert(edges_[random_line_2][0]);
          affected_nodes.insert(edges_[random_line_2][1]);

          std::set<unsigned int> affected_lines;
          affected_lines.insert(random_line_1);
          affected_lines.insert(random_line_2);

          // get all nodes that are connected to the respective four nodes
          std::vector<std::set<unsigned int>> nodes_to_nodes_1(2, std::set<unsigned int>());
          std::vector<std::set<unsigned int>> nodes_to_nodes_2(2, std::set<unsigned int>());

          for (unsigned int j = 0; j < 2; ++j)
          {
            for (unsigned int k = 0; k < node_to_edges_[nodes_line_1[j]].size(); ++k)
            {
              nodes_to_nodes_1[j].insert(edges_[k][0]);
              nodes_to_nodes_1[j].insert(edges_[k][1]);
            }
          }

          for (unsigned int j = 0; j < 2; ++j)
          {
            for (unsigned int k = 0; k < node_to_edges_[nodes_line_2[j]].size(); ++k)
            {
              nodes_to_nodes_2[j].insert(edges_[k][0]);
              nodes_to_nodes_2[j].insert(edges_[k][1]);
            }
          }

          // decide which case is attempted first
          //        int mov_case = dis_action(gen);

          // we always try movement one first
          bool move_one_sucess = true;

          // next, check if one of the two new lines already exists for case 1
          if ((nodes_to_nodes_1[0].count(nodes_line_2[1])) or (nodes_to_nodes_1[1].count(nodes_line_2[0])))
            move_one_sucess = false;

          // update new connectivity
          if (move_one_sucess == true)
          {
            this->edges_[random_line_1][1] = nodes_line_2[1];
            this->edges_[random_line_2][1] = nodes_line_1[1];
          }

          if (move_one_sucess == false)
          {
            // check if one of the two new lines already exists for case 2
            if (nodes_to_nodes_1[0].count(nodes_line_2[0]) or nodes_to_nodes_1[1].count(nodes_line_2[1]))
            {
              success = false;
              continue;
            }

            // update new connectivity
            this->edges_[random_line_1][1] = nodes_line_2[0];
            this->edges_[random_line_2][0] = nodes_line_1[1];
          }

          // update length distribution
          UpdateLengthDistributionOfLine(random_line_1, length_norm_fac, dir_vec_1,
                                         interval_size_lengths, edge_length_to_bin, length_distribution);
          UpdateLengthDistributionOfLine(random_line_2, length_norm_fac, dir_vec_1,
                                         interval_size_lengths, edge_length_to_bin, length_distribution);
          // recompute cosine distribution of affected nodes
          for (unsigned int j = 0; j < 2; ++j)
          {
            ComputeCosineDistributionOfNode(nodes_line_1[j], dir_vec_1, dir_vec_2,
                                            interval_size_cosines, node_cosine_to_bin, cosine_distribution);
            ComputeCosineDistributionOfNode(nodes_line_2[j], dir_vec_1, dir_vec_2,
                                            interval_size_cosines, node_cosine_to_bin, cosine_distribution);
          }

          // compute energies
          // 1.) line
          EnergyLineLindstrom(curr_energy_line, num_lines, interval_size_lengths, length_distribution);

          // 2.) cosine
          EnergyCosineLindstrom(last_energy_cosine, interval_size_cosines, cosine_distribution, num_cosines);

          // compute delta E
          delta_energy = weight_line * (curr_energy_line - last_energy_line) +
                         weight_cosine * (curr_energy_cosine - last_energy_cosine);

          if ((delta_energy < 0.0) or (dis_uni(gen) < std::exp(-delta_energy / temperature)))
          {
            last_energy_line = curr_energy_line;
            last_energy_cosine = curr_energy_cosine;
            UpdateBackUpOfEdges(affected_lines, uniqueVertexEdgePartners_backup);
            UpdateBackupOfCosineDistribution(affected_nodes, node_cosine_to_bin, node_cosine_to_bin_backup);
            UpdateBackupOfLineDistribution(affected_lines, edge_length_to_bin, edge_length_to_bin_backup);
            success = true;
          }
          else
          {
            RevertUpdateOfEdges(affected_lines, uniqueVertexEdgePartners_backup);
            RevertComputeCosineDistributionOfNode(affected_nodes, cosine_distribution, node_cosine_to_bin, node_cosine_to_bin_backup);
            RevertUpdateLengthDistributionOfLine(affected_lines, length_distribution, edge_length_to_bin, edge_length_to_bin_backup);
            success = false;
          }
        } while (success == false and subiter < max_subiter);

        // screen output
        if (iter % screen_output_every == 0)
        {
          std::cout << "line energy move 2 " << curr_energy_line << std::endl;
          std::cout << "cosine energy move 2 " << curr_energy_cosine << std::endl;
          std::cout << " iter 2 " << iter << std::endl;
        }
      }

      // neither movement one nor two
      if (action != 1 && action != 2 && action != 3)
      {
        throw "You should not be here";
      }

      // according to Nan2018 (power law cooling schedule)
      if (iter % 1000 == 0)
      {
        temperature = std::pow(decay_rate_temperature, iter / 1000.0) * temperature_inital;
        std::cout << "temperature " << temperature << std::endl;
      }

      ++iter;
    } while ((iter < max_iter) and ((last_energy_line > tolerance) or (last_energy_cosine > tolerance)));

    // print final cosine distribution
    std::ofstream filcoshisto_file(this->outputPrefix_.string() + "_cosine_histo.txt");
    filcoshisto_file << "cosine\n";
    for (unsigned int i_c = 0; i_c < cosine_distribution.size(); ++i_c)
      for (unsigned int j_c = 0; j_c < cosine_distribution[i_c]; ++j_c)
        filcoshisto_file << interval_size_cosines * i_c + interval_size_cosines * 0.5 - 1.0 << "\n";

    std::ofstream filcos_file(this->outputPrefix_.string() + "_cosine_normal.txt");
    filcos_file << "bin, cosine \n";
    for (unsigned int i_c = 0; i_c < cosine_distribution.size(); ++i_c)
    {
      filcos_file << interval_size_cosines * i_c + interval_size_cosines * 0.5;
      filcos_file << ", " << cosine_distribution[i_c] - 1.0 << "\n";
    }

    // time measurement end
    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = stop - start;

    std::cout << "\nSimulation annealing took " << elapsed.count() / 60 << " minutes for " << iter << " iterations" << std::endl;
    std::cout << "Final line energy:   " << last_energy_line << std::endl;
    std::cout << "Final cosine energy: " << last_energy_cosine << std::endl;
  }

  void OutputGeometry()
  {
    std::ofstream partnersOutFile(this->outputPrefix_.string() + "_partners.out");
    for (auto partner : this->edges_)
    {
      partnersOutFile << partner[0] << " "
                      << partner[1] << std::endl;
    }
    partnersOutFile.close();
    std::ofstream verticesOutFile(this->outputPrefix_.string() + "_vertices.out");
    for (auto vertexI = 0; vertexI < this->vertices_.size(); vertexI++)
    {
      verticesOutFile << boost::lexical_cast<std::string>(this->vertices_[vertexI][0]) << " "
                      << boost::lexical_cast<std::string>(this->vertices_[vertexI][1]) << " "
                      << boost::lexical_cast<std::string>(this->vertices_[vertexI][2]) << " "
                      << this->vertexEdgeCount_[vertexI] << std::endl;
    }
    verticesOutFile.close();

    std::ofstream nodesToEdgesOutFile(this->outputPrefix_.string() + "_nodes_to_edges.out");
    for (auto edges : this->node_to_edges_)
    {
      for (auto edge : edges)
        nodesToEdgesOutFile << edge << " ";
      nodesToEdgesOutFile << std::endl;
    }
    nodesToEdgesOutFile.close();
  }

  void ReadGeometry()
  {
    this->edges_.clear();
    this->vertexEdgeCount_.clear();
    this->vertices_map_.clear();
    this->vertices_.clear();
    this->node_to_edges_.clear();
    std::string row;
    std::ifstream verticesFile = std::ifstream(this->inputPrefix_.string() + "_vertices.out");
    auto rowI = 0;
    while (!verticesFile.eof())
    {
      std::getline(verticesFile, row);
      if (row == "")
        continue;
      if (verticesFile.bad() || verticesFile.fail())
        break;
      std::stringstream s(row);
      double d;
      std::vector<double> vertex;
      auto i = 0;
      while (s >> d)
        if (++i <= 3)
          vertex.push_back(d);
        else if (i > 3)
          this->vertexEdgeCount_.push_back(d);
        else
          throw "Error in voronoi geometry loading.";
      this->vertices_.push_back(vertex);
      if (this->vertexEdgeCount_[rowI] != 0)
        this->vertices_map_.emplace(rowI, vertex);
      rowI++;
    }
    std::ifstream partnersFile = std::ifstream(this->inputPrefix_.string() + "_partners.out");
    while (!partnersFile.eof())
    {
      std::getline(partnersFile, row);
      if (row == "")
        continue;
      if (partnersFile.bad() || partnersFile.fail())
        break;
      std::stringstream s(row);
      uint u;
      std::vector<uint> partners;
      while (s >> u)
        partners.push_back(u);
      this->edges_.push_back(partners);
    }

    std::ifstream nodesToEdgesFile = std::ifstream(this->inputPrefix_.string() + "_nodes_to_edges.out");
    while (!nodesToEdgesFile.eof())
    {
      std::getline(nodesToEdgesFile, row);
      if (row == "") //! DO NOT continue; empty lines mean no edges for nodes (aka dead nodes)
        ;
      if (nodesToEdgesFile.bad() || nodesToEdgesFile.fail())
        break;
      std::stringstream s(row);
      uint u;
      std::vector<uint> edges;
      while (s >> u)
        edges.push_back(u);
      this->node_to_edges_.push_back(edges);
    }
  }
};
