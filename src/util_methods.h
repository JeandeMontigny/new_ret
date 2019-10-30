#ifndef UTILS_METHODS
#define UTILS_METHODS

#include "extended_objects.h"
#include "rgc_soma_bm.h"
#include "rgc_dendrite_bm.h"

namespace bdm {
  using namespace std;

  // define my cell creator
  static void CellCreator(double min, double max, int num_cells, int cell_type) {
    auto* sim = Simulation::GetActive();
    auto* rm = sim->GetResourceManager();
    auto* random = sim->GetRandom();

    for (int i = 0; i < num_cells; i++) {
      double x = random->Uniform(min + 100, max - 100);
      double y = random->Uniform(min + 100, max - 100);
      // RGCL thickness before cell death ~24
      double z = random->Uniform(min + 20, min + 34);
      z = 27;

      MyCell* cell = new MyCell({x, y, z});
      cell->SetDiameter(random->Uniform(7, 8));
      cell->SetCellType(cell_type);
      cell->SetPreviousPosition({x, y, z});
      // cell->AddBiologyModule(new Substance_secretion_BM());
      // cell->AddBiologyModule(new RGC_mosaic_BM());
      cell->AddBiologyModule(new Internal_clock_BM());
      cell->AddBiologyModule(new Dendrite_creation_BM());
      rm->push_back(cell);
    }
  }  // end CellCreator


  inline void WritePositions(int i, int seed) {
    auto* sim = Simulation::GetActive();
    auto* rm = sim->GetResourceManager();
    auto* param = sim->GetParam();
    ofstream position_file;
    stringstream position_file_name;
    position_file_name << Concat(param->output_dir_, "/results",
                                 seed, "/cells_position/")
                       << i << "_seed" << seed << ".txt";
    position_file.open(position_file_name.str());

    rm->ApplyOnAllElements([&](SimObject* so, SoHandle) {
      auto* cell = dynamic_cast<MyCell*>(so);
      if (cell) {
        auto position = cell->GetPosition();
        // type x y z
        position_file << cell->GetCellType() << " " << position[0] << " "
                      << position[1] << " " << position[2] << "\n";
      }
    });  // end for cell in simulation

    position_file.close();
  } // end WritePositions


  template <typename T>
  inline string SwcNeurites(T ne_so_ptr, int label_parent,
                            const Double3& soma_position) {
    auto* ne = static_cast<MyNeurite*>(&*ne_so_ptr);
    Double3 ne_position = ne->GetPosition();
    ne_position[0] = ne_position[0] - soma_position[0];
    ne_position[1] = ne_position[1] - soma_position[1];
    ne_position[2] = ne_position[2] - soma_position[2];
    double radius = ne->GetDiameter() / 2;
    string temps;

    ne->GetMySoma()->IncreaseLabel();
    // set explicitly the value of GetLabel() to avoid reference
    int current_label = ne->GetMySoma()->GetLabel();

    // if branching point
    if (ne->GetDaughterRight() != nullptr) {
      temps = Concat(temps, "\n", current_label, " 3 ", ne_position[0], " ",
                ne_position[1], " ", ne_position[2], " ", radius, " ", label_parent,
                SwcNeurites(ne->GetDaughterRight(), current_label, soma_position))
              .c_str();
      ne->GetMySoma()->IncreaseLabel();
    }
    // if straigh dendrite
    current_label = ne->GetMySoma()->GetLabel();
    if (ne->GetDaughterLeft() != nullptr) {
      temps = Concat(temps, "\n", current_label, " 3 ", ne_position[0], " ",
                ne_position[1], " ", ne_position[2], " ", radius, " ", label_parent,
                SwcNeurites(ne->GetDaughterLeft(), current_label, soma_position))
              .c_str();
    }
    // if ending point
    if (ne->GetDaughterLeft() == nullptr && ne->GetDaughterRight() == nullptr) {
      temps = Concat(temps, "\n", current_label, " 6 ", ne_position[0], " ",
                ne_position[1], " ", ne_position[2], " ", radius, " ", label_parent)
              .c_str();
    }

    return temps;
  } // end SwcNeurites


  inline void WriteSwc(int i, int seed) {
    auto* sim = Simulation::GetActive();
    auto* rm = sim->GetResourceManager();
    auto* param = sim->GetParam();

    rm->ApplyOnAllElements([&](SimObject* so, SoHandle) {
      auto* cell = dynamic_cast<MyCell*>(so);
      if (cell) {
        int cell_type = cell->GetCellType();
        auto cell_position = cell->GetPosition();
        ofstream swc_file;
        string swc_fileName = Concat(param->output_dir_,
          "/results", seed, "/swc_files/cell", cell->GetUid(),
          "_type", cell_type, "_seed", seed, "_step", i, ".swc").c_str();
        swc_file.open(swc_fileName);
        cell->SetLabel(1);
        swc_file << cell->GetLabel() << " 1 0 0 0 " << cell->GetDiameter() / 2
                 << " -1";

        for (auto& ne : cell->GetDaughters()) {
          swc_file << SwcNeurites(ne, 1, cell_position);
        }
        swc_file.close();
      }
    });  // end for cell in simulation

  } // end WriteSwc


  // RI computation
  inline double ComputeRi(vector<Double3> coord_list) {
    if (coord_list.size() < 2) {
      return 0;
    }
    vector<double> shortest_dist_list;
    for (unsigned int i = 0; i < coord_list.size();
         i++) {  // for each cell of same type in the simulation
      Double3 cell_position = coord_list[i];

      vector<double> shortestDist;
      // for each other cell of same type in the simulation
      for (unsigned int j = 0; j < coord_list.size(); j++) {
        Double3 otherCellPosition = coord_list[j];

        // get the distance between those 2 cells (x-y plan only)
        double tempsDistance =
            sqrt(pow(cell_position[0] - otherCellPosition[0], 2) +
                 pow(cell_position[1] - otherCellPosition[1], 2));
        // if cell is closer and is not itself
        if (tempsDistance != 0) {
          shortestDist.push_back(tempsDistance);
        }
      }
      // save the shortest distance to neighbour
      shortest_dist_list.push_back(
        *min_element(shortestDist.begin(), shortestDist.end()));
    }
    // compute mean
    double temps_sum = 0;
    for (unsigned int i = 0; i < shortest_dist_list.size(); i++) {
      temps_sum += shortest_dist_list[i];
    }
    double ave_shortest_dist = temps_sum / (double)shortest_dist_list.size();
    // compute std
    double temp = 0;
    for (unsigned int i = 0; i < shortest_dist_list.size(); i++) {
      double val = shortest_dist_list[i];
      double squrDiffToMean = pow(val - ave_shortest_dist, 2);
      temp += squrDiffToMean;
    }
    double mean_of_diffs = temp / (double)(shortest_dist_list.size());
    double std = sqrt(mean_of_diffs);

    return ave_shortest_dist / std;  // return RI
  }  // end ComputeRi


  inline double GetRi(int desired_type) {
    auto* sim = Simulation::GetActive();
    auto* rm = sim->GetResourceManager();
    vector<Double3> coord_list; // list of coordinate
    rm->ApplyOnAllElements([&](SimObject* so, SoHandle) {
      auto* cell = dynamic_cast<MyCell*>(so);
      if (cell) {
        auto cell_type = cell->GetCellType();
        if (cell_type == desired_type) {  // if cell is of the desired type
          auto position = cell->GetPosition();  // get its position
          coord_list.push_back(position);  // put cell coord in the list
        }
      }
    });  // end for cell in simulation
    return ComputeRi(coord_list);  // return RI for desired cell type
  }  // end GetRi


  inline vector<array<double, 2>> GetAllRI() {
    auto* sim = Simulation::GetActive();
    auto* rm = sim->GetResourceManager();
    vector<double> types_list;
    rm->ApplyOnAllElements([&](SimObject* so, SoHandle) {
      auto* cell = dynamic_cast<MyCell*>(so);
      if (cell) {
        auto cell_type = cell->GetCellType();
        // if type is not -1 and not in list
        if (find(types_list.begin(), types_list.end(), cell_type) == types_list.end()) {
          types_list.push_back(cell_type); // put cell coord in the list
        }
      }
    });  // end for cell in simulation
    vector<array<double, 2>> listRi;
    for (unsigned int i = 0; i < types_list.size(); i++) {
      listRi.push_back({GetRi(types_list[i]), types_list[i]});
    }
    return listRi;
  }


  inline double GetDeathRate(int num_cells) {
    auto* sim = Simulation::GetActive();
    auto* rm = sim->GetResourceManager();
    int cell_in_simu = 0;
    rm->ApplyOnAllElements([&](SimObject* so, SoHandle) {
      auto* cell = dynamic_cast<MyCell*>(so);
      if (cell) {
        cell_in_simu++;
      }
    });  // end for cell in simulation

    return (1 - ((double)cell_in_simu / num_cells)) * 100;
  } // end GetDeathRate

}

#endif
