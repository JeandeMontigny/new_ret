#ifndef UTILS_METHODS
#define UTILS_METHODS

#include "extended_objects.h"
#include "rgc_soma_bm.h"

namespace bdm {
  using namespace std;

  // define my cell creator
  static void CellCreator(double min, double max, int num_cells, int cellType) {
    auto* sim = Simulation::GetActive();
    auto* rm = sim->GetResourceManager();
    auto* random = sim->GetRandom();

    for (int i = 0; i < num_cells; i++) {
      double x = random->Uniform(min + 10, max - 10);
      double y = random->Uniform(min + 10, max - 10);
      // RGCL thickness before cell death ~24
      double z = random->Uniform(min + 20, min + 34);

      MyCell* cell = new MyCell({x, y, z});
      cell->SetDiameter(random->Uniform(7, 8));
      cell->SetCellType(cellType);
      cell->SetInternalClock(0);
      cell->SetPreviousPosition({x, y, z});
      cell->AddBiologyModule(new Substance_secretion_BM());
      cell->AddBiologyModule(new RGC_mosaic_BM());
      // cell.AddBiologyModule(Neurite_creation_BM());
      rm->push_back(cell);
    }
  }  // end CellCreator


  // RI computation
  inline double computeRI(vector<Double3> coordList) {
    if (coordList.size() < 2) {
      return 0;
    }
    vector<double> shortestDistList;
    for (unsigned int i = 0; i < coordList.size();
         i++) {  // for each cell of same type in the simulation
      Double3 cellPosition = coordList[i];

      vector<double> shortestDist;
      // for each other cell of same type in the simulation
      for (unsigned int j = 0; j < coordList.size(); j++) {
        Double3 otherCellPosition = coordList[j];

        // get the distance between those 2 cells (x-y plan only)
        double tempsDistance =
            sqrt(pow(cellPosition[0] - otherCellPosition[0], 2) +
                 pow(cellPosition[1] - otherCellPosition[1], 2));
        // if cell is closer and is not itself
        if (tempsDistance != 0) {
          shortestDist.push_back(tempsDistance);
        }
      }
      // save the shortest distance to neighbour
      shortestDistList.push_back(
        *min_element(shortestDist.begin(), shortestDist.end()));
    }
    // compute mean
    double temps_sum = 0;
    for (unsigned int i = 0; i < shortestDistList.size(); i++) {
      temps_sum += shortestDistList[i];
    }
    double aveShortestDist = temps_sum / (double)shortestDistList.size();
    // compute std
    double temp = 0;
    for (unsigned int i = 0; i < shortestDistList.size(); i++) {
      double val = shortestDistList[i];
      double squrDiffToMean = pow(val - aveShortestDist, 2);
      temp += squrDiffToMean;
    }
    double meanOfDiffs = temp / (double)(shortestDistList.size());
    double std = sqrt(meanOfDiffs);

    return aveShortestDist / std;  // return RI
  }  // end computeRI


  inline double getRI(int desiredCellType) {
    auto* sim = Simulation::GetActive();
    auto* rm = sim->GetResourceManager();
    vector<Double3> coordList; // list of coordinate
    rm->ApplyOnAllElements([&](SimObject* so, SoHandle) {
      auto* cell = dynamic_cast<MyCell*>(so);
      if (cell) {
        auto thisCellType = cell->GetCellType();
        if (thisCellType == desiredCellType) {  // if cell is of the desired type
          auto position = cell->GetPosition();  // get its position
          coordList.push_back(position);  // put cell coord in the list
        }
      }
    });  // end for cell in simulation
    cout << coordList.size() << " cells of type " << desiredCellType << endl;
    return computeRI(coordList);  // return RI for desired cell type
  }  // end getRI


  inline vector<array<double, 2>> getAllRI() {
    auto* sim = Simulation::GetActive();
    auto* rm = sim->GetResourceManager();
    vector<double> typesList;
    rm->ApplyOnAllElements([&](SimObject* so, SoHandle) {
      auto* cell = dynamic_cast<MyCell*>(so);
      if (cell) {
        auto thisCellType = cell->GetCellType();
        // if type is not -1 and not in list
        if (find(typesList.begin(), typesList.end(), thisCellType) == typesList.end()) {
          typesList.push_back(thisCellType); // put cell coord in the list
        }
      }
    });  // end for cell in simulation
    vector<array<double, 2>> listRi;
    for (unsigned int i = 0; i < typesList.size(); i++) {
      listRi.push_back({getRI(typesList[i]), typesList[i]});
    }
    return listRi;
  }

  inline double getDeathRate(int num_cells) {
    auto* sim = Simulation::GetActive();
    auto* rm = sim->GetResourceManager();
    int cellInSimu = 0;
    rm->ApplyOnAllElements([&](SimObject* so, SoHandle) {
      auto* cell = dynamic_cast<MyCell*>(so);
      if (cell) {
        cellInSimu++;
      }
    });  // end for cell in simulation

    return (1 - ((double)cellInSimu / num_cells)) * 100;
  } // end getDeathRate

}

#endif
