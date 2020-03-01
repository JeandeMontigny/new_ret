// -----------------------------------------------------------------------------
//
// Copyright (C) The BioDynaMo Project.
// All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------
#ifndef NEW_RET_H_
#define NEW_RET_H_

#include "biodynamo.h"
#include "extended_objects.h"
#include "util_methods.h"

namespace bdm {

inline int Simulate(int argc, const char** argv) {
  bool cell_fate = false;

  bool verbose = false;
  bool write_ri = true;
  bool write_positions = true;
  bool write_swc = false;
  bool write_distance = true;
  bool clean_result_dir = true;

  int max_step = 2240; // 2240 = 13 days - 160 steps per day
  int cube_dim = 1000; // 1000
  int cell_density = 8600 ; // 8600 - 65% ~= 3000
  int num_cells = cell_density*((double)cube_dim/1000)*((double)cube_dim/1000);
  int grid_spacing = 4;

  auto set_param = [&](Param* param) {
    // Create an artificial bounds for the simulation space
    param->bound_space_ = true;
    param->min_bound_ = 0;
    param->max_bound_ = cube_dim + 300;
    param->run_mechanical_interactions_ = true;
  };

  // initialise neuroscience modlues
  experimental::neuroscience::InitModule();

  Simulation simulation(argc, argv, set_param);
  // auto* rm = simulation.GetResourceManager();
  auto* scheduler = simulation.GetScheduler();
  auto* param = simulation.GetParam();
  auto* random = simulation.GetRandom();

  int resolution = param->max_bound_/grid_spacing;
  // ratio diffusion_coef/spacing/spacing = 0.125
  double diffusion_coef = 0.125*grid_spacing*grid_spacing;
  double decay_const = 0;

  int my_seed = rand() % 10000;
  // my_seed = 56;
  random->SetSeed(my_seed);
  cout << "Start simulation with " << cell_density
       << " cells/mm^2 using seed " << my_seed << endl;

  // with cell fate
  if (cell_fate) {
    CellCreator(param->min_bound_, param->max_bound_, num_cells, -1);
  }
  // without cell fate
  else {
    CellCreator(param->min_bound_, param->max_bound_, 357, 0); // GCL
    CellCreator(param->min_bound_, param->max_bound_, 357, 1); // INL
  }

  // Order: substance_name, diffusion_coefficient, decay_constant, resolution
  ModelInitializer::DefineSubstance(dg_0_, "sac-gcl", diffusion_coef,
    decay_const, resolution);
  ModelInitializer::DefineSubstance(dg_1_, "sac-inl", diffusion_coef,
    decay_const, resolution);

  cout << "Cells created and substances initialised" << endl;

  // prepare export
  ofstream output_ri;
  if ((write_ri || write_positions || write_swc) && system(
    Concat("mkdir -p ", param->output_dir_,
           "/results", my_seed).c_str())) {
      cout << "error during " << param->output_dir_
           << "/results folder creation" << endl;
  }
  if (write_ri) {
    output_ri.open(Concat(param->output_dir_, "/results", my_seed,
                          "/RI_" + to_string(my_seed) + ".txt"));
  }
  if (write_positions && system(
    Concat("mkdir -p ", param->output_dir_,
           "/results", my_seed, "/cells_position").c_str())) {
      cout << "error during " << param->output_dir_
           << "/results"<< my_seed <<"cells_position folder creation" << endl;
  }
  if (write_swc && system(
    Concat("mkdir -p ", param->output_dir_,
           "/results", my_seed, "/swc_files").c_str())) {
      cout << "error during " << param->output_dir_
           << "/results"<< my_seed <<"/swc_files folder creation" << endl;
  }


  int max_pop_list [212] = { };

  // Run simulation
  cout << "Simulating.." << endl;
  for (int i = 0; i < max_step/160; i++) {
    // if we want to export data from simulation
    if (write_ri || write_positions || write_swc) {
      for (int repet = 0; repet < 10; repet++) {
        scheduler->Simulate(16);
        int current_step = 16+(16*repet)+(160*i);

        if (write_ri) {
          vector<array<double, 3>> all_ri = GetAllRI();
          double death_rate = GetDeathRate(num_cells);
          for (unsigned int ri_i = 0; ri_i < all_ri.size(); ri_i++) {
            // step ri type death
            output_ri << current_step << " " << all_ri[ri_i][0]
                      << " " << all_ri[ri_i][1]
	              << " " << all_ri[ri_i][2]
		      << " " << death_rate << "\n";

	    if (max_pop_list[(int)all_ri[ri_i][1]] < all_ri[ri_i][2]) {
	      max_pop_list[(int)all_ri[ri_i][1]] = all_ri[ri_i][2];
	    }

          }
        }

        if (write_positions) {
          WritePositions(current_step, my_seed);
        }
        if (false && write_swc) {
          WriteSwc(current_step, my_seed);
        }
      } // for step up to 160
    } // if export data

    else {
      scheduler->Simulate(160);
    }

    vector<array<double, 3>> all_ri = GetAllRI();
    double mean_ri = 0;
    for (unsigned int i = 0; i < all_ri.size(); i++) {
      if (verbose) {
        cout << "type: " << all_ri[i][1]
  	         << " - ri: " << all_ri[i][0]
  	         << " - population: " << all_ri[i][2]
  	         << " - max pop: " << max_pop_list[(int)all_ri[i][1]]
  	         << " - death: " << (1- ((double)all_ri[i][2] /
                max_pop_list[(int)all_ri[i][1]])) *100 << endl;
      }
      mean_ri += all_ri[i][0];
    }
    cout << setprecision(3)
	 << "Day " << i+1 << "/" << (int)max_step/160 << " simulated:\n"
	 << "Average ri = " << (double)mean_ri/all_ri.size() << " ; "
	 << GetDeathRate(num_cells) << "% of cell death"<< endl;
    //TODO delete all "mosaic" substances in simulation after mosaics are done
    // if (i > 2100) {
    //   delete [substances];
    // }
  }

  if (write_distance) {
    ExportMigrationDitance(my_seed);
    std::cout << "Migration distance exported" << std::endl;
  }

  if (write_swc) {
    WriteSwc(max_step, my_seed);
    std::cout << "Morphologies exported (swc files)" << std::endl;
  }

  cout << "Done" << endl;
  return 0;
} // end Simulate

}  // namespace bdm

#endif  // NEW_RET_H_
