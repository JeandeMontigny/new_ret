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
#include "core/substance_initializers.h"

namespace bdm {

inline int Simulate(int argc, const char** argv) {
  int max_step = 3360; // 2240 = 14 days - 160 steps per day
  int cube_dim = 500; // 1000
  int cell_density = 500;
  int num_cells = cell_density*((double)cube_dim/1000)*((double)cube_dim/1000);
  double diffusion_coef = 0.5;
  double decay_const = 0.1;

  bool write_ri = false;
  bool write_positions = false;
  bool write_swc = true;
  bool clean_output_dir = true;

  auto set_param = [&](Param* param) {
    // Create an artificial bounds for the simulation space
    param->bound_space_ = true;
    param->min_bound_ = 0;
    param->max_bound_ = cube_dim + 300;
  };

  // initialise neuroscience modlues
  experimental::neuroscience::InitModule();

  Simulation simulation(argc, argv, set_param);
  // auto* rm = simulation.GetResourceManager();
  auto* scheduler = simulation.GetScheduler();
  auto* param = simulation.GetParam();
  auto* random = simulation.GetRandom();

  int my_seed = rand() % 10000;
  // my_seed = 9408;
  random->SetSeed(my_seed);
  cout << "Start simulation with " << cell_density
       << " cells/mm^2 using seed " << my_seed << endl;

  // create cells
  // CellCreator(param->min_bound_, param->max_bound_, num_cells, -1);
  // CellCreator(param->min_bound_, param->max_bound_, 31, 0);
  // CellCreator(param->min_bound_, param->max_bound_, 31, 1);
  // CellCreator(param->min_bound_, param->max_bound_, 31, 2);
  // CellCreator(param->min_bound_, param->max_bound_, 31, 3);
  // CellCreator(param->min_bound_, param->max_bound_, 62, 5);
  // CellCreator(param->min_bound_, param->max_bound_, 50, 203);
  CellCreator(param->min_bound_, param->max_bound_, 87, 204);
  CellCreator(param->min_bound_, param->max_bound_, 0, 101);
  CellCreator(param->min_bound_, param->max_bound_, 0, 201);

  // Order: substance_name, diffusion_coefficient, decay_constant, resolution
  // ModelInitializer::DefineSubstance(dg_200_, "off_aplhaa", diffusion_coef, decay_const, param->max_bound_/4);
  // ModelInitializer::DefineSubstance(dg_201_, "off_aplhab", diffusion_coef, decay_const, param->max_bound_/4);
  // ModelInitializer::DefineSubstance(dg_202_, "off_m1", diffusion_coef, decay_const, param->max_bound_/4);
  // ModelInitializer::DefineSubstance(dg_203_, "off_j", diffusion_coef, decay_const, param->max_bound_/4);
  // ModelInitializer::DefineSubstance(dg_204_, "off_mini_j", diffusion_coef, decay_const, param->max_bound_/4);
  // ModelInitializer::DefineSubstance(dg_205_, "off_midi_j", diffusion_coef, decay_const, param->max_bound_/4);
  // ModelInitializer::DefineSubstance(dg_206_, "off_u", diffusion_coef, decay_const, param->max_bound_/4);
  // ModelInitializer::DefineSubstance(dg_207_, "off_v", diffusion_coef, decay_const, param->max_bound_/4);
  // ModelInitializer::DefineSubstance(dg_208_, "off_w", diffusion_coef, decay_const, param->max_bound_/4);
  // ModelInitializer::DefineSubstance(dg_209_, "off_x", diffusion_coef, decay_const, param->max_bound_/4);
  // ModelInitializer::DefineSubstance(dg_210_, "off_y", diffusion_coef, decay_const, param->max_bound_/4);
  // ModelInitializer::DefineSubstance(dg_211_, "off_z", diffusion_coef, decay_const, param->max_bound_/4);


  // create substance for neurite attraction
  ModelInitializer::DefineSubstance(on_dendrites, "on_dendrites", 0, 0,
                                    param->max_bound_ / 2 );
  ModelInitializer::DefineSubstance(off_dendrites, "off_dendrites", 0, 0,
                                    param->max_bound_ / 2);
  // average peak distance for ON cells: 15.959 with std of 5.297;
  ModelInitializer::InitializeSubstance(on_dendrites, "on_dendrites",
      GaussianBand(43, 6, Axis::kZAxis));
  // average peak distance for OFF cells: 40.405 with std of 8.39;
  ModelInitializer::InitializeSubstance(off_dendrites, "off_dendrites",
      GaussianBand(67, 8, Axis::kZAxis));

  cout << "Cells created and substances initialised" << endl;

  //clean output dir
  if (clean_output_dir) {
    if (!system(Concat("[ -d ", param->output_dir_ ," ]").c_str())) {
      system(Concat("rm -rf ", param->output_dir_ ,"/results*").c_str());
      cout << "\'" << param->output_dir_ << "\'" << " folder cleaned" << endl;
    }
  }
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

  // Run simulation
  cout << "Simulating.." << endl;
  for (int i = 0; i < max_step/160; i++) {
    // if we want to export data from simulation
    if (write_ri || write_positions || write_swc) {
      for (int repet = 0; repet < 10; repet++) {
        scheduler->Simulate(16);
        int current_step = 16+(16*repet)+(160*i);

        if (write_ri) {
          vector<array<double, 2>> all_ri = GetAllRI();
          double death_rate = GetDeathRate(num_cells);
          for (unsigned int ri_i = 0; ri_i < all_ri.size(); ri_i++) {
            // step ri type death
            output_ri << current_step << " " << all_ri[ri_i][0]
                      << " " << all_ri[ri_i][1] << " " << death_rate << "\n";
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

   vector<array<double, 2>> all_ri = GetAllRI();
   double mean_ri = 0;
   for (unsigned int i = 0; i < all_ri.size(); i++) {
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

  if (write_swc) {
    WriteSwc(max_step, my_seed);
    std::cout << "Morphologies exported (swc files)" << std::endl;
  }

  cout << "Done" << endl;
  return 0;
} // end Simulate

}  // namespace bdm

#endif  // NEW_RET_H_
