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
#include "util_methods.h"

namespace bdm {

inline int Simulate(int argc, const char** argv) {
  int max_step = 70; // 2000 = 12 days - 160 steps per day
  int cube_dim = 1000; // 1000
  int cell_density = 986;
  int num_cells = cell_density*((double)cube_dim/1000)*((double)cube_dim/1000);
  double diffusion_coef = 0.5;
  double decay_const = 0.1;

  auto set_param = [&](Param* param) {
    // Create an artificial bounds for the simulation space
    param->bound_space_ = true;
    param->min_bound_ = 0;
    param->max_bound_ = cube_dim + 20;
    param->run_mechanical_interactions_ = true;
  };

  Simulation simulation(argc, argv, set_param);
  // auto* rm = simulation.GetResourceManager();
  auto* scheduler = simulation.GetScheduler();
  auto* param = simulation.GetParam();
  auto* random = simulation.GetRandom();

  int my_seed = rand() % 10000;
  // my_seed = 1142; // 9670
  random->SetSeed(my_seed);
  cout << "Start simulation with " << cell_density
       << " cells/mm^2 using seed " << my_seed << endl;

  // create cells
  CellCreator(param->min_bound_, param->max_bound_, num_cells, -1);

  // Order: substance_name, diffusion_coefficient, decay_constant, resolution
  ModelInitializer::DefineSubstance(dg_200_, "off_aplhaa", diffusion_coef, decay_const, param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_201_, "off_aplhab", diffusion_coef, decay_const, param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_202_, "off_m1", diffusion_coef, decay_const, param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_203_, "off_j", diffusion_coef, decay_const, param->max_bound_/4);
  // ModelInitializer::DefineSubstance(dg_204_, "off_mini_j", diffusion_coef, decay_const, param->max_bound_/4);
  // ModelInitializer::DefineSubstance(dg_205_, "off_midi_j", diffusion_coef, decay_const, param->max_bound_/4);
  // ModelInitializer::DefineSubstance(dg_206_, "off_u", diffusion_coef, decay_const, param->max_bound_/4);
  // ModelInitializer::DefineSubstance(dg_207_, "off_v", diffusion_coef, decay_const, param->max_bound_/4);
  // ModelInitializer::DefineSubstance(dg_208_, "off_w", diffusion_coef, decay_const, param->max_bound_/4);
  // ModelInitializer::DefineSubstance(dg_209_, "off_x", diffusion_coef, decay_const, param->max_bound_/4);
  // ModelInitializer::DefineSubstance(dg_210_, "off_y", diffusion_coef, decay_const, param->max_bound_/4);
  // ModelInitializer::DefineSubstance(dg_211_, "off_z", diffusion_coef, decay_const, param->max_bound_/4);

  cout << "Cells created and substances initialised" << endl;

  // Run simulation
  // for (int i = 0; i <= max_step/160; i++) {
  //   scheduler->Simulate(160);
  //   cout << setprecision(3) << "day " << i << "/" << (int)max_step/160 << ": "
  //        << getDeathRate(num_cells) << "% of cell death" << endl;
  //  vector<array<double, 2>> all_ri = getAllRI();
  //  for (unsigned int i = 0; i < all_ri.size(); i++) {
  //    cout << "RI = " << all_ri[i][0] << " for cell type " << all_ri[i][1] << endl;
  //  }
  // }
  for (int i = 0; i <= max_step/10; i++) {
    scheduler->Simulate(10);
    cout << setprecision(3) << "step " << i*10 << "/" << (int)max_step<< endl;
     vector<array<double, 2>> all_ri = getAllRI();
     for (unsigned int i = 0; i < all_ri.size(); i++) {
       cout << "RI = " << all_ri[i][0] << " for cell type " << all_ri[i][1] << endl;
     }
  }

  return 0;
} // end Simulate

}  // namespace bdm

#endif  // NEW_RET_H_
