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
  int max_step = 2240; // 2240 = 13 days - 160 steps per day
  int cube_dim = 1000; // 1000
  int cell_density = 8600;
  int num_cells = cell_density*((double)cube_dim/1000)*((double)cube_dim/1000);
  double diffusion_coef = 0.5;
  double decay_const = 0.1;

  bool write_ri = true;
  bool write_positions = true;
  bool write_swc = false;
  bool clean_result_dir = true;

  auto set_param = [&](Param* param) {
    // Create an artificial bounds for the simulation space
    param->bound_space_ = true;
    param->min_bound_ = 0;
    param->max_bound_ = cube_dim + 100;
    param->run_mechanical_interactions_ = true;
  };

  // initialise neuroscience modlues
  experimental::neuroscience::InitModule();

  Simulation simulation(argc, argv, set_param);
  // auto* rm = simulation.GetResourceManager();
  auto* scheduler = simulation.GetScheduler();
  auto* param = simulation.GetParam();
  auto* random = simulation.GetRandom();

  int my_seed = rand() % 10000;
  // my_seed = 56;
  random->SetSeed(my_seed);
  cout << "Start simulation with " << cell_density
       << " cells/mm^2 using seed " << my_seed << endl;

  // create cells
  CellCreator(param->min_bound_, param->max_bound_, num_cells, -1);

  // Order: substance_name, diffusion_coefficient, decay_constant, resolution

  // on-off
  // 125, 125, 125, 125, 20, 250, 20, 20, 60, 50, 40, 40
  ModelInitializer::DefineSubstance(dg_000_, "on-off_dsgca", diffusion_coef,
    decay_const, param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_001_, "on-off_dsgcb", diffusion_coef,
    decay_const, param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_002_, "on-off_dsgcc", diffusion_coef,
    decay_const, param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_003_, "on-off_dsgcd", diffusion_coef,
    decay_const, param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_004_, "on-off_m3", diffusion_coef,
    decay_const, param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_005_, "on-off_led", diffusion_coef,
    decay_const, param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_006_, "on-off_u", diffusion_coef,
    decay_const, param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_007_, "on-off_v", diffusion_coef,
    decay_const, param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_008_, "on-off_w", diffusion_coef,
    decay_const, param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_009_, "on-off_x", diffusion_coef,
    decay_const, param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_010_, "on-off_y", diffusion_coef,
    decay_const, param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_011_, "on-off_z", diffusion_coef,
    decay_const, param->max_bound_/4);
  // on
  // 40, 40, 40, 40, 56, 20, 20, 150, 100, 100, 80, 60, 60, 50, 50, 34, 20, 20, 20
  ModelInitializer::DefineSubstance(dg_100_, "on_dsgca", diffusion_coef,
    decay_const, param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_101_, "on_dsgcb", diffusion_coef,
    decay_const, param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_102_, "on_dsgcc", diffusion_coef,
    decay_const, param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_103_, "on_aplha", diffusion_coef,
    decay_const, param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_104_, "on_m2", diffusion_coef, decay_const,
    param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_105_, "on_m4", diffusion_coef, decay_const,
    param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_106_, "on_m5", diffusion_coef, decay_const,
    param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_107_, "on_o", diffusion_coef, decay_const,
    param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_108_, "on_p", diffusion_coef, decay_const,
    param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_109_, "on_q", diffusion_coef, decay_const,
    param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_110_, "on_r", diffusion_coef, decay_const,
    param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_111_, "on_s", diffusion_coef, decay_const,
    param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_112_, "on_t", diffusion_coef, decay_const,
    param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_113_, "on_u", diffusion_coef, decay_const,
    param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_114_, "on_v", diffusion_coef, decay_const,
    param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_115_, "on_w", diffusion_coef, decay_const,
    param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_116_, "on_x", diffusion_coef, decay_const,
    param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_117_, "on_y", diffusion_coef, decay_const,
    param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_118_, "on_z", diffusion_coef, decay_const,
    param->max_bound_/4);
  // off
  // 40, 40, 63, 200, 350, 80, 20, 20, 60, 50, 40, 37
  ModelInitializer::DefineSubstance(dg_200_, "off_aplhaa", diffusion_coef,
    decay_const, param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_201_, "off_aplhab", diffusion_coef,
    decay_const, param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_202_, "off_m1", diffusion_coef, decay_const,
    param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_203_, "off_j", diffusion_coef, decay_const,
    param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_204_, "off_mini_j", diffusion_coef,
    decay_const, param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_205_, "off_midi_j", diffusion_coef,
    decay_const, param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_206_, "off_u", diffusion_coef, decay_const,
    param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_207_, "off_v", diffusion_coef, decay_const,
    param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_208_, "off_w", diffusion_coef, decay_const,
    param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_209_, "off_x", diffusion_coef, decay_const,
    param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_210_, "off_y", diffusion_coef, decay_const,
    param->max_bound_/4);
  ModelInitializer::DefineSubstance(dg_211_, "off_z", diffusion_coef, decay_const,
    param->max_bound_/4);

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
