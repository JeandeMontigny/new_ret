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
  bool clean_result_dir = false;

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
    // on-off
    CellCreator(param->min_bound_, param->max_bound_, 357, 0);
    CellCreator(param->min_bound_, param->max_bound_, 357, 1);
    CellCreator(param->min_bound_, param->max_bound_, 357, 2);
    CellCreator(param->min_bound_, param->max_bound_, 357, 3);
    CellCreator(param->min_bound_, param->max_bound_, 57, 4);
    CellCreator(param->min_bound_, param->max_bound_, 714, 5);
    CellCreator(param->min_bound_, param->max_bound_, 57, 6);
    CellCreator(param->min_bound_, param->max_bound_, 57, 7);
    CellCreator(param->min_bound_, param->max_bound_, 171, 8);
    CellCreator(param->min_bound_, param->max_bound_, 143, 9);
    CellCreator(param->min_bound_, param->max_bound_, 114, 10);
    CellCreator(param->min_bound_, param->max_bound_, 114, 11);
    // on
    CellCreator(param->min_bound_, param->max_bound_, 114, 100);
    CellCreator(param->min_bound_, param->max_bound_, 114, 101);
    CellCreator(param->min_bound_, param->max_bound_, 114, 102);
    CellCreator(param->min_bound_, param->max_bound_, 114, 103);
    CellCreator(param->min_bound_, param->max_bound_, 160, 104);
    CellCreator(param->min_bound_, param->max_bound_, 57, 105);
    CellCreator(param->min_bound_, param->max_bound_, 57, 106);
    CellCreator(param->min_bound_, param->max_bound_, 428, 107);
    CellCreator(param->min_bound_, param->max_bound_, 286, 108);
    CellCreator(param->min_bound_, param->max_bound_, 286, 109);
    CellCreator(param->min_bound_, param->max_bound_, 228, 110);
    CellCreator(param->min_bound_, param->max_bound_, 171, 111);
    CellCreator(param->min_bound_, param->max_bound_, 171, 112);
    CellCreator(param->min_bound_, param->max_bound_, 143, 113);
    CellCreator(param->min_bound_, param->max_bound_, 143, 114);
    CellCreator(param->min_bound_, param->max_bound_, 97, 115);
    CellCreator(param->min_bound_, param->max_bound_, 57, 116);
    CellCreator(param->min_bound_, param->max_bound_, 57, 117);
    CellCreator(param->min_bound_, param->max_bound_, 57, 118);
    // off
    CellCreator(param->min_bound_, param->max_bound_, 114, 200);
    CellCreator(param->min_bound_, param->max_bound_, 114, 201);
    CellCreator(param->min_bound_, param->max_bound_, 180, 202);
    CellCreator(param->min_bound_, param->max_bound_, 571, 203);
    CellCreator(param->min_bound_, param->max_bound_, 1000, 204);
    CellCreator(param->min_bound_, param->max_bound_, 228, 205);
    CellCreator(param->min_bound_, param->max_bound_, 57, 206);
    CellCreator(param->min_bound_, param->max_bound_, 57, 207);
    CellCreator(param->min_bound_, param->max_bound_, 171, 208);
    CellCreator(param->min_bound_, param->max_bound_, 143, 209);
    CellCreator(param->min_bound_, param->max_bound_, 114, 210);
    CellCreator(param->min_bound_, param->max_bound_, 106, 211);
  }
  // Order: substance_name, diffusion_coefficient, decay_constant, resolution

  // on-off
  // 125, 125, 125, 125, 20, 250, 20, 20, 60, 50, 40, 40
  ModelInitializer::DefineSubstance(dg_000_, "on-off_dsgca", diffusion_coef,
    decay_const, resolution);
  ModelInitializer::DefineSubstance(dg_001_, "on-off_dsgcb", diffusion_coef,
    decay_const, resolution);
  ModelInitializer::DefineSubstance(dg_002_, "on-off_dsgcc", diffusion_coef,
    decay_const, resolution);
  ModelInitializer::DefineSubstance(dg_003_, "on-off_dsgcd", diffusion_coef,
    decay_const, resolution);
  ModelInitializer::DefineSubstance(dg_004_, "on-off_m3", diffusion_coef,
    decay_const, resolution);
  ModelInitializer::DefineSubstance(dg_005_, "on-off_led", diffusion_coef,
    decay_const, resolution);
  ModelInitializer::DefineSubstance(dg_006_, "on-off_u", diffusion_coef,
    decay_const, resolution);
  ModelInitializer::DefineSubstance(dg_007_, "on-off_v", diffusion_coef,
    decay_const, resolution);
  ModelInitializer::DefineSubstance(dg_008_, "on-off_w", diffusion_coef,
    decay_const, resolution);
  ModelInitializer::DefineSubstance(dg_009_, "on-off_x", diffusion_coef,
    decay_const, resolution);
  ModelInitializer::DefineSubstance(dg_010_, "on-off_y", diffusion_coef,
    decay_const, resolution);
  ModelInitializer::DefineSubstance(dg_011_, "on-off_z", diffusion_coef,
    decay_const, resolution);
  // on
  // 40, 40, 40, 40, 56, 20, 20, 150, 100, 100, 80, 60, 60, 50, 50, 34, 20, 20, 20
  ModelInitializer::DefineSubstance(dg_100_, "on_dsgca", diffusion_coef,
    decay_const, resolution);
  ModelInitializer::DefineSubstance(dg_101_, "on_dsgcb", diffusion_coef,
    decay_const, resolution);
  ModelInitializer::DefineSubstance(dg_102_, "on_dsgcc", diffusion_coef,
    decay_const, resolution);
  ModelInitializer::DefineSubstance(dg_103_, "on_aplha", diffusion_coef,
    decay_const, resolution);
  ModelInitializer::DefineSubstance(dg_104_, "on_m2", diffusion_coef, decay_const,
    resolution);
  ModelInitializer::DefineSubstance(dg_105_, "on_m4", diffusion_coef, decay_const,
    resolution);
  ModelInitializer::DefineSubstance(dg_106_, "on_m5", diffusion_coef, decay_const,
    resolution);
  ModelInitializer::DefineSubstance(dg_107_, "on_o", diffusion_coef, decay_const,
    resolution);
  ModelInitializer::DefineSubstance(dg_108_, "on_p", diffusion_coef, decay_const,
    resolution);
  ModelInitializer::DefineSubstance(dg_109_, "on_q", diffusion_coef, decay_const,
    resolution);
  ModelInitializer::DefineSubstance(dg_110_, "on_r", diffusion_coef, decay_const,
    resolution);
  ModelInitializer::DefineSubstance(dg_111_, "on_s", diffusion_coef, decay_const,
    resolution);
  ModelInitializer::DefineSubstance(dg_112_, "on_t", diffusion_coef, decay_const,
    resolution);
  ModelInitializer::DefineSubstance(dg_113_, "on_u", diffusion_coef, decay_const,
    resolution);
  ModelInitializer::DefineSubstance(dg_114_, "on_v", diffusion_coef, decay_const,
    resolution);
  ModelInitializer::DefineSubstance(dg_115_, "on_w", diffusion_coef, decay_const,
    resolution);
  ModelInitializer::DefineSubstance(dg_116_, "on_x", diffusion_coef, decay_const,
    resolution);
  ModelInitializer::DefineSubstance(dg_117_, "on_y", diffusion_coef, decay_const,
    resolution);
  ModelInitializer::DefineSubstance(dg_118_, "on_z", diffusion_coef, decay_const,
    resolution);
  // off
  // 40, 40, 63, 200, 350, 80, 20, 20, 60, 50, 40, 37
  ModelInitializer::DefineSubstance(dg_200_, "off_aplhaa", diffusion_coef,
    decay_const, resolution);
  ModelInitializer::DefineSubstance(dg_201_, "off_aplhab", diffusion_coef,
    decay_const, resolution);
  ModelInitializer::DefineSubstance(dg_202_, "off_m1", diffusion_coef, decay_const,
    resolution);
  ModelInitializer::DefineSubstance(dg_203_, "off_j", diffusion_coef, decay_const,
    resolution);
  ModelInitializer::DefineSubstance(dg_204_, "off_mini_j", diffusion_coef,
    decay_const, resolution);
  ModelInitializer::DefineSubstance(dg_205_, "off_midi_j", diffusion_coef,
    decay_const, resolution);
  ModelInitializer::DefineSubstance(dg_206_, "off_u", diffusion_coef, decay_const,
    resolution);
  ModelInitializer::DefineSubstance(dg_207_, "off_v", diffusion_coef, decay_const,
    resolution);
  ModelInitializer::DefineSubstance(dg_208_, "off_w", diffusion_coef, decay_const,
    resolution);
  ModelInitializer::DefineSubstance(dg_209_, "off_x", diffusion_coef, decay_const,
    resolution);
  ModelInitializer::DefineSubstance(dg_210_, "off_y", diffusion_coef, decay_const,
    resolution);
  ModelInitializer::DefineSubstance(dg_211_, "off_z", diffusion_coef, decay_const,
    resolution);

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

  WritePositions(0, my_seed);

  scheduler->Simulate(1);
  AnkurDeath();
  scheduler->Simulate(1);

  WritePositions(1, my_seed);
  
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

  WritePositions(10, my_seed);

  cout << "Done" << endl;
  return 0;
} // end Simulate

}  // namespace bdm

#endif  // NEW_RET_H_
