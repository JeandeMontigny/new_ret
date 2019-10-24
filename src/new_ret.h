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
  int maxStep = 2000; // 12 days - 160 steps per day
  int cubeDim = 1000;
  int cell_density = 800;
  int num_cells = cell_density*((double)cubeDim/1000)*((double)cubeDim/1000);
  double diffusion_coef = 0.5;
  double decay_const = 0.1;

  auto set_param = [&](Param* param) {
    // Create an artificial bounds for the simulation space
    param->bound_space_ = true;
    param->min_bound_ = 0;
    param->max_bound_ = cubeDim + 20;
    param->run_mechanical_interactions_ = true;
  };

  Simulation simulation(argc, argv, set_param);
  // auto* rm = simulation.GetResourceManager();
  auto* scheduler = simulation.GetScheduler();
  auto* param = simulation.GetParam();
  auto* random = simulation.GetRandom();

  int mySeed = rand() % 10000;
  // mySeed = 1142; // 9670
  random->SetSeed(mySeed);
  cout << "Start simulation with " << cell_density
       << " cells/mm^2 using seed " << mySeed << endl;

  // create cells
  CellCreator(param->min_bound_, param->max_bound_, num_cells, 0);

  // Order: substance_name, diffusion_coefficient, decay_constant, resolution
  ModelInitializer::DefineSubstance(dg_0_, "on", diffusion_coef, decay_const,
                                    param->max_bound_/4);

  cout << "Cells created and substances initialised" << endl;

  // Run simulation
  for (int i = 0; i <= maxStep/160; i++) {
    scheduler->Simulate(160);
    cout << setprecision(3) << "day " << i << "/" << (int)maxStep/160 << ": "
         << getDeathRate(num_cells) << "% of cell death\tRI = " << getRI(0) << endl;
  }

  return 0;
} // end Simulate

}  // namespace bdm

#endif  // NEW_RET_H_
