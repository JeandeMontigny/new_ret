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
  int maxStep = 2900;
  int cubeDim = 500;
  int num_cells = 80; // x4 to have c/mm2 density
  double diffusion_coef = 0.65;
  double decay_const = 0.1;

  double cellDensity = (double)num_cells * 1e6 / (cubeDim * cubeDim);
  cout << "cell density: " << cellDensity << " cells per mm^2" << endl;

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
  // mySeed = 2089; // 9670
  random->SetSeed(mySeed);
  cout << "modelling with seed " << mySeed << endl;

  // create cells
  CellCreator(0, 500, num_cells, 0);

  // Order: substance_name, diffusion_coefficient, decay_constant, resolution
  ModelInitializer::DefineSubstance(dg_0_, "on", diffusion_coef, decay_const,
                                    param->max_bound_/2);

  // Run simulation for one timestep
  for (int i = 0; i <= maxStep/240; i++) {
    scheduler->Simulate(240);
    cout << "day " << i << "/" << (int)maxStep/240
         << " ; " << getDeathRate(num_cells) << "% of cell death"
         << " ; RI = " << getRI(0) << endl;
  }

  return 0;
}

}  // namespace bdm

#endif  // NEW_RET_H_
