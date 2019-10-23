#ifndef RGC_SOMA_BM_
#define RGC_SOMA_BM_

#include "biodynamo.h"
#include "extended_objects.h"

namespace bdm {

  enum Substances { dg_0_ };

  // Define cell behavior for mosaic formation
  struct RGC_mosaic_BM : public BaseBiologyModule {
    BDM_STATELESS_BM_HEADER(RGC_mosaic_BM, BaseBiologyModule, 1);

  public:
    RGC_mosaic_BM() : BaseBiologyModule(gAllEventIds) {}

    void Run(SimObject* so) override {
      if (auto* cell = dynamic_cast<MyCell*>(so)) {
        auto* sim = Simulation::GetActive();
        auto* rm = sim->GetResourceManager();
        auto* random = sim->GetRandom();

        auto& position = cell->GetPosition();
        int cellClock = cell->GetInternalClock();
        double concentration = 0;
        Double3 gradient, diff_gradient, gradient_z;
        DiffusionGrid* dg = nullptr;

        // if cell is type 0, concentration and gradient are substance 0
        if (cell->GetCellType() == 0) {
          dg = rm->GetDiffusionGrid(dg_0_);
          dg->GetGradient(position, &gradient);
          gradient_z = gradient * 0.2;
          gradient_z[0] = 0; gradient_z[1] = 0;
          diff_gradient = gradient * -0.1;
          diff_gradient[2] = 0;
          concentration = dg->GetConcentration(position);
        }

        bool withMovement = true;
        double movementThreshold = 2.0; // 03707 // 0371 - 03705
        bool withDeath = true;
        double deathThreshold = 3.35; // 1.118 // 0532 // 0524 - 0523
        // 1.155 // 1.154

        if (cellClock < 1900 && cellClock%3==0 ) {
          // // add small random movements
          cell->UpdatePosition(
              {random->Uniform(-0.01, 0.01), random->Uniform(-0.01, 0.01), 0});
          // cell growth
          if (cell->GetDiameter() < 14 && random->Uniform(0, 1) < 0.02) {
            cell->ChangeVolume(2000);
          }
          // add vertical migration as the multi layer colapse in just on layer
          // cell->UpdatePosition(gradient_z);
        }

        /* -- cell movement -- */
        if (cellClock > 1200) {movementThreshold=1.9;}
        if (withMovement && cellClock >= 100 && cellClock < 2800
          && concentration >= movementThreshold && cellClock%3==0) {
            // cell movement based on homotype substance gradient
            cell->UpdatePosition(diff_gradient);
            // update distance travelled by this cell
            auto previousPosition = cell->GetPreviousPosition();
            auto currentPosition = cell->GetPosition();
            cell->SetDistanceTravelled(cell->GetDistanceTravelled() +
            (sqrt(pow(currentPosition[0] - previousPosition[0], 2) +
            pow(currentPosition[1] - previousPosition[1], 2))));
            cell->SetPreviousPosition(cell->GetPosition());
          }  // end tangential migration

          /* -- cell death -- */
          if (withDeath && cellClock >= 100 && cellClock < 1200
            && cellClock%4==0) {
              // add vertical migration as the multi layer colapse in just on layer
              cell->UpdatePosition(gradient_z);
              // cell death depending on homotype substance concentration
              if (concentration > deathThreshold
                  && random->Uniform(0, 1) < 0.1) { // 0.25
                cell->RemoveFromSimulation();
              }
            } // end cell death

            /* -- internal clock -- */
            // probability to increase internal clock
            if (random->Uniform(0, 1) < 0.96) {
              // update cell internal clock
              cell->SetInternalClock(cell->GetInternalClock() + 1);
            } // end update cell internal clock

          }
        } // end Run()
  }; // end biologyModule RGC_mosaic_BM


  // Define cell behavior for substance secretion
  struct Substance_secretion_BM : public BaseBiologyModule {
    BDM_STATELESS_BM_HEADER(Substance_secretion_BM, BaseBiologyModule, 1);

  public:
    Substance_secretion_BM() : BaseBiologyModule(gAllEventIds) {}

    void Run(SimObject* so) override {
      if (auto* cell = dynamic_cast<MyCell*>(so)) {
        auto* rm = Simulation::GetActive()->GetResourceManager();

        DiffusionGrid* dg = nullptr;
        if (cell->GetCellType() == 0) {
          dg = rm->GetDiffusionGrid(dg_0_);
        }
        int interval = 4;
        // if (cell->GetInternalClock() > 1200) {interval = 2;}
        if (cell->GetInternalClock()%interval==0) {
          auto& secretion_position = cell->GetPosition();
          dg->IncreaseConcentrationBy(secretion_position, 1);
          }
        }
      } // end Run()
  }; // end biologyModule Substance_secretion_BM


// cell->RemoveBiologyModule(this);

} // namespace bdm

#endif
