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
          if (position[2]>27) {gradient_z={0, 0, -0.01};}
          else {gradient_z={0, 0, 0.01};}
          // gradient_z = gradient * 0.2;
          // gradient_z[0] = 0; gradient_z[1] = 0;
          diff_gradient = gradient * -0.1;
          diff_gradient[2] = 0;
          concentration = dg->GetConcentration(position);
        }

        bool withMovement = true;
        double movementThreshold = 1.735;
        bool withDeath = true;
        double deathThreshold = 1.78;

        // deathThreshold depending on initial density to obtain ~65% death rate
        if (false) {
          // 1000 (350 final) -- RI ~6-7
          movementThreshold = 1.75;
          deathThreshold = 1.82;
          // 800 (280 final) -- RI ~5-7
          movementThreshold = 1.735;
          deathThreshold = 1.78;
          // 600 (210 final) -- RI ~4-7
          movementThreshold = 1.73;
          deathThreshold = 1.8;
          // 400 (140 final) -- RI ~4-5
          movementThreshold = 1.72;
          deathThreshold = 1.8;
          // 200 (70 final) -- RI ~2.5-4
          movementThreshold = 1.71;
          deathThreshold = 1.81;
          // 100 (35 final) -- RI ~1.8-3
          movementThreshold = 1.65;
          deathThreshold = 1.85;
          // 60 (20 final) -- RI ~1.8-2.5
          movementThreshold = 1.66;
          deathThreshold = 1.86;
        }

        if (cellClock < 960 && cellClock%3==0) {
          // // add small random movements
          cell->UpdatePosition(
              {random->Uniform(-0.01, 0.01), random->Uniform(-0.01, 0.01), 0});
          // cell growth
          if (cell->GetDiameter() < 14 && random->Uniform(0, 1) < 0.02) {
            cell->ChangeVolume(3000);
          }
          // add vertical migration as the multi layer colapse in just on layer
          // cell->UpdatePosition(gradient_z);
        }

        /* -- cell movement -- */
        if (withMovement && cellClock >= 100 && cellClock < 1920
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
          if (withDeath && cellClock >= 100 && cellClock < 960
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
        if (cell->GetInternalClock()%3==0) {
          auto& secretion_position = cell->GetPosition();
          dg->IncreaseConcentrationBy(secretion_position, 1);
          }
        }
      } // end Run()
  }; // end biologyModule Substance_secretion_BM


// cell->RemoveBiologyModule(this);

} // namespace bdm

#endif
