#ifndef RGC_SOMA_BM_
#define RGC_SOMA_BM_

#include "biodynamo.h"
#include "extended_objects.h"

namespace bdm {

  // enumerate substances in simulation
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

        // use corresponding diffusion grid
        if (cell->GetCellType() == 0) {
          dg = rm->GetDiffusionGrid(dg_0_);
          dg->GetGradient(position, &gradient);
          concentration = dg->GetConcentration(position);
          if (position[2]>27) {gradient_z={0, 0, -0.01};}
          else {gradient_z={0, 0, 0.01};}
          diff_gradient = gradient * -0.1; diff_gradient[2] = 0;
        }

        bool withMovement = true;
        double movementThreshold = 1.735;
        bool withDeath = true;
        double deathThreshold = 1.76;

        // thresholds depending on initial density to obtain ~65% death rate
        if (false) {
          // 1000 (350 final) -- RI ~6-7
          movementThreshold = 1.745;
          deathThreshold = 1.765;
          // 800 (280 final) -- RI ~5-7
          movementThreshold = 1.735;
          deathThreshold = 1.76;
          // 600 (210 final) -- RI ~5-7
          movementThreshold = 1.73;
          deathThreshold = 1.77;
          // 400 (140 final) -- RI ~4-5
          movementThreshold = 1.72;
          deathThreshold = 1.775;
          // 200 (70 final) -- RI ~2.5-3
          movementThreshold = 1.71;
          deathThreshold = 1.78;
          // 100 (35 final) -- RI ~2-3
          movementThreshold = 1.7;
          deathThreshold = 1.79;
          // 60 (20 final) -- RI ~1.8-2.5
          movementThreshold = 1.7;
          deathThreshold = 1.78;
        }

        /* -- cell growth -- */
        if (cellClock < 960 && cellClock%3==0) {
          // // add small random movements
          cell->UpdatePosition(
              {random->Uniform(-0.01, 0.01), random->Uniform(-0.01, 0.01), 0});
          // cell growth
          if (cell->GetDiameter() < 14 && random->Uniform(0, 1) < 0.02) {
            cell->ChangeVolume(3000);
          }
          // layer colapse if no cell death
          if (!withDeath) {
            cell->UpdatePosition(gradient_z);
          }
        } // end cell growth

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

          } // end if MyCell
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
      } // end if MyCell
    } // end Run()
  }; // end biologyModule Substance_secretion_BM


// cell->RemoveBiologyModule(this);

} // namespace bdm

#endif
