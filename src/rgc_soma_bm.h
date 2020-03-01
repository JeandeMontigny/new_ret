#ifndef RGC_SOMA_BM_
#define RGC_SOMA_BM_

#include "biodynamo.h"
#include "extended_objects.h"
#include "rgc_dendrite_bm.h"

namespace bdm {

  // enumerate substances in simulation
  enum Substances { dg_0_, dg_1_ };

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
        int cell_clock = cell->GetInternalClock();
        int cell_type = cell->GetCellType();
        double concentration = 0;
        Double3 gradient_gcl, gradient_inl, diff_gradient, gradient_z;
        DiffusionGrid* dg = nullptr;

        bool with_fate = false;
        bool with_movement = true;
        double movement_threshold = 1.735;
        bool with_death = true;
        double death_threshold = 1.76;

	dg_gcl = rm->GetDiffusionGrid("sac-gcl");
	double concentration_gcl = dg_gcl->GetConcentration(position);
        dg_gcl->GetGradient(position, &gradient_gcl);

	dg_inl = rm->GetDiffusionGrid("sac-inl");
	double concentration_inl = dg_inl->GetConcentration(position);
        dg_inl->GetGradient(position, &gradient_inl);


        /* -- cell fate -- */
        if (cell_type == -1) {
          if (cell_clock%2!=0 || random->Uniform(0, 1) < 0.9) { return; }
	  // if no significant substances concentration, return
	  if (concentration_gcl < 1e-5 || concentration_inl < 1e-5) { return; }
	  // choose type
	  if (concentration_gcl < concentration_inl) {
	    cell->SetCellType(0);
	  }
	  else {
	    cell->SetCellType(1);
	  }

          return;
        } // end cell fate


        /* -- cell growth -- */
        if (cell_clock >= 100 && cell_clock < 1060 && cell_clock%3==0) {
          // add small random movements
          cell->UpdatePosition(
              {random->Uniform(-0.01, 0.01), random->Uniform(-0.01, 0.01), 0});
          // cell growth
          if (cell->GetDiameter() < 12 && random->Uniform(0, 1) < 0.02) {
            cell->ChangeVolume(2500);
          }
        } // end cell growth


        /* -- cell movement -- */
        if (with_movement && cell_clock >= 200 && cell_clock < 2020
          && concentration >= movement_threshold && cell_clock%3==0) {
            // cell movement based on homotype substance gradient
            cell->UpdatePosition(diff_gradient);
            // update distance travelled by this cell
            auto previous_position = cell->GetPreviousPosition();
            auto current_position = cell->GetPosition();
            cell->SetDistanceTravelled(cell->GetDistanceTravelled() +
            (sqrt(pow(current_position[0] - previous_position[0], 2) +
            pow(current_position[1] - previous_position[1], 2))));
            cell->SetPreviousPosition(cell->GetPosition());
          }  // end tangential migration

          /* -- cell death -- */
          if (with_death && cell_clock >= 200 && cell_clock < 840
            && cell_clock%4==0) {
	    // add vertical migration as the multi layer colapse in just on layer
	    cell->UpdatePosition(gradient_z);
	    // cell death depending on homotype substance concentration
	    if (concentration > death_threshold
		&& random->Uniform(0, 1) < 0.05) { // 0.25
	      cell->RemoveFromSimulation();
	    }
	  } // end cell death

            // remove RGC_mosaic_BM when mosaics are over
	  if (cell->GetInternalClock() > 2020) {
	    cell->RemoveBiologyModule(this);
	  }

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

        if (cell->GetCellType() == -1) { return; }
        int cell_type = cell->GetCellType();
        DiffusionGrid* dg = nullptr;
        // use corresponding diffusion grid
        if (cell_type == 0) {
          dg = rm->GetDiffusionGrid("sac-gcl");
        }
        else if (cell_type == 1) {
          dg = rm->GetDiffusionGrid("sac-inl");
        }

	auto& secretion_position = cell->GetPosition();
	dg->IncreaseConcentrationBy(secretion_position, 1);

        // remove Substance_secretion_BM when mosaics are over
        if (cell->GetInternalClock() > 2020) {
          cell->RemoveBiologyModule(this);
        }

      } // end if MyCell
    } // end Run()
  }; // end biologyModule Substance_secretion_BM


  // Define cell behavior for dendrites creation
  struct Internal_clock_BM : public BaseBiologyModule {
    BDM_STATELESS_BM_HEADER(Internal_clock_BM, BaseBiologyModule, 1);

  public:
    Internal_clock_BM() : BaseBiologyModule(gAllEventIds) {}

    void Run(SimObject* so) override {
      if (auto* cell = dynamic_cast<MyCell*>(so)) {
        auto* random = Simulation::GetActive()->GetRandom();

        // probability to increase internal clock
        if (random->Uniform(0, 1) < 0.96) {
          // update cell internal clock
          cell->SetInternalClock(cell->GetInternalClock() + 1);
        }

        if (cell->GetInternalClock() > 2021) {
          // remove Internal_clock_BM when not needed anymore
          cell->RemoveBiologyModule(this);
        }

      } // end if MyCell
    } // end Run()
  }; // end biologyModule Internal_clock_BM


  // Define cell behavior for dendrites creation
  struct Dendrite_creation_BM : public BaseBiologyModule {
    BDM_STATELESS_BM_HEADER(Dendrite_creation_BM, BaseBiologyModule, 1);

  public:
    Dendrite_creation_BM() : BaseBiologyModule(gAllEventIds) {}

    void Run(SimObject* so) override {
      if (auto* cell = dynamic_cast<MyCell*>(so)) {

        bool createDendrites = true;

        if (createDendrites && cell->GetInternalClock() > 2021) {
          auto* random = Simulation::GetActive()->GetRandom();

          int cell_type = cell->GetCellType();

          int dendrite_nb = 0;
          // dendrites number depending on cell type
          //NOTE: average dendrites number = 4.5; std = 1.2
          if (cell_type == 200) {
            dendrite_nb = 2 + (int)random->Uniform(1, 3);
          }
          else if (cell_type == 201) {
            dendrite_nb = 2 + (int)random->Uniform(1, 3);
          }
          else if (cell_type == 202) {
            dendrite_nb = 2 + (int)random->Uniform(1, 3);
          }
          else if (cell_type == 203) {
            dendrite_nb = 3 + (int)random->Uniform(1, 3);
          }
          else {
            dendrite_nb = 2;
          }

          for (int i = 0; i <= dendrite_nb; i++) {
            // root location - TODO: no overlap
            Double3 dendrite_root = {0,0,1};
            // create dendrites
            MyNeurite my_neurite;
            auto* ne = bdm_static_cast<MyNeurite*>(
              cell->ExtendNewNeurite(dendrite_root, &my_neurite));
            ne->AddBiologyModule(new RGC_dendrite_BM());
            ne->SetHasToRetract(false);
            ne->SetBeyondThreshold(false);
            ne->SetSubtype(cell_type);
          }

          // remove Dendrite_creation_BM when dendrites are created
          cell->RemoveBiologyModule(this);

        } // end dendrites creation
      } // end if MyCell
    } // end Run()
  }; // end biologyModule Dendrite_creation_BM


} // namespace bdm

#endif
