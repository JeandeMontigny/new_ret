#ifndef RGC_SOMA_BM_
#define RGC_SOMA_BM_

#include "biodynamo.h"
#include "extended_objects.h"
#include "rgc_dendrite_bm.h"
#include "util_methods.h"

namespace bdm {

  // enumerate substances in simulation
  // enum Substances { dg_200_, dg_201_, dg_202_, dg_203_, dg_204_, dg_205_,
  //                   dg_206_, dg_207_, dg_208_, dg_209_, dg_210_, dg_211_ };
  // enum Substances { dg_200_, dg_201_, dg_202_, dg_203_}
  enum Substances { on_dendrites, off_dendrites};

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
        Double3 gradient, diff_gradient, gradient_z;
        DiffusionGrid* dg = nullptr;

        bool with_movement = true;
        double movement_threshold = 1.735;
        bool with_death = true;
        double death_threshold = 1.76;

        if (false) {
          // 1000 (350 final) -- RI ~6-7
          movement_threshold = 1.745;
          death_threshold = 1.765;
          // 800 (280 final) -- RI ~5-7
          movement_threshold = 1.735;
          death_threshold = 1.76;
          // 600 (210 final) -- RI ~5-7
          movement_threshold = 1.73;
          death_threshold = 1.77;
          // 400 (140 final) -- RI ~4-5
          movement_threshold = 1.72;
          death_threshold = 1.775;
          // 200 (70 final) -- RI ~2.5-3
          movement_threshold = 1.71;
          death_threshold = 1.78;
          // 100 (35 final) -- RI ~2-3
          movement_threshold = 1.7;
          death_threshold = 1.79;
          // 60 (20 final) -- RI ~1.8-2.5
          movement_threshold = 1.7;
          death_threshold = 1.78;
        }

        /* -- cell fate -- */
        if (cell_type == -1) {
          if (cell_clock%2!=0 || random->Uniform(0, 1) < 0.94) { return; }

          struct conc_type {
            double concentration;
            int type;
            double probability;
            conc_type(double c, int t, double p) {
              concentration = c;
              type = t;
              probability = p;
            }
          };

          // density to obtain: 114, 114, 185, 571
          array<string, 4> substances_list = { "off_aplhaa", "off_aplhab", "off_m1", "off_j" };
          array<int, 4> cells_types = { 200, 201, 202, 203 };
          array<double, 4> proba = { 0.115, 0.115, 0.188, 0.58 };
          vector<conc_type> conc_type_list;
          for (size_t i=0; i < substances_list.size(); i++) {
            dg = rm->GetDiffusionGrid(substances_list[i]);
            double concentration = dg->GetConcentration(position);
            conc_type_list.push_back(conc_type(concentration, cells_types[i], proba[i]) );
          }

          double concentration_threshold = 1e-3;
          vector<conc_type> conc_type_list_potential;
          size_t nb_zero; double sum_proba;
          do {
            nb_zero = 0; sum_proba = 0;
            conc_type_list_potential.clear();
            for (size_t i = 0; i < conc_type_list.size(); i++) {
              if (conc_type_list[i].concentration < concentration_threshold * pow(conc_type_list[i].probability, 3)) {
              // if (conc_type_list[i].concentration < concentration_threshold) {
                conc_type_list_potential.push_back(conc_type_list[i]);
                sum_proba += conc_type_list[i].probability;
                if (conc_type_list[i].concentration == 0) {
                  nb_zero++;
                }
              }
            }
            concentration_threshold *= 10;
          } while (conc_type_list_potential.size() == 0);

          // if no substances around
          if (nb_zero == conc_type_list.size() && random->Uniform(0, 1) < 0.9) {
            return;
          }

          vector<double> cumulative_proba; double previous_proba = 0;
          for (size_t i = 0; i < conc_type_list_potential.size(); i++) {
            cumulative_proba.push_back(conc_type_list_potential[i].probability+previous_proba);
            previous_proba = cumulative_proba[i];
          }

          double random_double = random->Uniform(0, sum_proba);
          size_t j = 0;
          while (random_double > cumulative_proba[j]) {
            j++;
          }

          cell->SetCellType(conc_type_list_potential[j].type);
          return;
        } // end cell fate

        /* -- initialisation -- */
        // use corresponding diffusion grid
        // set thresholds depending on initial density to obtain ~65% death rate
        if (cell_type == 200) {
          dg = rm->GetDiffusionGrid("off_aplhaa");
          movement_threshold = 1.7;
          death_threshold = 1.79;
        }
        else if (cell_type == 201) {
          dg = rm->GetDiffusionGrid("off_aplhab");
          movement_threshold = 1.7;
          death_threshold = 1.79;
        }
        else if (cell_type == 202) {
          dg = rm->GetDiffusionGrid("off_m1");
          movement_threshold = 1.71;
          death_threshold = 1.78;
        }
        else if (cell_type == 203) {
          dg = rm->GetDiffusionGrid("off_j");
          movement_threshold = 1.727;
          death_threshold = 1.772;
        }

        dg->GetGradient(position, &gradient);
        concentration = dg->GetConcentration(position);
        if (position[2]>27) {gradient_z={0, 0, -0.01};}
        else {gradient_z={0, 0, 0.01};}
        diff_gradient = gradient * -0.1; diff_gradient[2] = 0;

        /* -- cell growth -- */
        if (cell_clock >= 100 && cell_clock < 1060 && cell_clock%3==0) {
          // // add small random movements
          cell->UpdatePosition(
              {random->Uniform(-0.01, 0.01), random->Uniform(-0.01, 0.01), 0});
          // cell growth
          if (cell->GetDiameter() < 14 && random->Uniform(0, 1) < 0.02) {
            cell->ChangeVolume(3000);
          }
          // layer colapse if no cell death
          if (!with_death) {
            cell->UpdatePosition(gradient_z);
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
          if (with_death && cell_clock >= 200 && cell_clock < 1060
            && cell_clock%4==0) {
              // add vertical migration as the multi layer colapse in just on layer
              cell->UpdatePosition(gradient_z);
              // cell death depending on homotype substance concentration
              if (concentration > death_threshold
                  && random->Uniform(0, 1) < 0.1) { // 0.25
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
        if (cell_type == 200) {
          dg = rm->GetDiffusionGrid("off_aplhaa");
        }
        if (cell_type == 201) {
          dg = rm->GetDiffusionGrid("off_aplhab");
        }
        if (cell_type == 202) {
          dg = rm->GetDiffusionGrid("off_m1");
        }
        if (cell_type == 203) {
          dg = rm->GetDiffusionGrid("off_j");
        }

        if (cell->GetInternalClock()%3==0) {
          auto& secretion_position = cell->GetPosition();
          dg->IncreaseConcentrationBy(secretion_position, 1);
        }

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
          // average dendrites number = 4.5; std = 1.2
          int dendrite_nb = (int)random->Uniform(2, 7);
          // dendrites number depending on cell type
          if (cell_type == 0 || cell_type == 1
            || cell_type == 2 || cell_type == 3) {
            // dendrite_nb = RandomPoisson(4.7);
            double L = exp(-4.7);
            int k = 0;
            double p = 1;
            do {
              k = k + 1;
              double u = random->Uniform(0.15, 1);
              p = p * u;
            } while (p > L);
            dendrite_nb = k - 1;
            if (dendrite_nb < 3) { dendrite_nb = 3; }
            else if (dendrite_nb > 6) { dendrite_nb = 3; }
          }
	  if (cell_type/100 == 1) {
            dendrite_nb = (int)random->Uniform(3.5, 7.5);
          }
	  if (cell_type/100 == 2) {
            dendrite_nb = (int)random->Uniform(2.8, 6.5);
          }
          if (cell_type == 5) {
            dendrite_nb = (int)random->Uniform(3, 6);
          }
          if (cell_type == 203) {
            dendrite_nb = (int)random->Uniform(2.8, 5.8);
          }

          // if (cell_type == 200) {
          //   dendrite_nb = 2 + (int)random->Uniform(1, 3);
          // }
          // else if (cell_type == 201) {
          //   dendrite_nb = 2 + (int)random->Uniform(1, 3);
          // }
          // else if (cell_type == 202) {
          //   dendrite_nb = 2 + (int)random->Uniform(1, 3);
          // }
          // else if (cell_type == 203) {
          //   dendrite_nb = 3 + (int)random->Uniform(1, 3);
          // }

          for (int i = 0; i < dendrite_nb; i++) {
            // root location - TODO: no overlap
            Double3 dendrite_root = {0,0,1};
            // create dendrites
            MyNeurite my_neurite;
            auto* ne = bdm_static_cast<MyNeurite*>(
              cell->ExtendNewNeurite(dendrite_root, &my_neurite));
            ne->SetDiameter(1);
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
