#ifndef RGC_SOMA_BM_
#define RGC_SOMA_BM_

#include "biodynamo.h"
#include "extended_objects.h"
#include "rgc_dendrite_bm.h"

namespace bdm {

  // enumerate substances in simulation
  enum Substances { dg_000_, dg_001_, dg_002_, dg_003_, dg_004_, dg_005_,
                    dg_006_, dg_007_, dg_008_, dg_009_, dg_010_, dg_011_,
                    dg_100_, dg_101_, dg_102_, dg_103_, dg_104_, dg_105_,
                    dg_106_, dg_107_, dg_108_, dg_109_, dg_110_, dg_111_,
                    dg_112_, dg_113_, dg_114_, dg_115_, dg_116_, dg_117_,
                    dg_118_,
                    dg_200_, dg_201_, dg_202_, dg_203_, dg_204_, dg_205_,
                    dg_206_, dg_207_, dg_208_, dg_209_, dg_210_, dg_211_ };

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

        bool with_movement = false;
        double movement_threshold = 1.735;
        bool with_death = true;
        double death_threshold = 1.76;

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
          array<string, 43> substances_list =
            { "on-off_dsgca", "on-off_dsgcb", "on-off_dsgcc", "on-off_dsgcd",
              "on-off_m3", "on-off_led", "on-off_u", "on-off_v", "on-off_w",
              "on-off_x", "on-off_y", "on-off_z",
              "on_dsgca", "on_dsgcb", "on_dsgcc", "on_aplha", "on_m2", "on_m4",
              "on_m5", "on_o", "on_p","on_q", "on_r", "on_s", "on_t", "on_u",
              "on_v", "on_w", "on_x", "on_y", "on_z",
              "off_aplhaa", "off_aplhab", "off_m1", "off_j", "off_mini_j",
              "off_midi_j", "off_u", "off_v", "off_w", "off_x", "off_y", "off_z" };
          array<int, 43> cells_types =
            { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
              100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
              112, 113, 114, 115, 116, 117, 118,
              200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211 };
          array<double, 43> proba =
            { 0.04166, 0.04166, 0.04166, 0.04166, 0.00666, 0.08333, 0.00666,
              0.00666, 0.02, 0.01666, 0.01333, 0.01333,
              0.01333, 0.01333, 0.01333, 0.01333, 0.01866, 0.00666, 0.00666, 0.05,
              0.03333, 0.03333, 0.02666, 0.02, 0.02, 0.01666, 0.01666, 0.01133,
              0.00666, 0.00666, 0.00666,
              0.01333, 0.01333, 0.021, 0.06666, 0.11666, 0.02666, 0.00666, 0.00666,
              0.02, 0.01666, 0.01333, 0.01233 };
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
        // on-off
        if (cell_type == 0) {
          dg = rm->GetDiffusionGrid("on-off_dsgca");
          movement_threshold = 1.717;
          death_threshold = 1.776;
        }
        else if (cell_type == 1) {
          dg = rm->GetDiffusionGrid("on-off_dsgcb");
          movement_threshold = 1.717;
          death_threshold = 1.776;
        }
        else if (cell_type == 2) {
          dg = rm->GetDiffusionGrid("on-off_dsgcc");
          movement_threshold = 1.717;
          death_threshold = 1.776;
        }
        else if (cell_type == 3) {
          dg = rm->GetDiffusionGrid("on-off_dsgcd");
          movement_threshold = 1.717;
          death_threshold = 1.776;
        }
        else if (cell_type == 4) {
          dg = rm->GetDiffusionGrid("on-off_m3");
          movement_threshold = 1.7;
          death_threshold = 1.78;
        }
        else if (cell_type == 5) {
          dg = rm->GetDiffusionGrid("on-off_led");
          movement_threshold = 1.735;
          death_threshold = 2.05; // 1.77
        }
        else if (cell_type == 6) {
          dg = rm->GetDiffusionGrid("on-off_u");
          movement_threshold = 1.7;
          death_threshold = 1.78;
        }
        else if (cell_type == 7) {
          dg = rm->GetDiffusionGrid("on-off_v");
          movement_threshold = 1.7;
          death_threshold = 1.78;
        }
        else if (cell_type == 8) {
          dg = rm->GetDiffusionGrid("on-off_w");
          movement_threshold = 1.708;
          death_threshold = 1.782;
        }
        else if (cell_type == 9) {
          dg = rm->GetDiffusionGrid("on-off_x");
          movement_threshold = 1.705;
          death_threshold = 1.785;
        }
        else if (cell_type == 10) {
          dg = rm->GetDiffusionGrid("on-off_y");
          movement_threshold = 1.701;
          death_threshold = 1.789;
        }
        else if (cell_type == 11) {
          dg = rm->GetDiffusionGrid("on-off_z");
          movement_threshold = 1.701;
          death_threshold = 1.789;
        }
        // on
        else if (cell_type == 100) {
          dg = rm->GetDiffusionGrid("on_dsgca");
          movement_threshold = 1.701;
          death_threshold = 1.789;
        }
        else if (cell_type == 101) {
          dg = rm->GetDiffusionGrid("on_dsgcb");
          movement_threshold = 1.701;
          death_threshold = 1.789;
        }
        else if (cell_type == 102) {
          dg = rm->GetDiffusionGrid("on_dsgcc");
          movement_threshold = 1.701;
          death_threshold = 1.789;
        }
        else if (cell_type == 103) {
          dg = rm->GetDiffusionGrid("on_aplha");
          movement_threshold = 1.701;
          death_threshold = 1.789;
        }
        else if (cell_type == 104) {
          dg = rm->GetDiffusionGrid("on_m2");
          movement_threshold = 1.705;
          death_threshold = 1.785;
        }
        else if (cell_type == 105) {
          dg = rm->GetDiffusionGrid("on_m4");
          movement_threshold = 1.7;
          death_threshold = 1.78;
        }
        else if (cell_type == 106) {
          dg = rm->GetDiffusionGrid("on_m5");
          movement_threshold = 1.7;
          death_threshold = 1.78;
        }
        else if (cell_type == 107) {
          dg = rm->GetDiffusionGrid("on_o");
          movement_threshold = 1.72;
          death_threshold = 1.774;
        }
        else if (cell_type == 108) {
          dg = rm->GetDiffusionGrid("on_p");
          movement_threshold = 1.714;
          death_threshold = 1.778;
        }
        else if (cell_type == 109) {
          dg = rm->GetDiffusionGrid("on_q");
          movement_threshold = 1.714;
          death_threshold = 1.778;
        }
        else if (cell_type == 110) {
          dg = rm->GetDiffusionGrid("on_r");
          movement_threshold = 1.711;
          death_threshold = 1.779;
        }
        else if (cell_type == 111) {
          dg = rm->GetDiffusionGrid("on_s");
          movement_threshold = 1.708;
          death_threshold = 1.782;
        }
        else if (cell_type == 112) {
          dg = rm->GetDiffusionGrid("on_t");
          movement_threshold = 1.708;
          death_threshold = 1.782;
        }
        else if (cell_type == 113) {
          dg = rm->GetDiffusionGrid("on_u");
          movement_threshold = 1.705;
          death_threshold = 1.785;
        }
        else if (cell_type == 114) {
          dg = rm->GetDiffusionGrid("on_v");
          movement_threshold = 1.705;
          death_threshold = 1.785;
        }
        else if (cell_type == 115) {
          dg = rm->GetDiffusionGrid("on_w");
          movement_threshold = 1.7;
          death_threshold = 1.79;
        }
        else if (cell_type == 116) {
          dg = rm->GetDiffusionGrid("on_x");
          movement_threshold = 1.7;
          death_threshold = 1.78;
        }
        else if (cell_type == 117) {
          dg = rm->GetDiffusionGrid("on_y");
          movement_threshold = 1.7;
          death_threshold = 1.78;
        }
        else if (cell_type == 118) {
          dg = rm->GetDiffusionGrid("on_z");
          movement_threshold = 1.7;
          death_threshold = 1.78;
        }
        // off
        else if (cell_type == 200) {
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
          movement_threshold = 1.709;
          death_threshold = 1.78;
        }
        else if (cell_type == 203) {
          dg = rm->GetDiffusionGrid("off_j");
          movement_threshold = 1.726;
          death_threshold = 1.772;
        }
        else if (cell_type == 204) {
          dg = rm->GetDiffusionGrid("off_mini_j");
          movement_threshold = 1.744;
          death_threshold = 1.765;
        }
        else if (cell_type == 205) {
          dg = rm->GetDiffusionGrid("off_midi_j");
          movement_threshold = 1.71;
          death_threshold = 1.779;
        }
        else if (cell_type == 206) {
          dg = rm->GetDiffusionGrid("off_u");
          movement_threshold = 1.7;
          death_threshold = 1.78;
        }
        else if (cell_type == 207) {
          dg = rm->GetDiffusionGrid("off_v");
          movement_threshold = 1.7;
          death_threshold = 1.78;
        }
        else if (cell_type == 208) {
          dg = rm->GetDiffusionGrid("off_w");
          movement_threshold = 1.708;
          death_threshold = 1.782;
        }
        else if (cell_type == 209) {
          dg = rm->GetDiffusionGrid("off_x");
          movement_threshold = 1.705;
          death_threshold = 1.785;
        }
        else if (cell_type == 210) {
          dg = rm->GetDiffusionGrid("off_y");
          movement_threshold = 1.701;
          death_threshold = 1.789;
        }
        else if (cell_type == 211) {
          dg = rm->GetDiffusionGrid("off_z");
          movement_threshold = 1.701;
          death_threshold = 1.789;
        }
        else {
          cout << "error: no valid cell type" << endl;
          return;
        }

        dg->GetGradient(position, &gradient);
        concentration = dg->GetConcentration(position);
        if (position[2]>27) {gradient_z={0, 0, -0.02};}
        else {gradient_z={0, 0, 0.02};}
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
        if (cell_type == 0) {
          dg = rm->GetDiffusionGrid("on-off_dsgca");
        }
        else if (cell_type == 1) {
          dg = rm->GetDiffusionGrid("on-off_dsgcb");
        }
        else if (cell_type == 2) {
          dg = rm->GetDiffusionGrid("on-off_dsgcc");
        }
        else if (cell_type == 3) {
          dg = rm->GetDiffusionGrid("on-off_dsgcd");
        }
        else if (cell_type == 4) {
          dg = rm->GetDiffusionGrid("on-off_m3");
        }
        else if (cell_type == 5) {
          dg = rm->GetDiffusionGrid("on-off_led");
        }
        else if (cell_type == 6) {
          dg = rm->GetDiffusionGrid("on-off_u");
        }
        else if (cell_type == 7) {
          dg = rm->GetDiffusionGrid("on-off_v");
        }
        else if (cell_type == 8) {
          dg = rm->GetDiffusionGrid("on-off_w");
        }
        else if (cell_type == 9) {
          dg = rm->GetDiffusionGrid("on-off_x");
        }
        else if (cell_type == 10) {
          dg = rm->GetDiffusionGrid("on-off_y");
        }
        else if (cell_type == 11) {
          dg = rm->GetDiffusionGrid("on-off_z");
        }
        // on
        else if (cell_type == 100) {
          dg = rm->GetDiffusionGrid("on_dsgca");
        }
        else if (cell_type == 101) {
          dg = rm->GetDiffusionGrid("on_dsgcb");
        }
        else if (cell_type == 102) {
          dg = rm->GetDiffusionGrid("on_dsgcc");
        }
        else if (cell_type == 103) {
          dg = rm->GetDiffusionGrid("on_aplha");
        }
        else if (cell_type == 104) {
          dg = rm->GetDiffusionGrid("on_m2");
        }
        else if (cell_type == 105) {
          dg = rm->GetDiffusionGrid("on_m4");
        }
        else if (cell_type == 106) {
          dg = rm->GetDiffusionGrid("on_m5");
        }
        else if (cell_type == 107) {
          dg = rm->GetDiffusionGrid("on_o");
        }
        else if (cell_type == 108) {
          dg = rm->GetDiffusionGrid("on_p");
        }
        else if (cell_type == 109) {
          dg = rm->GetDiffusionGrid("on_q");
        }
        else if (cell_type == 110) {
          dg = rm->GetDiffusionGrid("on_r");
        }
        else if (cell_type == 111) {
          dg = rm->GetDiffusionGrid("on_s");
        }
        else if (cell_type == 112) {
          dg = rm->GetDiffusionGrid("on_t");
        }
        else if (cell_type == 113) {
          dg = rm->GetDiffusionGrid("on_u");
        }
        else if (cell_type == 114) {
          dg = rm->GetDiffusionGrid("on_v");
        }
        else if (cell_type == 115) {
          dg = rm->GetDiffusionGrid("on_w");
        }
        else if (cell_type == 116) {
          dg = rm->GetDiffusionGrid("on_x");
        }
        else if (cell_type == 117) {
          dg = rm->GetDiffusionGrid("on_y");
        }
        else if (cell_type == 118) {
          dg = rm->GetDiffusionGrid("on_z");
        }
        // off
        else if (cell_type == 200) {
          dg = rm->GetDiffusionGrid("off_aplhaa");
        }
        else if (cell_type == 201) {
          dg = rm->GetDiffusionGrid("off_aplhab");
        }
        else if (cell_type == 202) {
          dg = rm->GetDiffusionGrid("off_m1");
        }
        else if (cell_type == 203) {
          dg = rm->GetDiffusionGrid("off_j");
        }
        else if (cell_type == 204) {
          dg = rm->GetDiffusionGrid("off_mini_j");
        }
        else if (cell_type == 205) {
          dg = rm->GetDiffusionGrid("off_midi_j");
        }
        else if (cell_type == 206) {
          dg = rm->GetDiffusionGrid("off_u");
        }
        else if (cell_type == 207) {
          dg = rm->GetDiffusionGrid("off_v");
        }
        else if (cell_type == 208) {
          dg = rm->GetDiffusionGrid("off_w");
        }
        else if (cell_type == 209) {
          dg = rm->GetDiffusionGrid("off_x");
        }
        else if (cell_type == 210) {
          dg = rm->GetDiffusionGrid("off_y");
        }
        else if (cell_type == 211) {
          dg = rm->GetDiffusionGrid("off_z");
        }
        else {
          cout << "error: no valid cell type" << endl;
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
