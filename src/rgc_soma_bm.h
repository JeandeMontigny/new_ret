#ifndef RGC_SOMA_BM_
#define RGC_SOMA_BM_

#include "biodynamo.h"
#include "extended_objects.h"

namespace bdm {

  // enumerate substances in simulation
  // enum Substances { dg_200_, dg_201_, dg_202_, dg_203_, dg_204_, dg_205_,
  //                   dg_206_, dg_207_, dg_208_, dg_209_, dg_210_, dg_211_ };
  enum Substances { dg_200_, dg_201_, dg_202_, dg_203_};

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


        /* -- cell fate -- */
        if (cell_type == -1 && cell_clock%3==0) {

          struct conc_type {
            double concentration;
            int type;
            string name;
            double probability;
            conc_type(double c, int t, string n, double p) {
              concentration = c;
              type = t;
              name = n;
              probability = p;
            }
          };
          //NOTE: density to obtain: 40, 40, 65, 200
          array<string, 4> substances_list = { "off_aplhaa", "off_aplhab", "off_m1", "off_j" };
          array<int, 4> cells_types = { 200, 201, 202, 203 };
          array<double, 4> proba = { 0.115, 0.115, 0.188, 0.579 };
          vector<conc_type> conc_type_list;
          for (size_t i=0; i < substances_list.size(); i++) {
            dg = rm->GetDiffusionGrid(substances_list[i]);
            double concentration = dg->GetConcentration(position);
              conc_type_list.push_back(conc_type(concentration, cells_types[i], substances_list[i], proba[i]) );
          }

          size_t nbOfZero = 0;
          for (size_t i = 0; i < conc_type_list.size(); i++) {
            if (conc_type_list[i].concentration==0) {
              nbOfZero++;
            }
          }

          if (nbOfZero == conc_type_list.size()) {
            if (random->Uniform(0, 1) < 0.001) {
              // TODO: non linear type attribution
              int selected_type = conc_type_list[random->Uniform(0, conc_type_list.size())].type;
              cell->SetCellType(selected_type);
            }
            return;
          }

          if (random->Uniform(0, 1) < 0.05) { return; }

          // TODO: non linear type attribution
          // lowest concentration?
          // if just lowest -> homogeneous distribution of types
          // divide concentration by wanted density? (lower concentration of denser type)


          //TODO: special case
          // if (conc_type_list.size() == 0 ) { return; }

          // int selected_type = conc_type_list[random->Uniform(0, conc_type_list.size())].type;
          // cell->SetCellType(selected_type);




          //TODO: count concentration of 0
          // if at least a 0 :
          //  select from 0 concentration (depending on conc_type_list[i].probability)
          // else :
          //  remove high concentration ( > concentration_threshold)
          //      NOTE: increase high proba cells concentration_threshold?

          // double concentration_threshold = 0.1;
          // // remove concentration from list if it is higher than threshold
          // for (size_t i = 0; i < conc_type_list.size(); i++) {
          //   if (conc_type_list[i].concentration > concentration_threshold) {
          //   if (conc_type_list[i].concentration > concentration_threshold * conc_type_list[i].probability) {
          //     conc_type_list.erase(conc_type_list.begin()+i);
          //   }
          // }


          //  select from the remaining concentrations (depending on conc_type_list[i].probability)

          // cell->SetCellType(selected_type);

          return;
        } // end cell fate


        // use corresponding diffusion grid
        if (cell_type == 200) {
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

        dg->GetGradient(position, &gradient);
        concentration = dg->GetConcentration(position);
        if (position[2]>27) {gradient_z={0, 0, -0.01};}
        else {gradient_z={0, 0, 0.01};}
        diff_gradient = gradient * -0.1; diff_gradient[2] = 0;

        bool with_movement = true;
        double movement_threshold = 1.735;
        bool with_death = true;
        double death_threshold = 1.76;

        // thresholds depending on initial density to obtain ~65% death rate
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

        /* -- cell growth -- */
        if (cell_clock < 960 && cell_clock%3==0) {
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
        if (with_movement && cell_clock >= 100 && cell_clock < 1920
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
          if (with_death && cell_clock >= 100 && cell_clock < 960
            && cell_clock%4==0) {
              // add vertical migration as the multi layer colapse in just on layer
              cell->UpdatePosition(gradient_z);
              // cell death depending on homotype substance concentration
              if (concentration > death_threshold
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
      } // end if MyCell
    } // end Run()
  }; // end biologyModule Substance_secretion_BM


// cell->RemoveBiologyModule(this);

} // namespace bdm

#endif
