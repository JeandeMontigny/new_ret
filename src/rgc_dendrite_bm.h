#ifndef RGC_DENDRITE_BM_
#define RGC_DENDRITE_BM_

#include "biodynamo.h"
#include "rgc_soma_bm.h"

namespace bdm {
using namespace std;


// Define rgc dendrite behavior
struct RGC_dendrite_BM : public BaseBiologyModule {
  RGC_dendrite_BM() : BaseBiologyModule(gAllEventIds) {}
  // Default event constructor
  RGC_dendrite_BM(const Event& event, BaseBiologyModule* other,
                  uint64_t new_oid = 0)
      : BaseBiologyModule(event, other, new_oid) {}

  // Create a new instance of this object using the default constructor.
  BaseBiologyModule* GetInstance(const Event& event, BaseBiologyModule* other,
                                 uint64_t new_oid = 0) const override {
    return new RGC_dendrite_BM(event, other, new_oid);
  }

  // Create a copy of this biology module.
  BaseBiologyModule* GetCopy() const override {
    return new RGC_dendrite_BM(*this);
  }

  void Run(SimObject* so) override {
    if (auto* ne = dynamic_cast<MyNeurite*>(so)) {
      auto* sim = Simulation::GetActive();
      auto* rm = sim->GetResourceManager();
      auto* random = sim->GetRandom();

      if (!init_) {
        dg_guide_on_ =
        rm->GetDiffusionGrid("on_dendrites");
        dg_guide_off_ =
        rm->GetDiffusionGrid("off_dendrites");
        init_ = true;
      }

      if (ne->IsTerminal() && ne->GetDiameter() >= 0.5) {

        int cell_type = ne->GetMySoma()->GetCellType();
        Double3 gradient_guide;
        Double3 prefered_dir = {0, 0, 0};
        double concentration = 0;
        double on_off_factor = 2;
        bool homotypic_competition = true;
        double gradient_weight = 0.2;
        double randomness_weight = 0.4;
        double old_direction_weight = 4.5;
        double conc_retract_threshold = 0.01;
        double diam_retract_threshold = 0.85;
        double shrinkage = 0.00058;
        double bifurc_proba = 0.01*ne->GetDiameter();
        double bifurc_threshold = 0.04;

        // if on-off cells
        if (cell_type/100 == 0) {
          shrinkage = 0.0005;
          randomness_weight = 0.5;
          bifurc_proba = 0.0058 * ne->GetDiameter();
        }
        // if on cells
        if (cell_type/100 == 1) {
          shrinkage = 0.0007;
          randomness_weight = 0.4;
          bifurc_proba = 0.003 * ne->GetDiameter();
        }
        // if off cells
        if (cell_type/100 == 2) {
          shrinkage = 0.00038;
          randomness_weight = 0.5;
          bifurc_proba = 0.0044 * ne->GetDiameter();
        }

        if (ne->GetSubtype() == 0 || ne->GetSubtype() == 1
        || ne->GetSubtype() == 2 || ne->GetSubtype() == 3) {
          on_off_factor = 3;
          shrinkage = 0.00042;
          randomness_weight = 0.7;
          bifurc_proba = 0.0066 * ne->GetDiameter();
        }
        if (cell_type == 5) {
          shrinkage = 0.00074;
          randomness_weight = 0.5;
          bifurc_proba = 0.018 * ne->GetDiameter();
          bifurc_threshold = 0.038;
          diam_retract_threshold = 0.8;
        }

        if (cell_type == 100 || cell_type == 101 || cell_type == 102) {
          shrinkage = 0.00055;
          randomness_weight = 0.4;
          bifurc_proba = 0.0082 * ne->GetDiameter();
        }
        if (cell_type == 103) {
          shrinkage = 0.0005;
          randomness_weight = 0.4;
          bifurc_proba = 0.0061 * ne->GetDiameter();
        }
        if (cell_type == 104 || cell_type == 105 || cell_type == 106) {
          shrinkage = 0.00044;
          randomness_weight = 0.4;
          bifurc_proba = 0.00489 * ne->GetDiameter();
        }

        if (cell_type == 200 || cell_type == 201) {
          shrinkage = 0.00044;
          randomness_weight = 0.4;
          bifurc_proba = 0.007 * ne->GetDiameter();
        }
        if (cell_type == 202) {
          shrinkage = 0.00045;
          randomness_weight = 0.4;
          bifurc_proba = 0.0058 * ne->GetDiameter();
        }
        if (cell_type == 203) {
          shrinkage = 0.00036;
          randomness_weight = 0.6;
          bifurc_proba = 0.00682 * ne->GetDiameter();
          prefered_dir = {0.05, 0.05, 0};
        }
        if (cell_type == 204) {
          shrinkage = 0.00064;
          randomness_weight = 0.6;
          bifurc_proba = 0.011 * ne->GetDiameter();
          prefered_dir = {0.05, 0.05, 0};
        }
        if (cell_type == 205) {
          shrinkage = 0.00064;
          randomness_weight = 0.6;
          bifurc_proba = 0.011 * ne->GetDiameter();
          prefered_dir = {0.04, 0.04, 0};
        }

        // set correct concentration and gradient
        if (cell_type/100 == 0) {
          double conc_on = dg_guide_on_->GetConcentration(ne->GetPosition());
          double conc_off = dg_guide_off_->GetConcentration(ne->GetPosition());
          if (conc_on > conc_off * on_off_factor) {
            concentration = conc_on;
            dg_guide_on_->GetGradient(ne->GetPosition(), &gradient_guide);
          }
          else {
            concentration = conc_off;
            dg_guide_off_->GetGradient(ne->GetPosition(), &gradient_guide);
          }
        }
        if (cell_type/100 == 1) {
          dg_guide_on_->GetGradient(ne->GetPosition(), &gradient_guide);
          concentration = dg_guide_on_->GetConcentration(ne->GetPosition());
        }
        if (cell_type/100 == 2 || cell_type == 5) {
          dg_guide_off_->GetGradient(ne->GetPosition(), &gradient_guide);
          concentration = dg_guide_off_->GetConcentration(ne->GetPosition());
        }


        // if neurite doesn't have to retract
        if (!ne->GetHasToRetract()) {
          Double3 random_axis = {random->Uniform(-1, 1),
                                 random->Uniform(-1, 1),
                                 random->Uniform(-1, 1)};
          auto old_direction = ne->GetSpringAxis() * old_direction_weight;
          auto grad_direction = gradient_guide * gradient_weight;
          auto random_direction = random_axis * randomness_weight;
          Double3 new_step_direction =
            old_direction + grad_direction + random_direction + prefered_dir;

          ne->ElongateTerminalEnd(25, new_step_direction);
          ne->SetDiameter(ne->GetDiameter() - shrinkage);

          if (concentration > bifurc_threshold && random->Uniform() < bifurc_proba) {
            ne->SetDiameter(ne->GetDiameter() - 0.005);
            ne->Bifurcate();
          }

          // homo-type interaction
          if (homotypic_competition) {
            double squared_radius = 1.44;
            int homotypic_arbour = 0, my_arbour = 0;
            // counters for neurites neighbours
            auto count_homotypic_neighbours = [&ne, &homotypic_arbour, &my_arbour, &cell_type](const auto* neighbor) {
              if (neighbor->GetShape() == Shape::kCylinder) {
                auto* ne_neighbor = bdm_static_cast<const MyNeurite*>(neighbor);
                // if not the same soma
                if (!(ne_neighbor->GetMySoma() == ne->GetMySoma())) {
                  // if same cell type
                  if (ne_neighbor->GetMySoma()->GetCellType() == cell_type) {
                    homotypic_arbour++;
                  }
                } // end if not same soma
                // if same arbour
                else {
                  my_arbour++;
                }
              } // end if neighbor
            };

            auto* ctxt = sim->GetExecutionContext();
            ctxt->ForEachNeighborWithinRadius(count_homotypic_neighbours, *ne, squared_radius);

            // if is surrounded by homotype dendrites
            if (homotypic_arbour > my_arbour && ne->GetDiameter() < 0.85) {
              ne->SetHasToRetract(true);
              ne->SetDiamBeforeRetraction(ne->GetDiameter());
              ne->SetBeyondThreshold(false);
            }
          } // end homo-type interaction

          // if neurite is going too far away from guide
          if (concentration < conc_retract_threshold
            && ne->GetDiameter() < diam_retract_threshold) {
            ne->SetHasToRetract(true);
            ne->SetBeyondThreshold(true);
          }

        } // if ! has to retract

        // if neurite has to retract
        else {
          ne->RetractTerminalEnd(40);
          ne->SetDiameter(ne->GetDiameter()+0.00005);
          // if neurite has retracted enough because of interactions
          if (!ne->GetBeyondThreshold()
            && ne->GetDiameter() > ne->GetDiamBeforeRetraction()+0.001) {
            ne->SetHasToRetract(false);
	          // ne->RemoveBiologyModule(this);
	         }
          // if neurite is back to higher concentration
      	  if (ne->GetBeyondThreshold() && concentration > 0.02) {
      	    ne->SetBeyondThreshold(false);
      	    ne->SetHasToRetract(false);
      	  }
        }

      } // end if terminal

      // if not terminal or diameter is too small
      else {
        ne->RemoveBiologyModule(this);
      }

    } // end if MyNeurite
  } // end Run()

  private:
    bool init_ = false;
    DiffusionGrid* dg_guide_on_ = nullptr;
    DiffusionGrid* dg_guide_off_ = nullptr;
    BDM_CLASS_DEF_OVERRIDE(RGC_dendrite_BM, 1);
}; // end biologyModule RGC_dendrite_BM

}

#endif
