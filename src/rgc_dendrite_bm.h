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
        double concentration = 0;
        double on_off_factor = 1;
        bool homotypic_competition = true;
        double gradient_weight = 0.2;
        double randomness_weight = 0.5;
        double old_direction_weight = 4.5;
        double conc_retract_threshold = 0.01;
        double shrinkage = 0.0007;
        double bifurc_proba = 0.01*ne->GetDiameter();
        double bifurc_threshold = 0.04;

        if (ne->GetSubtype() == 0 || ne->GetSubtype() == 1
        || ne->GetSubtype() == 2 || ne->GetSubtype() == 3) {
          on_off_factor = 3;
          randomness_weight = 0.8;
          shrinkage = 0.0004;
          bifurc_proba = 0.008 * ne->GetDiameter(); // 0.00675
        }

        if (ne->GetSubtype()/100 == 0) {
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
        if (cell_type/100 == 2) {
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
            old_direction + grad_direction + random_direction;

          ne->ElongateTerminalEnd(25, new_step_direction);
          ne->SetDiameter(ne->GetDiameter() - shrinkage);

          if (concentration > bifurc_threshold && random->Uniform() < bifurc_proba) {
            ne->SetDiameter(ne->GetDiameter() - 0.005);
            ne->Bifurcate();
          }

          // homo-type interaction
          if (homotypic_competition) {
            double squared_radius = 1.21;
            int sameType = 0;
            int otherType = 0;
            // counters for neurites neighbours
            rm->ApplyOnAllElements([&](SimObject* so, SoHandle) {
              auto* neighbor = dynamic_cast<MyNeurite*>(so);
              if (neighbor) {
                Double3 neighbor_position = neighbor->GetPosition();
                Double3 ne_position = ne->GetPosition();
                double distance =
                  pow(ne_position[0] - neighbor_position[0], 2) +
                  pow(ne_position[1] - neighbor_position[1], 2) +
                  pow(ne_position[2] - neighbor_position[2], 2);
                // if not the same soma
                if (distance < squared_radius) {
                  if (!(neighbor->GetMySoma() == ne->GetMySoma())) {
                    if (neighbor->GetMySoma()->GetCellType() == cell_type) {
                      sameType++;
                    }
                    else {
                      otherType++;
                    }
                  }
                } // if within radius
              }
            });
            // if is surrounded by homotype dendrites
            if (sameType > otherType) {
              ne->SetHasToRetract(true);
              ne->SetDiamBeforeRetraction(ne->GetDiameter());
            }
          } // end homo-type interaction

          // if neurite is going too far away from guide
          if (concentration < conc_retract_threshold && ne->GetDiameter() < 0.85) {
            ne->SetHasToRetract(true);
            ne->SetBeyondThreshold(true);
          }

        } // if ! has to retract

        // if neurite has to retract
        else {
          ne->RetractTerminalEnd(40);
          // if neurite has retracted enough because of interactions
        	if (!ne->GetBeyondThreshold()
            && ne->GetDiameter() > ne->GetDiamBeforeRetraction()+0.02) {
            ne->SetHasToRetract(false);
        	  ne->RemoveBiologyModule(this);
        	}
          // if neurite is back to higher concentration
        	if (ne->GetBeyondThreshold() && concentration>0.02) {
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
