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

        int cell_type = ne->GetMySoma()->GetLabel();
        Double3 gradient_guide;
        double concentration = 0;

        double gradient_weight = 0.2;
        double randomness_weight = 0.5;
        double old_direction_weight = 4.5;
        double concentration_threshold = 0.01;

        if (ne->GetSubtype()/100 == 0) {
          double conc_on = dg_guide_on_->GetConcentration(ne->GetPosition());
          double conc_off = dg_guide_off_->GetConcentration(ne->GetPosition());
          if (conc_on > conc_off) {
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


        // ---- NOTE: FOR TEST ONLY ---- //
        if (cell_type == 202) {
          dg_guide_on_->GetGradient(ne->GetPosition(), &gradient_guide);
          concentration = dg_guide_on_->GetConcentration(ne->GetPosition());
        }
        if (cell_type == 201 || cell_type == 200) {
          double conc_on = dg_guide_on_->GetConcentration(ne->GetPosition());
          double conc_off = dg_guide_off_->GetConcentration(ne->GetPosition());
          if (conc_on > conc_off) {
            concentration = conc_on;
            dg_guide_on_->GetGradient(ne->GetPosition(), &gradient_guide);
          }
          else {
            concentration = conc_off;
            dg_guide_off_->GetGradient(ne->GetPosition(), &gradient_guide);
          }
        }
        // ---- END FOR TEST ONLY ---- //

        // if neurite doesn't have to retract
        if (!ne->GetHasToRetract()) {
          double bifurcProba = 0.01*ne->GetDiameter();

          Double3 random_axis = {random->Uniform(-1, 1),
                                 random->Uniform(-1, 1),
                                 random->Uniform(-1, 1)};
          auto old_direction = ne->GetSpringAxis() * old_direction_weight;
          auto grad_direction = gradient_weight * gradient_weight;
          auto random_direction = random_axis * randomness_weight;
          Double3 new_step_direction =
            old_direction + grad_direction + random_direction;

          ne->ElongateTerminalEnd(25, new_step_direction);
          ne->SetDiameter(ne->GetDiameter()-0.0007);

          if (concentration > 0.04 && random->Uniform() < bifurcProba) {
            ne->SetDiameter(ne->GetDiameter()-0.005);
            ne->Bifurcate();
          }

          // homo-type interaction
          double squared_radius = 1.1; // 1.2
          int sameType = 0;
          int otherType = 0;
          // counters for neurites neighbours
          rm->ApplyOnAllElements([&](SimObject* so, SoHandle) {
            auto* neighbor = dynamic_cast<MyNeurite*>(so);
            if (neighbor) {
              // if not the same soma
              if (!(neighbor->GetMySoma() == ne->GetMySoma())) {
                if (neighbor->GetMySoma()->GetCellType() == cell_type) {
                  sameType++;
                }
                else {
                  otherType++;
                }
              }
            }
          });

          // if is surrounded by homotype dendrites
          if (sameType >= otherType) {
            ne->SetHasToRetract(true);
            ne->SetDiamBeforeRetraction(ne->GetDiameter());
          }

          // if neurite is going too far away from guide
          if (concentration < 0.01 && ne->GetDiameter() < 0.9) {
  					ne->SetHasToRetract(true);
  					ne->SetBeyondThreshold(true);
  				}

          // if neurite is going too far away from guide
          if (concentration < concentration_threshold && ne->GetDiameter() < 0.9) {
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
