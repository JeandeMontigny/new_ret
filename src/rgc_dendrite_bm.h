#ifndef RGC_DENDRITE_BM_
#define RGC_DENDRITE_BM_

#include "biodynamo.h"

namespace bdm {
using namespace std;


// Define rgc dendrite behavior
struct RGC_dendrite_BM : public BaseBiologyModule {
  BDM_STATELESS_BM_HEADER(RGC_dendrite_BM, BaseBiologyModule, 1);

public:
  RGC_dendrite_BM() : BaseBiologyModule(gAllEventIds) {}

  void Run(SimObject* so) override {
    if (auto* ne = dynamic_cast<MyNeurite*>(so)) {
      auto* sim = Simulation::GetActive();
      auto* random = sim->GetRandom();

      ne->GetMySoma()->GetLabel();

    } // end if MyNeurite
  } // end Run()
}; // end biologyModule RGC_dendrite_BM

}

#endif
