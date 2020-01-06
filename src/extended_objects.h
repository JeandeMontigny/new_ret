#ifndef EXTENDED_OBJECTS_
#define EXTENDED_OBJECTS_

#include "core/sim_object/sim_object.h"
#include "neuroscience/neuroscience.h"

namespace bdm {

  // Define custom cell MyCell extending NeuronSoma
  class MyCell : public experimental::neuroscience::NeuronSoma {
    BDM_SIM_OBJECT_HEADER(MyCell, experimental::neuroscience::NeuronSoma, 1,
                          cell_type_, internal_clock_, swc_label_,
                          previous_position_, distance_travelled_);

   public:
    MyCell() : Base() {}

    virtual ~MyCell() {}

    MyCell(const Double3& position) : Base(position) {}

    /// Default event constructor
    MyCell(const Event& event, SimObject* other, uint64_t new_oid = 0)
        : Base(event, other, new_oid) {}

    void SetCellType(int t) { cell_type_ = t; }
    int GetCellType() const { return cell_type_; }

    void SetInternalClock(int t) { internal_clock_ = t; }
    int GetInternalClock() const { return internal_clock_; }

    inline void SetLabel(int label) { swc_label_ = label; }
    inline int GetLabel() const { return swc_label_; }
    inline void IncreaseLabel() { swc_label_ += 1; }

    void SetPreviousPosition(Double3 position) { previous_position_ = position; }
    const Double3& GetPreviousPosition() const { return previous_position_; }

    void SetDistanceTravelled(double distance) { distance_travelled_ = distance; }
    double GetDistanceTravelled() const {return distance_travelled_; }

   private:
     int cell_type_ = -1;
     int internal_clock_ = 0;
     int swc_label_ = 0;
     Double3 previous_position_;
     double distance_travelled_ = 0;
  }; // end MyCell definition


  // Define custom neurite MyNeurite extending NeuriteElement
  class MyNeurite : public experimental::neuroscience::NeuriteElement {
    BDM_SIM_OBJECT_HEADER(MyNeurite, experimental::neuroscience::NeuriteElement, 1,
                          has_to_retract_, beyond_threshold_,
                          diam_before_retract_, subtype_, its_soma_);

   public:
    MyNeurite() : Base() {}

    virtual ~MyNeurite() {}

    // Default event constructor
    MyNeurite(const Event& event, SimObject* other, uint64_t new_oid = 0)
        : Base(event, other, new_oid) {
      if (event.GetId() ==
      experimental::neuroscience::NewNeuriteExtensionEvent::kEventId) {
        its_soma_ = static_cast<MyCell*>(other)->GetSoPtr<MyCell>();
      } else {
        its_soma_ = static_cast<MyNeurite*>(other)->its_soma_;
      }
    }

    // Default event handler
    void EventHandler(const Event& event, SimObject* other1,
                      SimObject* other2 = nullptr) {
      Base::EventHandler(event, other1, other2);
    }

    void SetHasToRetract(int r) { has_to_retract_ = r; }
    bool GetHasToRetract() const { return has_to_retract_; }

    void SetBeyondThreshold(int r) { beyond_threshold_ = r; }
    bool GetBeyondThreshold() const { return beyond_threshold_; }

    void SetDiamBeforeRetraction(double d) { diam_before_retract_ = d; }
    double GetDiamBeforeRetraction() const { return diam_before_retract_; }

    void SetSubtype(int st) { subtype_ = st; }
    int GetSubtype() { return subtype_; }

    void SetMySoma(SoPointer<MyCell> soma) { its_soma_ = soma; }
    SoPointer<MyCell> GetMySoma() { return its_soma_; }
    SoPointer<MyCell> GetMySoma() const { return its_soma_; }

   private:
     bool has_to_retract_;
     bool beyond_threshold_;
     double diam_before_retract_;
     int subtype_;
     SoPointer<MyCell> its_soma_;
  }; // end MyNeurite definition


} // namespace bdm

#endif
