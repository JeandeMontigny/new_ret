#ifndef EXTENDED_OBJECTS_
#define EXTENDED_OBJECTS_

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
     double distance_travelled_;
  };

}

#endif
