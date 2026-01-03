#ifndef STRESSINBOX_H_
#define STRESSINBOX_H_

#include "BaseObservable.h"
#include <fstream>
#include <string>

class StressInBox : public BaseObservable {
protected:
    // Box definition
    LR_vector _center = LR_vector(0., 0., 0.);
    double _L = 10.0;          // side length
    double _halfL = 5.0;       // computed
    double _Vbox = 1000.0;     // L^3, computed

    // Output
    bool _dump_tensor = true;
    std::string _tensor_filename = "stress_in_box_3x3.dat";
    std::ofstream _out;

    // If you want minimum-image relative to center
    inline LR_vector _pbc_displacement_to_center(const LR_vector& r) const;

    inline bool _inside_box(const LR_vector& r) const;

public:
    StressInBox();
    virtual ~StressInBox();

    void get_settings(input_file &my_inp, input_file &sim_inp) override;
    void init() override;

    // We will need particle positions/forces on CPU
    bool require_data_on_CPU() override;

    void update_data(llint curr_step) override;

    std::string get_output_string(llint curr_step) override;
};

#endif
