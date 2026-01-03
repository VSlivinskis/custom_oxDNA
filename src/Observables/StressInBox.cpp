#include "StressInBox.h"

#include <sstream>
#include <cmath>
#include <fstream>   // for header-if-empty check

StressInBox::StressInBox() : BaseObservable() {}

StressInBox::~StressInBox() {
    if(_out.is_open()) _out.close();
}

void StressInBox::get_settings(input_file &my_inp, input_file &sim_inp) {
    BaseObservable::get_settings(my_inp, sim_inp);

    // side length
    getInputDouble(&my_inp, "L", &_L, 0);

    // center (default 0,0,0)
    double cx = 0., cy = 0., cz = 0.;
    getInputDouble(&my_inp, "cx", &cx, 0);
    getInputDouble(&my_inp, "cy", &cy, 0);
    getInputDouble(&my_inp, "cz", &cz, 0);
    _center = LR_vector(cx, cy, cz);

    // output control
    getInputBool(&my_inp, "dump_tensor", &_dump_tensor, 0);
    getInputString(&my_inp, "tensor_file", _tensor_filename, 0);

    _halfL = 0.5 * _L;
    _Vbox  = _L * _L * _L;

    OX_LOG(Logger::LOG_INFO,
           "Initialising StressInBox with center=(%.6f,%.6f,%.6f), L=%.6f, tensor_file=%s",
           _center.x, _center.y, _center.z, _L, _tensor_filename.c_str());
}

void StressInBox::init() {
    if(_dump_tensor) {
        // Check if file exists + empty (to write header once)
        bool is_empty = true;
        {
            std::ifstream fin(_tensor_filename.c_str());
            if(fin.good()) {
                is_empty = (fin.peek() == std::ifstream::traits_type::eof());
            }
        }

        _out.open(_tensor_filename.c_str(), std::ios::out | std::ios::app);
        if(!_out.good()) {
            throw oxDNAException("StressInBox: cannot open '%s' for writing", _tensor_filename.c_str());
        }

        if(is_empty) {
            _out << "# step  "
                 << "xx xy xz  "
                 << "yx yy yz  "
                 << "zx zy zz  "
                 << "N_in_box"
                 << std::endl;
            _out.flush();
        }
    }
}

bool StressInBox::require_data_on_CPU() {
    // Keep true (BaseObservable default is true too), but explicit makes intent clear.
    return true;
}

// In this oxDNA version BaseBox::min_image takes TWO args and returns the minimum-image displacement
LR_vector StressInBox::_pbc_displacement_to_center(const LR_vector& r) const {
    return _config_info->box->min_image(r, _center);
}

bool StressInBox::_inside_box(const LR_vector& r) const {
    const LR_vector dr = _pbc_displacement_to_center(r);
    return (std::fabs(dr.x) <= _halfL &&
            std::fabs(dr.y) <= _halfL &&
            std::fabs(dr.z) <= _halfL);
}

void StressInBox::update_data(llint curr_step) {
    // IMPORTANT: in your codebase, particles is a METHOD
    auto &particles = _config_info->particles();

    // Tensor accumulator
    double xx = 0., xy = 0., xz = 0.;
    double yx = 0., yy = 0., yz = 0.;
    double zx = 0., zy = 0., zz = 0.;
    llint n_in = 0;

    for(size_t i = 0; i < particles.size(); i++) {
        BaseParticle *p = particles[i];

        const LR_vector r = p->pos;
        if(!_inside_box(r)) continue;

        n_in++;

        // total force on particle
        const LR_vector F  = p->force;

        // displacement from center with MIC
        const LR_vector dr = _pbc_displacement_to_center(r);

        // virial proxy: dr ⊗ F
        xx += dr.x * F.x;  xy += dr.x * F.y;  xz += dr.x * F.z;
        yx += dr.y * F.x;  yy += dr.y * F.y;  yz += dr.y * F.z;
        zx += dr.z * F.x;  zy += dr.z * F.y;  zz += dr.z * F.z;
    }

    // Convert to "stress" by dividing by Vbox
    const double invV = 1.0 / _Vbox;

    xx *= invV; xy *= invV; xz *= invV;
    yx *= invV; yy *= invV; yz *= invV;
    zx *= invV; zy *= invV; zz *= invV;

    if(_dump_tensor && _out.is_open()) {
        _out << curr_step << "  "
             << xx << " " << xy << " " << xz << "  "
             << yx << " " << yy << " " << yz << "  "
             << zx << " " << zy << " " << zz << "  "
             << n_in
             << "\n";
        _out.flush();
    }

    _times_updated++;
}

std::string StressInBox::get_output_string(llint curr_step) {
    // This observable is file-writing. Returning empty avoids clutter in any generic stream.
    (void) curr_step;
    return std::string();
}
