#include "StressInBox.h"
#include <sstream>
#include <cmath>

StressInBox::StressInBox() : BaseObservable() {}

StressInBox::~StressInBox() {
    if(_out.is_open()) _out.close();
}

void StressInBox::get_settings(input_file &my_inp, input_file &sim_inp) {
    BaseObservable::get_settings(my_inp, sim_inp);

    // side length
    getInputDouble(&my_inp, "L", &_L, 0);

    // center (default 0,0,0). Depending on your parse_input helpers, you may have:
    // - getInputDouble for each component, or
    // - a vector parser. Here’s the safe approach:
    double cx=0., cy=0., cz=0.;
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
           "Initialising StressInBox with center=(%.3f,%.3f,%.3f), L=%.3f, file=%s",
           _center.x, _center.y, _center.z, _L, _tensor_filename.c_str());
}

void StressInBox::init() {
    if(_dump_tensor) {
        _out.open(_tensor_filename, std::ios::out | std::ios::app);
        if(!_out.good()) {
            throw oxDNAException("StressInBox: cannot open '%s' for writing", _tensor_filename.c_str());
        }
        // Optional header
        _out << "# step  xx xy xz  yx yy yz  zx zy zz  N_in_box\n";
        _out.flush();
    }
}

bool StressInBox::require_data_on_CPU() {
    // IMPORTANT:
    // If you run GPU backend, this should force positions/forces to be accessible on CPU.
    // (Exact behavior depends on oxDNA version; but this is the intent.)
    return true;
}

// Minimum-image displacement from center -> particle
LR_vector StressInBox::_pbc_displacement_to_center(const LR_vector& r) const {
    // If your Box class provides minimum image, use it.
    // Many oxDNA versions have something like: _config_info->box->min_image(dr)
    LR_vector dr = r - _center;
    _config_info->box->min_image(dr); // if available in your build
    return dr;
}

bool StressInBox::_inside_box(const LR_vector& r) const {
    const LR_vector dr = _pbc_displacement_to_center(r);
    return (std::fabs(dr.x) <= _halfL &&
            std::fabs(dr.y) <= _halfL &&
            std::fabs(dr.z) <= _halfL);
}

void StressInBox::update_data(llint curr_step) {
    // Access particles
    // In many oxDNA versions CONFIG_INFO has particles in _config_info->particles
    // If yours differs, adapt this block.
    auto &particles = _config_info->particles;

    // Tensor accumulator
    double xx=0., xy=0., xz=0.;
    double yx=0., yy=0., yz=0.;
    double zx=0., zy=0., zz=0.;
    llint n_in = 0;

    for(size_t i = 0; i < particles.size(); i++) {
        BaseParticle *p = particles[i];
        const LR_vector r = p->pos;

        if(!_inside_box(r)) continue;
        n_in++;

        // Total force on particle (already computed by integrator each step)
        const LR_vector F = p->force;

        // Use displacement from center, not absolute coordinate, to avoid origin dependence
        const LR_vector dr = _pbc_displacement_to_center(r);

        // Virial proxy: dr ⊗ F
        xx += dr.x * F.x;  xy += dr.x * F.y;  xz += dr.x * F.z;
        yx += dr.y * F.x;  yy += dr.y * F.y;  yz += dr.y * F.z;
        zx += dr.z * F.x;  zy += dr.z * F.y;  zz += dr.z * F.z;
    }

    // Convert to "stress" by dividing by volume.
    // Sign conventions vary; many MD codes define stress = -virial/V.
    // Pick one and be consistent with your existing global tensor.
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
    // If you want this observable to print into the normal data_output stream instead of a file,
    // you can return a one-line summary here.
    std::stringstream ss;
    ss << curr_step << " # stress_in_box dumped to " << _tensor_filename << "\n";
    return ss.str();
}
