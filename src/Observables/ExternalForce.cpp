/*
 * ExternalForce.cpp
 *
 *  Created on: Feb 2, 2024
 *      Author: rovigatti
 */

#include "ExternalForce.h"

ExternalForce::ExternalForce() {
}

ExternalForce::~ExternalForce() {
}

void ExternalForce::get_settings(input_file &my_inp, input_file &sim_inp) {
    BaseObservable::get_settings(my_inp, sim_inp);
    getInputString(&my_inp, "particles", _particle_ids, 0);
}

void ExternalForce::init() {
    std::vector<int> ids = Utils::get_particles_from_string(
        _config_info->particles(), _particle_ids, "ExternalForce observable"
    );

    // Treat "-1" as "all particles"
    if (ids.size() == 1 && ids[0] == -1) {
        ids.clear();
    }

    _particles.clear();

    if (ids.size() == 0) {
        _particles = _config_info->particles();
    } else {
        for (auto p : _config_info->particles()) {
            if (std::find(ids.begin(), ids.end(), p->index) != std::end(ids)) {
                _particles.push_back(p);
            }
        }
    }

    _force_averages.assign(_particles.size(), LR_vector());
    _times_updated = 0;
    _printed_header = false;
}

void ExternalForce::update_data(llint /*curr_step*/) {
    // IMPORTANT:
    // Do NOT call p->set_initial_forces() here.
    // Observables must not mutate the simulation state, and calling it would reset forces.
    //
    // We read p->force as computed by the simulation at the time observables are updated.
    for (uint i = 0; i < _particles.size(); i++) {
        _force_averages[i] += _particles[i]->force;
    }
    _times_updated++;
}

std::string ExternalForce::get_output_string(llint curr_step) {
    std::string output;

    if (!_printed_header) {
        output += "# step particle_id fx fy fz\n";
        _printed_header = true;
    }

    // If update_every was not set (or update never triggered), avoid divide-by-zero.
    // In this situation, we print the instantaneous forces currently stored.
    const double denom = (_times_updated > 0) ? (double)_times_updated : 1.0;

    for (uint i = 0; i < _particles.size(); i++) {
        LR_vector avg = _force_averages[i] / denom;

        output += Utils::sformat(
            "%lld %d %lf %lf %lf\n",
            curr_step,
            _particles[i]->index,
            avg.x, avg.y, avg.z
        );

        _force_averages[i] = LR_vector();
    }

    _times_updated = 0;
    return output;
}
