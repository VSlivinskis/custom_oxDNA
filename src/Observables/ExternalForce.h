/*
 * ExternalForce.h
 *
 *  Created on: Feb 02, 2024
 *      Author: rovigatti
 */

#ifndef EXTERNALFORCE_H_
#define EXTERNALFORCE_H_

#include "BaseObservable.h"

/**
 * @brief Outputs per-particle forces for selected particles.
 *
 * Output format (one line per particle):
 *   step  particle_id  fx  fy  fz
 *
 * Notes:
 * - Set `particles = -1` to select all particles (handled explicitly).
 * - `update_every` controls the sampling cadence; `print_every` controls output cadence.
 * - This observable reads p->force as computed by the simulation (does not modify state).
 */
class ExternalForce : public BaseObservable {
protected:
    std::string _particle_ids;
    std::vector<BaseParticle *> _particles;

    // Per-particle accumulated forces across updates between prints.
    std::vector<LR_vector> _force_averages;

    bool _printed_header = false;

public:
    ExternalForce();
    virtual ~ExternalForce();

    void get_settings(input_file &my_inp, input_file &sim_inp) override;
    void init() override;

    void update_data(llint curr_step) override;

    std::string get_output_string(llint curr_step) override;
};

#endif /* EXTERNALFORCE_H_ */
