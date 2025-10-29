/*
 * custom_force.cpp
 *
 *  Created on: 22/oct/2025
 *      Author: Victor
 */

#include "RepulsiveSphereMoving.h"

#include <tuple>
#include <vector>
#include <string>
#include <cstdio>     // sscanf
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

RepulsiveSphereMoving::RepulsiveSphereMoving() : BaseForce() {
    _r0     = -1.;
    _r_ext  = 1e10;
    _center = LR_vector(0., 0., 0.);
    _rate   = 0.;
    _steps  = 0;
    _origin = LR_vector(0., 0., 0.);
    _target = LR_vector(0., 0., 0.);
}

std::tuple<std::vector<int>, std::string>
RepulsiveSphereMoving::init(input_file &inp) {
    BaseForce::init(inp);

    getInputNumber(&inp, "stiff", &_stiff, 1);
    getInputNumber(&inp, "r0",    &_r0,    1);
    getInputNumber(&inp, "rate",  &_rate,  0);
    getInputNumber(&inp, "r_ext", &_r_ext, 0);

    // Parse center (optional)
    std::string str_center;
    if (getInputString(&inp, "center", str_center, 0) == KEY_FOUND) {
        double tmpf[3];
        if (std::sscanf(str_center.c_str(), "%lf,%lf,%lf", &tmpf[0], &tmpf[1], &tmpf[2]) != 3)
            throw oxDNAException("Could not parse center '%s' in external forces file.", str_center.c_str());
        _center = LR_vector((number)tmpf[0], (number)tmpf[1], (number)tmpf[2]);
    }

    // Parse origin (optional; defaults to _center if not provided)
    std::string str_origin;
    if (getInputString(&inp, "origin", str_origin, 0) == KEY_FOUND) {
        double tmpf[3];
        if (std::sscanf(str_origin.c_str(), "%lf,%lf,%lf", &tmpf[0], &tmpf[1], &tmpf[2]) != 3)
            throw oxDNAException("Could not parse origin '%s' in external forces file.", str_origin.c_str());
        _origin = LR_vector((number)tmpf[0], (number)tmpf[1], (number)tmpf[2]);
    } else {
        _origin = _center;
    }

    // Parse target (optional)
    std::string str_target;
    if (getInputString(&inp, "target", str_target, 0) == KEY_FOUND) {
        double tmpf[3];
        if (std::sscanf(str_target.c_str(), "%lf,%lf,%lf", &tmpf[0], &tmpf[1], &tmpf[2]) != 3)
            throw oxDNAException("Could not parse target '%s' in external forces file.", str_target.c_str());
        _target = LR_vector((number)tmpf[0], (number)tmpf[1], (number)tmpf[2]);
    } // else stays at default (0,0,0)

    // Steps (how many MD steps to go from origin -> target)
    // Read as `number` to use the compiled getInputNumber<number> overload, then cast.
    number steps_num;
    if (getInputNumber(&inp, "steps", &steps_num, 0) == KEY_FOUND) {
        _steps = static_cast<llint>(steps_num);
    } else if (getInputNumber(&inp, "move_steps", &steps_num, 0) == KEY_FOUND) {
        _steps = static_cast<llint>(steps_num);
    } else {
        _steps = 0; // static if not set
    }


    std::string particles_string;
    getInputString(&inp, "particle", particles_string, 1);

    std::string description = Utils::sformat(
        "RepulsiveSphereMoving (stiff=%g, r0=%g, rate=%g, center=%g,%g,%g; origin=%g,%g,%g -> target=%g,%g,%g; steps=%lld)",
        _stiff, _r0, _rate,
        _center.x, _center.y, _center.z,
        _origin.x, _origin.y, _origin.z,
        _target.x, _target.y, _target.z,
        (long long)_steps
    );

    auto particle_ids = Utils::get_particles_from_string(
        CONFIG_INFO->particles(), particles_string, "RepulsiveSphereMoving");

    return std::make_tuple(particle_ids, description);
}

// Helper: linear-interpolated center for this step
LR_vector RepulsiveSphereMoving::center_for_step(llint step) const {
    if (_steps <= 0) return _origin; // static at origin if steps not set
    number t = (number)step / (number)_steps;
    if (t < 0.) t = 0.;
    if (t > 1.) t = 1.;
    return _origin + (_target - _origin) * t;
}

LR_vector RepulsiveSphereMoving::value(llint step, LR_vector &pos) {
    LR_vector c = center_for_step(step);
    LR_vector dist = CONFIG_INFO->box->min_image(c, pos);
    number mdist = dist.module();
	number radius = _r0 + _rate * (number) step;
    double sigma = 0.5;

	if(mdist <= radius || mdist >= _r_ext) return LR_vector(0., 0., 0.);
        else return dist * (_stiff * exp(-pow((mdist - radius), 2) / (2 * pow(sigma, 2))));

}

number RepulsiveSphereMoving::potential(llint step, LR_vector &pos) {
    LR_vector c = center_for_step(step);
    LR_vector dist = CONFIG_INFO->box->min_image(c, pos);
    number mdist = dist.module();
    number radius = _r0 + _rate * (number) step;

    number sigma = 0.5;  // nm
    if(mdist <= radius || mdist >= _r_ext)
        return 0.;

    number term1 = -_stiff * (-pow(sigma, 2) * exp(-pow((mdist - radius), 2) / (2 * pow(sigma, 2))));
    number term2 = -_stiff * (radius * sigma * sqrt(M_PI / 2.0) * erf((mdist - radius) / (sqrt(2.0) * sigma)));
    return term1 + term2;
}
