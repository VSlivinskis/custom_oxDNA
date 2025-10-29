/*
 * RepulsiveSphere.cpp
 *
 *  Created on: 28/nov/2014
 *      Author: Flavio
 */

#include "RepulsiveSphere.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

RepulsiveSphere::RepulsiveSphere() :
				BaseForce() {
	_r0 = -1.;
	_r_ext = 1e10;
	_center = LR_vector(0., 0., 0.);
	_rate = 0.;
}

std::tuple<std::vector<int>, std::string> RepulsiveSphere::init(input_file &inp) {
	BaseForce::init(inp);

	std::cerr << "[DEBUG] RepulsiveSphere::init() called from "
          << __FILE__ << ":" << __LINE__ << " in " << __func__ << std::endl;

	getInputNumber(&inp, "stiff", &_stiff, 1);
	getInputNumber(&inp, "r0", &_r0, 1);
	getInputNumber(&inp, "rate", &_rate, 0);
	getInputNumber(&inp, "r_ext", &_r_ext, 0);

	std::string particles_string;
	getInputString(&inp, "particle", particles_string, 1);

	std::string strdir;
	if(getInputString(&inp, "center", strdir, 0) == KEY_FOUND) {
		double tmpf[3];
		int tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
		if(tmpi != 3) throw oxDNAException("Could not parse center %s in external forces file. Aborting", strdir.c_str());
		_center = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	}

	std::string description = Utils::sformat("RepulsiveSphere (stiff=%g, r0=%g, rate=%g, center=%g,%g,%g)", _stiff, _r0, _rate, _center.x, _center.y, _center.z);
	auto particle_ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles_string, "RepulsiveSphere");

	return std::make_tuple(particle_ids, description);
}

LR_vector RepulsiveSphere::value(llint step, LR_vector &pos) {
    // Vector from center to particle
    LR_vector dist = CONFIG_INFO->box->min_image(_center, pos);
    number d = dist.module(); // distance from tip center

    // Current tip radius (can move/grow with time)
    number R = _r0 + _rate * (number) step;

    // No force if we're outside (or exactly on) the tip surface
    if (d >= R || d <= 0.0) {
        return LR_vector(0., 0., 0.);
    }

    // penetration fraction x in [0,1] for d in [R,0]
    number x = 1.0 - d / R; // how far "inside" the tip you are

    // cubic penalty with zero slope at boundary:
    // U = K * x^3  =>  F_mag = (3K/R) * x^2
    number K = _stiff;  // rename conceptually
    number F_mag = (3.0 * K / R) * x * x;

    // Direction: outward from center (repulsive core)
    LR_vector dir = dist * (1.0 / d);

    return dir * F_mag;
}

number RepulsiveSphere::potential(llint step, LR_vector &pos) {
    LR_vector dist = CONFIG_INFO->box->min_image(_center, pos);
    number d = dist.module();

    number R = _r0 + _rate * (number) step;

    if (d >= R || d <= 0.0) {
        return 0.0;
    }

    number x = 1.0 - d / R;
    number K = _stiff;

    // U = K * x^3, continuous and 0 at boundary
    number U = K * x * x * x;
    return U;
}
