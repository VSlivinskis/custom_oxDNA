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
    // Vector center -> particle with PBC
    LR_vector dist = CONFIG_INFO->box->min_image(_center, pos);
    number d = dist.module();

    // Desired spherical cutoff radius where U=F=0 (can grow with "rate")
    const number Rc = _r0 + _rate * (number) step;
    if (d <= 0.0 || d >= Rc) {
        return LR_vector(0., 0., 0.);
    }

    // Map cutoff to WCA sigma: rc = 2^(1/6) * sigma  =>  sigma = Rc / 2^(1/6)
    static const number two_to_1_over_6 = pow(2.0, 1.0/6.0);
    const number sigma = Rc / two_to_1_over_6;

    // Precompute powers
    const number inv_d  = 1.0 / d;
    const number s_over_d = sigma * inv_d;
    const number s6 = pow(s_over_d, 6);
    const number s12 = s6 * s6;

    // epsilon from "stiff"
    const number eps = _stiff;

    // Force magnitude for LJ (same as WCA inside cutoff), outward (repulsive)
    // F = 24*eps*(2*sigma^12/d^13 - sigma^6/d^7)
    const number Fmag = 24.0 * eps * (2.0 * s12 - s6) * inv_d;

    // Direction (normalize dist)
    LR_vector dir = dist * inv_d;
    return dir * Fmag;
}

number RepulsiveSphere::potential(llint step, LR_vector &pos) {
    LR_vector dist = CONFIG_INFO->box->min_image(_center, pos);
    number d = dist.module();

    const number Rc = _r0 + _rate * (number) step;
    if (d <= 0.0 || d >= Rc) {
        return 0.0;
    }

    static const number two_to_1_over_6 = pow(2.0, 1.0/6.0);
    const number sigma = Rc / two_to_1_over_6;

    const number inv_d  = 1.0 / d;
    const number s_over_d = sigma * inv_d;
    const number s6  = pow(s_over_d, 6);
    const number s12 = s6 * s6;

    const number eps = _stiff;

    // U_LJ shifted so U(Rc)=0: U = 4*eps*(s12 - s6) + eps   (for d < Rc)
    const number U = 4.0 * eps * (s12 - s6) + eps;
    return U;
}
