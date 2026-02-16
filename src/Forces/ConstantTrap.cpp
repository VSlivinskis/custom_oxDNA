/*
 * ConstantTrap.cpp
 *
 *  Created on: 18/oct/2011
 *      Author: Flavio
 */

#include "ConstantTrap.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

// constant force between a particle and either:
//  - another particle (ref_particle >= 0), or
//  - a fixed point in space (ref_particle = -1, with center_x/y/z provided)

ConstantTrap::ConstantTrap() :
				BaseForce() {
	_p_ptr = NULL;
	PBC = false;
	_r0 = -1.;
	_ref_id = -2;

	_use_fixed_point = false;
	_ref_pos = LR_vector(0., 0., 0.);
}

std::tuple<std::vector<int>, std::string> ConstantTrap::init(input_file &inp) {
	BaseForce::init(inp);

	int particle;
	getInputInt(&inp, "particle", &particle, 1);

	// Make ref_particle optional:
	//   ref_particle >= 0 : reference is that particle
	//   ref_particle = -1 : reference is a fixed point (center_x/y/z)
	_ref_id = -1;
	getInputInt(&inp, "ref_particle", &_ref_id, 0);

	getInputNumber(&inp, "r0", &_r0, 1);
	getInputNumber(&inp, "stiff", &_stiff, 1);
	getInputBool(&inp, "PBC", &PBC, 0);

	int N = CONFIG_INFO->particles().size();

	if (_ref_id == -1) {
		_use_fixed_point = true;

		number cx, cy, cz;
		getInputNumber(&inp, "center_x", &cx, 1);
		getInputNumber(&inp, "center_y", &cy, 1);
		getInputNumber(&inp, "center_z", &cz, 1);
		_ref_pos = LR_vector(cx, cy, cz);
	}
	else {
		_use_fixed_point = false;

		if(_ref_id < 0 || _ref_id >= N)
			throw oxDNAException("Invalid reference particle %d for ConstantTrap", _ref_id);

		_p_ptr = CONFIG_INFO->particles()[_ref_id];
	}

	if(particle >= N || particle < -1)
		throw oxDNAException("Trying to add a ConstantTrap on non-existent particle %d. Aborting", particle);
	if(particle == -1)
		throw oxDNAException("Cannot apply ConstantTrap to all particles. Aborting");

	std::string description;
	if (_use_fixed_point) {
		description = Utils::sformat(
			"ConstantTrap (stiff=%g, r0=%g, ref_particle=-1, center=(%g,%g,%g), PBC=%d)",
			_stiff, _r0, _ref_pos.x, _ref_pos.y, _ref_pos.z, PBC
		);
	} else {
		description = Utils::sformat(
			"ConstantTrap (stiff=%g, r0=%g, ref_particle=%d, PBC=%d)",
			_stiff, _r0, _ref_id, PBC
		);
	}

	return std::make_tuple(std::vector<int> {particle}, description);
}

LR_vector ConstantTrap::_distance(LR_vector u, LR_vector v) {
	if(PBC) {
		return CONFIG_INFO->box->min_image(u, v);
	}
	else {
		return v - u;
	}
}

LR_vector ConstantTrap::value(llint step, LR_vector &pos) {
	LR_vector ref = _use_fixed_point ? _ref_pos : CONFIG_INFO->box->get_abs_pos(_p_ptr);
	LR_vector dr = _distance(pos, ref); // ref - self

	number r = dr.module();
	if (r == 0.) {
		// undefined direction; return no force to avoid NaNs
		return LR_vector(0., 0., 0.);
	}

	number sign = copysign(1., (double) (r - _r0));
	return (_stiff * sign) * (dr / r);
}

number ConstantTrap::potential(llint step, LR_vector &pos) {
	LR_vector ref = _use_fixed_point ? _ref_pos : CONFIG_INFO->box->get_abs_pos(_p_ptr);
	LR_vector dr = _distance(pos, ref); // ref - self
	return _stiff * fabs((dr.module() - _r0));
}
