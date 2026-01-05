#include "RepulsiveEllipsoid.h"

#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

RepulsiveEllipsoid::RepulsiveEllipsoid() : BaseForce() {
	_r_2    = LR_vector(1., 1., 1.);
	_centre = LR_vector(0., 0., 0.);
	_rate   = 0.;
}

/**
 * Parses:
 *  - stiff (required)
 *  - rate  (optional, default 0)
 *  - r_2   (required): ax,ay,az
 *  - center (optional, default 0,0,0)
 *  - particle (required): id or list or -1
 */
std::tuple<std::vector<int>, std::string> RepulsiveEllipsoid::init(input_file &inp) {
	BaseForce::init(inp);

	getInputNumber(&inp, "stiff", &_stiff, 1);
	getInputNumber(&inp, "rate",  &_rate,  0);

	std::string str;
	double tmpf[3];

	// r_2 (semi-axes at step 0)
	getInputString(&inp, "r_2", str, 1);
	int tmpi = sscanf(str.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) throw oxDNAException("Could not parse r_2 %s", str.c_str());
	_r_2 = LR_vector(tmpf[0], tmpf[1], tmpf[2]);

	// center
	if(getInputString(&inp, "center", str, 0) == KEY_FOUND) {
		tmpi = sscanf(str.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
		if(tmpi != 3) throw oxDNAException("Could not parse center %s", str.c_str());
		_centre = LR_vector(tmpf[0], tmpf[1], tmpf[2]);
	}

	std::string particles_string;
	getInputString(&inp, "particle", particles_string, 1);

	auto particle_ids =
		Utils::get_particles_from_string(
			CONFIG_INFO->particles(),
			particles_string,
			"RepulsiveEllipsoid"
		);

	std::string desc = Utils::sformat(
		"RepulsiveEllipsoid(stiff=%g, rate=%g, axes0=%g,%g,%g, center=%g,%g,%g)",
		_stiff, _rate, _r_2.x, _r_2.y, _r_2.z, _centre.x, _centre.y, _centre.z
	);

	return std::make_tuple(particle_ids, desc);
}

/**
 * CPU mirror of CUDA_REPULSIVE_ELLIPSOID:
 *
 * growth = 1 + rate * step
 * a = a0 * growth, b = b0 * growth, c = c0 * growth
 *
 * r'^2 = (dx/a)^2 + (dy/b)^2 + (dz/c)^2
 * apply WCA only if 0 < r' < 1
 *
 * sigma = 1 / 2^(1/6)
 *
 * pref = 24*eps*(2*s12 - s6) / r'^2
 * Fx = pref * dx / a^2, etc.
 */
LR_vector RepulsiveEllipsoid::value(llint step, LR_vector &pos) {
	// Uniform axis growth, identical to CUDA
	const number growth = (number)1.0 + _rate * (number)step;
	if(growth <= (number)0.0) return LR_vector(0., 0., 0.);

	const number ax = _r_2.x * growth;
	const number ay = _r_2.y * growth;
	const number az = _r_2.z * growth;

	// Guard against degenerate axes
	if(ax <= (number)0.0 || ay <= (number)0.0 || az <= (number)0.0) return LR_vector(0., 0., 0.);

	// Minimum-image displacement from centre to particle (same ordering as CUDA: minimum_image(centre, ppos))
	LR_vector dist = CONFIG_INFO->box->min_image(_centre, pos);

	const number inv_ax2 = (number)1.0 / (ax * ax);
	const number inv_ay2 = (number)1.0 / (ay * ay);
	const number inv_az2 = (number)1.0 / (az * az);

	const number dx = dist.x;
	const number dy = dist.y;
	const number dz = dist.z;

	// r'^2 in ellipsoidal metric
	const number rp2 = dx*dx*inv_ax2 + dy*dy*inv_ay2 + dz*dz*inv_az2;
	if(rp2 <= (number)0.0) return LR_vector(0., 0., 0.);

	const number rp = sqrt(rp2);

	// CUDA cutoff: r' must be < 1.0
	const number Rc_d = (number)1.0;
	if(rp <= (number)0.0 || rp >= Rc_d) return LR_vector(0., 0., 0.);

	// sigma = Rc_d / 2^(1/6) = 1 / 2^(1/6)
	static const number two_to_1_over_6 = pow((number)2.0, (number)(1.0/6.0));
	const number sigma = Rc_d / two_to_1_over_6;

	const number inv_rp    = (number)1.0 / rp;
	const number s_over_rp = sigma * inv_rp;

	const number s2  = s_over_rp * s_over_rp;
	const number s4  = s2 * s2;
	const number s6  = s4 * s2;
	const number s12 = s6 * s6;

	const number eps = _stiff;

	// pref = 24*eps*(2*s12 - s6) / r'^2
	const number pref = (number)24.0 * eps * ((number)2.0 * s12 - s6) / rp2;

	// Chain-rule back to unscaled coordinates (identical to CUDA)
	const number Fx = pref * dx * inv_ax2;
	const number Fy = pref * dy * inv_ay2;
	const number Fz = pref * dz * inv_az2;

	return LR_vector(Fx, Fy, Fz);
}

number RepulsiveEllipsoid::potential(llint step, LR_vector &pos) {
	const number growth = (number)1.0 + _rate * (number)step;
	if(growth <= (number)0.0) return (number)0.0;

	const number ax = _r_2.x * growth;
	const number ay = _r_2.y * growth;
	const number az = _r_2.z * growth;

	if(ax <= (number)0.0 || ay <= (number)0.0 || az <= (number)0.0) return (number)0.0;

	LR_vector dist = CONFIG_INFO->box->min_image(_centre, pos);

	const number inv_ax2 = (number)1.0 / (ax * ax);
	const number inv_ay2 = (number)1.0 / (ay * ay);
	const number inv_az2 = (number)1.0 / (az * az);

	const number dx = dist.x;
	const number dy = dist.y;
	const number dz = dist.z;

	const number rp2 = dx*dx*inv_ax2 + dy*dy*inv_ay2 + dz*dz*inv_az2;
	if(rp2 <= (number)0.0) return (number)0.0;

	const number rp = sqrt(rp2);

	const number Rc_d = (number)1.0;
	if(rp <= (number)0.0 || rp >= Rc_d) return (number)0.0;

	static const number two_to_1_over_6 = pow((number)2.0, (number)(1.0/6.0));
	const number sigma = Rc_d / two_to_1_over_6;

	const number inv_rp    = (number)1.0 / rp;
	const number s_over_rp = sigma * inv_rp;

	const number s2  = s_over_rp * s_over_rp;
	const number s4  = s2 * s2;
	const number s6  = s4 * s2;
	const number s12 = s6 * s6;

	const number eps = _stiff;

	// WCA potential shift (+eps) exactly like your CPU code and standard WCA
	const number U = (number)4.0 * eps * (s12 - s6) + eps;
	return U;
}
