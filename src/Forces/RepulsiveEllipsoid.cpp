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

LR_vector RepulsiveEllipsoid::value(llint step, LR_vector &pos) {
    // Uniform axis growth
    const number growth = (number)1.0 + _rate * (number)step;
    if(growth <= (number)0.0) return LR_vector(0., 0., 0.);

    const number ax = _r_2.x * growth;
    const number ay = _r_2.y * growth;
    const number az = _r_2.z * growth;

    if(ax <= (number)0.0 || ay <= (number)0.0 || az <= (number)0.0) return LR_vector(0., 0., 0.);

    // Minimum-image displacement from centre to particle
    LR_vector dist = CONFIG_INFO->box->min_image(_centre, pos);

    const number dx = dist.x;
    const number dy = dist.y;
    const number dz = dist.z;

    const number inv_ax2 = (number)1.0 / (ax * ax);
    const number inv_ay2 = (number)1.0 / (ay * ay);
    const number inv_az2 = (number)1.0 / (az * az);

    // Ellipsoidal metric r'^2
    const number rp2 = dx*dx*inv_ax2 + dy*dy*inv_ay2 + dz*dz*inv_az2;
    if(rp2 <= (number)0.0) return LR_vector(0., 0., 0.);

    const number rp = sqrt(rp2);
    if(rp <= (number)0.0) return LR_vector(0., 0., 0.);

    // --- WCA in r' ---
    const number x        = (number)2;         // same as your moving sphere default
    const number sigma    = (number)1;         // same as moving sphere default
    const number epsilon  = _stiff;

    // cutoff in r' (WCA cutoff)
    const number rc = pow((number)2.0, (number)(1.0/x)) * sigma;

    // If you want the same "only inside r'<1" behavior as your CUDA mirror, clamp cutoff:
    // const number rc = (number)1.0;

    if(rp >= rc) return LR_vector(0., 0., 0.);

    // Clamp rp to avoid blowup at rp->0
    const number rp_safe = (rp > (number)1e-9) ? rp : (number)1e-9;

    // U(r') = 4ε[(σ/r')^(2x) - (σ/r')^x] + ε  for r'<rc
    // Let A = (σ/r')^x
    const number A = pow(sigma / rp_safe, x);

    // dU/dr' = 4ε(2A - 1) dA/dr' , dA/dr' = -(x/r') A
    const number dUdrp = (number)4.0 * epsilon * ((number)2.0 * A - (number)1.0) * (-(x / rp_safe) * A);

    // Force: F = -dU/dr' * grad(r')
    // grad(r') = (1/r') * (dx/a^2, dy/b^2, dz/c^2)
    const number scale = (-dUdrp) / rp_safe;

    const number Fx = scale * dx * inv_ax2;
    const number Fy = scale * dy * inv_ay2;
    const number Fz = scale * dz * inv_az2;

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

    const number dx = dist.x;
    const number dy = dist.y;
    const number dz = dist.z;

    const number inv_ax2 = (number)1.0 / (ax * ax);
    const number inv_ay2 = (number)1.0 / (ay * ay);
    const number inv_az2 = (number)1.0 / (az * az);

    const number rp2 = dx*dx*inv_ax2 + dy*dy*inv_ay2 + dz*dz*inv_az2;
    if(rp2 <= (number)0.0) return (number)0.0;

    const number rp = sqrt(rp2);
    if(rp <= (number)0.0) return (number)0.0;

    // --- WCA in r' ---
    const number x        = (number)2;
    const number sigma    = (number)1;
    const number epsilon  = _stiff;

    const number rc = pow((number)2.0, (number)(1.0/x)) * sigma;

    // If you want the same "only inside r'<1" behavior as your CUDA mirror, clamp cutoff:
    // const number rc = (number)1.0;

    if(rp >= rc) return (number)0.0;

    const number rp_safe = (rp > (number)1e-9) ? rp : (number)1e-9;

    const number A = pow(sigma / rp_safe, x);
    const number U = (number)4.0 * epsilon * (A*A - A) + A*sigma;
    return U;
}
