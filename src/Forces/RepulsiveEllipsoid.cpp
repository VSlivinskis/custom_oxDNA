#include "RepulsiveEllipsoid.h"

#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

RepulsiveEllipsoid::RepulsiveEllipsoid() :
    BaseForce()
{
    _r_2 = LR_vector(1.,1.,1.);
    _centre = LR_vector(0.,0.,0.);
}

std::tuple<std::vector<int>, std::string> RepulsiveEllipsoid::init(input_file &inp)
{
    BaseForce::init(inp);

    getInputNumber(&inp, "stiff", &_stiff, 1);

    std::string strdir;
    double tmpf[3];

    // r_2 (semi-axes)
    getInputString(&inp, "r_2", strdir, 1);
    int tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf",
                      tmpf, tmpf+1, tmpf+2);
    if(tmpi != 3)
        throw oxDNAException("Could not parse r_2 %s", strdir.c_str());
    _r_2 = LR_vector(tmpf[0], tmpf[1], tmpf[2]);

    // center
    if(getInputString(&inp, "center", strdir, 0) == KEY_FOUND){
        tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf",
                      tmpf, tmpf+1, tmpf+2);
        if(tmpi != 3)
            throw oxDNAException("Could not parse center %s", strdir.c_str());
        _centre = LR_vector(tmpf[0], tmpf[1], tmpf[2]);
    }

    std::string particles_string;
    getInputString(&inp, "particle", particles_string, 1);

    auto particle_ids =
        Utils::get_particles_from_string(CONFIG_INFO->particles(),
                                         particles_string,
                                         "RepulsiveEllipsoid");

    std::string desc = Utils::sformat(
        "RepulsiveEllipsoid(stiff=%g, axes=%g,%g,%g, center=%g,%g,%g)",
        _stiff, _r_2.x, _r_2.y, _r_2.z, _centre.x, _centre.y, _centre.z
    );

    return std::make_tuple(particle_ids, desc);
}


/********************************************************************
 *  Compute ellipsoidal radius:
 *      d_ellip = sqrt((dx/ax)^2 + (dy/ay)^2 + (dz/az)^2)
 *
 *  Then apply *exact* WCA/LJ 6–12 repulsion as in RepulsiveSphere:
 *      sigma = Rc / 2^(1/6)
 *      U = 4 eps (s12 - s6) + eps
 *      F = derivative (cannot use the del symbol)
 ********************************************************************/

// derivative of ellipsoidal radius direction
static inline LR_vector ellip_normal(const LR_vector &dist, const LR_vector &axes)
{
    // gradient of sqrt( (x/ax)^2 + ... )
    number gx = dist.x / (axes.x*axes.x);
    number gy = dist.y / (axes.y*axes.y);
    number gz = dist.z / (axes.z*axes.z);

    LR_vector g(gx, gy, gz);
    number gm = g.module();
    if(gm > 0) g /= gm;  // unit normal
    return g;
}


LR_vector RepulsiveEllipsoid::value(llint step, LR_vector &pos)
{
    LR_vector dist = CONFIG_INFO->box->min_image(_centre, pos);

    // ellipsoidal scaled distance
    number d =
        sqrt(SQR(dist.x/_r_2.x) +
             SQR(dist.y/_r_2.y) +
             SQR(dist.z/_r_2.z));

    // the "cutoff" in ellipsoidal space is 1.0
    const number Rc = 1.0;
    if(d <= 0.0 || d >= Rc)
        return LR_vector(0.,0.,0.);

    // WCA parameters identical to RepulsiveSphere
    static const number two_to_1_over_6 = pow(2.0, 1.0/6.0);
    number sigma = Rc / two_to_1_over_6;

    const number inv_d = 1.0/d;
    const number s_over_d = sigma*inv_d;
    const number s6 = pow(s_over_d, 6);
    const number s12 = s6*s6;

    const number eps = _stiff;

    // magnitude of derivative dU/dd
    number Fmag = 24.0 * eps * (2.0*s12 - s6) * inv_d;

    // convert scalar radial derivative into force vector:
    LR_vector normal = ellip_normal(dist, _r_2);

    return normal * Fmag;
}


number RepulsiveEllipsoid::potential(llint step, LR_vector &pos)
{
    LR_vector dist = CONFIG_INFO->box->min_image(_centre, pos);

    number d =
        sqrt(SQR(dist.x/_r_2.x) +
             SQR(dist.y/_r_2.y) +
             SQR(dist.z/_r_2.z));

    const number Rc = 1.0;
    if(d <= 0.0 || d >= Rc)
        return 0.0;

    static const number two_to_1_over_6 = pow(2.0,1.0/6.0);
    number sigma = Rc / two_to_1_over_6;

    number inv_d = 1.0/d;
    number s_over_d = sigma*inv_d;
    number s6 = pow(s_over_d,6);
    number s12 = s6*s6;

    const number eps = _stiff;

    number U = 4.0*eps*(s12 - s6) + eps;  // WCA shift
    return U;
}
