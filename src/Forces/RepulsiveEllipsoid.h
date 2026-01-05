/**
 * @file    RepulsiveEllipsoid.h
 * @date    26/mag/2021
 * @author  Andrea
 *
 */

#ifndef REPULSIVEELLIPSOID_H_
#define REPULSIVEELLIPSOID_H_

#include "BaseForce.h"

/**
 * @brief External force field that confines particles into a (possibly growing) ellipsoid.
 *
 * This CPU implementation mirrors the CUDA implementation (CUDA_REPULSIVE_ELLIPSOID):
 *  - Semi-axes grow uniformly with: growth = 1 + rate * step
 *  - Ellipsoidal radius r' = sqrt((dx/a)^2 + (dy/b)^2 + (dz/c)^2)
 *  - WCA repulsion is applied ONLY for r' < 1.0 (purely repulsive, truncated)
 *
 * Input (external forces file):
 *
 * {
 *   particle = -1
 *   type = repulsive_ellipsoid   # or whatever token you map to this class on CPU
 *   r_2 = ax,ay,az               # base semi-axes at step 0
 *   stiff = 10.0                 # epsilon
 *   rate = 0.0                   # linear growth factor per step (dimensionless)
 *   center = 0,0,0
 * }
 */
class RepulsiveEllipsoid: public BaseForce {
public:
	/// center of the ellipsoid
	LR_vector _centre;

	/// base semi-axes at step 0
	LR_vector _r_2;

	/// (kept for compatibility; unused)
	LR_vector _r_1;

	/// linear growth rate (dimensionless per step): growth = 1 + rate * step
	number _rate;

	RepulsiveEllipsoid();
	virtual ~RepulsiveEllipsoid() {}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	LR_vector value(llint step, LR_vector &pos) override;
	number potential(llint step, LR_vector &pos) override;
};

#endif // REPULSIVEELLIPSOID_H_
