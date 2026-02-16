/**
 * @file    ConstantTrap.h
 * @date    18/oct/2014
 * @author  Flavio
 *
 */

#ifndef CONSTANTTRAP_H_
#define CONSTANTTRAP_H_

#include "BaseForce.h"

class BaseParticle;

/**
 * ConstantTrap
 *
 * Applies a constant-magnitude force along the line connecting the target particle to a reference.
 *
 * Original behavior: reference is another particle (ref_particle).
 * Extended behavior: if ref_particle = -1, use a fixed reference point in space given by
 *   center_x, center_y, center_z
 *
 * The force magnitude is |_stiff| and its sign is determined by (r - r0):
 *   - if r > r0, force points toward the reference
 *   - if r < r0, force points away from the reference
 */
class ConstantTrap: public BaseForce {
private:
	int _ref_id;
	bool _use_fixed_point;
	LR_vector _ref_pos;

public:
	BaseParticle *_p_ptr;
	number _r0;
	bool PBC;

	ConstantTrap();
	virtual ~ConstantTrap() {
	}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);

protected:
	LR_vector _distance(LR_vector u, LR_vector v);
};

#endif // CONSTANTTRAP_H
