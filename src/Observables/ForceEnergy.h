/*
 * ForceEnergy.h
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#ifndef FORCEENERGY_H_
#define FORCEENERGY_H_

#include "BaseObservable.h"

#include <string>
#include <vector>

/**
 * @brief Outputs the energy of the system due to external forces acting upon particles.
 *
 * Optional parameters:
 * @verbatim
 [print_group = <string> (limits the energy computation to forces belonging to a specific group)]
 [per_particle = <bool> (if true, output "particle_id:U,particle_id:U,..." for each selected particle)]
 [particles = <string> (particle selection. "-1" means all. Otherwise supports "0,1,2,10-15")]
 * @endverbatim
 */

class ForceEnergy: public BaseObservable {
protected:
	std::string _group_name;
	bool _per_particle = false;

	// selection: if empty => all
	std::vector<int> _particle_ids;

	// helpers
	static std::vector<int> _parse_particles_spec(const std::string &spec);
	static std::string _trim(const std::string &s);

public:
	ForceEnergy();
	virtual ~ForceEnergy();

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init();

	std::string get_output_string(llint curr_step);
};

#endif /* FORCEENERGY_H_ */
