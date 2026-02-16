/*
 * ForceEnergy.cpp
 *
 *  Created on: Oct 1, 2013
 *      Author: rovigatti
 */

#include "ForceEnergy.h"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>

ForceEnergy::ForceEnergy() {
}

ForceEnergy::~ForceEnergy() {
}

std::string ForceEnergy::_trim(const std::string &s) {
	size_t i = 0;
	while(i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) i++;
	size_t j = s.size();
	while(j > i && std::isspace(static_cast<unsigned char>(s[j - 1]))) j--;
	return s.substr(i, j - i);
}

std::vector<int> ForceEnergy::_parse_particles_spec(const std::string &spec_in) {
	// Supports:
	//   "-1" => all (represented by empty vector here; caller handles it)
	//   "0,1,2"
	//   "10-15"
	//   "0, 1, 2, 10-15"
	std::string spec = _trim(spec_in);
	if(spec.empty() || spec == "-1") return {};

	std::vector<int> out;
	std::stringstream ss(spec);
	std::string token;

	while(std::getline(ss, token, ',')) {
		token = _trim(token);
		if(token.empty()) continue;

		auto dash = token.find('-');
		if(dash == std::string::npos) {
			// single int
			int id = std::stoi(token);
			out.push_back(id);
		}
		else {
			// range a-b
			std::string a_s = _trim(token.substr(0, dash));
			std::string b_s = _trim(token.substr(dash + 1));
			int a = std::stoi(a_s);
			int b = std::stoi(b_s);
			if(b < a) std::swap(a, b);
			for(int id = a; id <= b; id++) out.push_back(id);
		}
	}

	// de-dup + sort
	std::sort(out.begin(), out.end());
	out.erase(std::unique(out.begin(), out.end()), out.end());
	return out;
}

void ForceEnergy::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputString(&my_inp, "print_group", _group_name, 0);

	int tmp = 0;
	getInputBoolAsInt(&my_inp, "per_particle", &tmp, 0);
	_per_particle = bool(tmp);

	// particles selection: read as string so we can support ranges and CSV
	std::string particles_spec = "-1";
	getInputString(&my_inp, "particles", particles_spec, 0);
	_particle_ids = _parse_particles_spec(particles_spec);
}

void ForceEnergy::init() {
	BaseObservable::init();

	// validate selected ids (if any)
	if(!_particle_ids.empty()) {
		for(int id : _particle_ids) {
			if(id < 0 || id >= _config_info->N()) {
				throw oxDNAException("ForceEnergy: invalid particle id %d", id);
			}
		}
	}
}

std::string ForceEnergy::get_output_string(llint curr_step) {
	// Decide which particles we are iterating
	const bool all = _particle_ids.empty();
	const int n_selected = all ? _config_info->N() : (int)_particle_ids.size();

	if(n_selected <= 0) {
		// nothing selected; return 0 in legacy mode or empty in per-particle
		return _per_particle ? std::string("") : Utils::sformat(_number_formatter, (number)0.f);
	}

	number U_sum = (number)0.f;

	if(_per_particle) {
		std::ostringstream out;
		bool first = true;

		auto emit_one = [&](int pid, number Ui) {
			std::string Ui_s = Utils::sformat(_number_formatter, Ui);
			Ui_s = _trim(Ui_s);

			if(!first) out << ",";
			first = false;
			out << pid << ":" << Ui_s;
		};

		if(all) {
			for(int pid = 0; pid < _config_info->N(); pid++) {
				auto p = _config_info->particles()[pid];
				number Ui = (number)0.f;

				if(_group_name == "") {
					p->set_ext_potential(curr_step, _config_info->box);
					Ui = p->ext_potential;
				}
				else {
					LR_vector abs_pos = _config_info->box->get_abs_pos(p);
					for(auto ext_force : p->ext_forces) {
						if(ext_force->get_group_name() == _group_name) {
							Ui += ext_force->potential(curr_step, abs_pos);
						}
					}
				}

				emit_one(pid, Ui);
			}
		}
		else {
			for(int pid : _particle_ids) {
				auto p = _config_info->particles()[pid];
				number Ui = (number)0.f;

				if(_group_name == "") {
					p->set_ext_potential(curr_step, _config_info->box);
					Ui = p->ext_potential;
				}
				else {
					LR_vector abs_pos = _config_info->box->get_abs_pos(p);
					for(auto ext_force : p->ext_forces) {
						if(ext_force->get_group_name() == _group_name) {
							Ui += ext_force->potential(curr_step, abs_pos);
						}
					}
				}

				emit_one(pid, Ui);
			}
		}

		return out.str();
	}

	// legacy scalar mode: average over selected particles (not always N())
	if(all) {
		for(auto p : _config_info->particles()) {
			if(_group_name == "") {
				p->set_ext_potential(curr_step, _config_info->box);
				U_sum += p->ext_potential;
			}
			else {
				LR_vector abs_pos = _config_info->box->get_abs_pos(p);
				for(auto ext_force : p->ext_forces) {
					if(ext_force->get_group_name() == _group_name) {
						U_sum += ext_force->potential(curr_step, abs_pos);
					}
				}
			}
		}
	}
	else {
		for(int pid : _particle_ids) {
			auto p = _config_info->particles()[pid];
			if(_group_name == "") {
				p->set_ext_potential(curr_step, _config_info->box);
				U_sum += p->ext_potential;
			}
			else {
				LR_vector abs_pos = _config_info->box->get_abs_pos(p);
				for(auto ext_force : p->ext_forces) {
					if(ext_force->get_group_name() == _group_name) {
						U_sum += ext_force->potential(curr_step, abs_pos);
					}
				}
			}
		}
	}

	number U_avg = U_sum / (number)n_selected;
	return Utils::sformat(_number_formatter, U_avg);
}
