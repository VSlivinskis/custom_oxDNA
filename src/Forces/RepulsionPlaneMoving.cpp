/*
 * RepulsionPlaneMoving.cpp
 *
 *  Created on: 18/oct/2011
 *      Author: Flavio 
 */

#include "RepulsionPlaneMoving.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"
#include "../Utilities/Utils.h"

#include <algorithm>

RepulsionPlaneMoving::RepulsionPlaneMoving() :
                BaseForce() {
    _particles_string = "-1";
    _ref_particles_string = "-1";
    low_idx = high_idx = -1;
    _origin = LR_vector(0.,0.,0.);
    _target = LR_vector(0.,0.,0.);
    _plane_point = LR_vector(0.,0.,0.);
    _move_steps = 0;
}

std::tuple<std::vector<int>, std::string> RepulsionPlaneMoving::init(input_file &inp) {
    BaseForce::init(inp);

    getInputString(&inp, "particle", _particles_string, 1);
    // ref_particle is now optional; when absent we use ref_position/target
    getInputString(&inp, "ref_particle", _ref_particles_string, 0);

    getInputNumber(&inp, "stiff", &_stiff, 1);

    int tmpi;
    double tmpf[3];
    std::string strdir;
    getInputString(&inp, "dir", strdir, 1);
    tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
    if(tmpi != 3) throw oxDNAException("Could not parse dir %s in external forces file. Aborting", strdir.c_str());
    // _direction = LR_vector(-1, 0, 0);
    _direction = LR_vector(tmpf[0], tmpf[1], tmpf[2]);
    _direction.normalize();

    auto particle_indices_vector = Utils::get_particles_from_string(CONFIG_INFO->particles(), _particles_string, "moving repulsion plane force (RepulsionPlaneMoving.cpp)");
    // Decide which path we use: ref particles vs ref_position/target
    std::vector<int> ref_particle_indices_vector;
    if(!_ref_particles_string.empty() && _ref_particles_string != "-1") {
      ref_particle_indices_vector = Utils::get_particles_from_string(CONFIG_INFO->particles(), _ref_particles_string, "moving repulsion plane force (RepulsionPlaneMoving.cpp)");
    }
    _use_ref = !ref_particle_indices_vector.empty();
  
    if(_use_ref) {
      sort(ref_particle_indices_vector.begin(), ref_particle_indices_vector.end());
      low_idx = ref_particle_indices_vector.front();
      high_idx = ref_particle_indices_vector.back();
      if((high_idx - low_idx + 1) != (int) ref_particle_indices_vector.size())
        throw oxDNAException("RepulsionPlaneMoving requires the list of ref_particle indices to be contiguous");
      for(std::vector<int>::size_type i = 0; i < ref_particle_indices_vector.size(); i++) {
        ref_p_ptr.push_back(CONFIG_INFO->particles()[ref_particle_indices_vector[i]]);
      }
   } else {
      // Parse ref_position (optional, defaults to 0,0,0), target (optional), move_steps (optional)
      std::string strpos, strtgt;
      if(getInputString(&inp, "ref_position", strpos, 0) == 0) {
        if(sscanf(strpos.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2) != 3)
          throw oxDNAException("Could not parse ref_position %s in external forces file. Aborting", strpos.c_str());
        _origin = LR_vector(tmpf[0], tmpf[1], tmpf[2]);
      }
      if(getInputString(&inp, "target", strtgt, 0) == 0) {
        if(sscanf(strtgt.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2) != 3)
          throw oxDNAException("Could not parse target %s in external forces file. Aborting", strtgt.c_str());
        _target = LR_vector(tmpf[0], tmpf[1], tmpf[2]);
      } else {
        _target = _origin; // default: static plane
      }
      number moves_num = 0.0;
      getInputNumber(&inp, "move_steps", &moves_num, 0);   // oxDNA expects number*
      // clamp to [0, +inf) and cast to llint
      long long moves_ll = (long long) llround(moves_num);
      if(moves_ll < 0) moves_ll = 0;
      _move_steps = (llint) moves_ll;
      _plane_point = _origin; // initialize
   }

   std::string description;
   if(_use_ref) {
      description = Utils::sformat("RepulsionPlaneMoving (ref_particle) stiff=%g, n=(%g,%g,%g)", _stiff, _direction.x, _direction.y, _direction.z);
   } else {
      description = Utils::sformat("RepulsionPlaneMoving (ref_position→target) stiff=%g, n=(%g,%g,%g), origin=(%g,%g,%g), target=(%g,%g,%g), steps=%lld",
                                   _stiff, _direction.x, _direction.y, _direction.z,
                                   _origin.x, _origin.y, _origin.z,
                                   _target.x, _target.y, _target.z,
                                   (long long)_move_steps);
   }

    return std::make_tuple(particle_indices_vector, description);
}
inline void RepulsionPlaneMoving::update_plane_point(llint step) {
    if(_use_ref) {
    // Use average absolute position of the reference particles as the plane anchor
    if(ref_p_ptr.empty()) return;
    LR_vector c(0.,0.,0.);
    for(auto *p : ref_p_ptr) c += CONFIG_INFO->box->get_abs_pos(p);
    c /= (number)ref_p_ptr.size();
    _plane_point = c;
    return;
  }
  // Linear interpolation from _origin to _target over _move_steps
  if(_move_steps <= 0) {
    _plane_point = _origin;
    return;
  }
  llint s = std::max<llint>(0, std::min(step, _move_steps));
  number alpha = (number)s / (number)_move_steps;
  _plane_point = _origin + (_target - _origin) * alpha;
}

LR_vector RepulsionPlaneMoving::value(llint step, LR_vector &pos) {
    update_plane_point(step);
    LR_vector force(0., 0., 0.);
    number distance = (pos - _plane_point) * _direction;      // signed distance along normal
    if(distance < (number)0.) {
        force = -_stiff * _direction * distance;                // push back into allowed half-space
    }
    return force;
}

number RepulsionPlaneMoving::potential(llint step, LR_vector &pos) {
    update_plane_point(step);
    number distance = (pos - _plane_point) * _direction;
    if(distance < (number)0.) {
        return (number)0.5 * _stiff * distance * distance;
    }
    return (number)0.;
}