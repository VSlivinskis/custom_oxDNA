/*
 * custom_force.cpp
 *
 *  Created on: 22/oct/2025
 *      Author: Victor
 */

#include "AFMMovingSphere.h"

#include <cmath>   // exp, pow, sqrt, erf, M_PI
#include <tuple>
#include <vector>
#include <string>
#include <numeric> // <— add at the top for std::iota
#include <cstdio>     // sscanf
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

AFMMovingSphere::AFMMovingSphere() :
                BaseForce() {
    _r0     = -1.;
    _r_ext  = 1e10;
    _center = LR_vector(0., 0., 0.);
    _rate   = 0.;
    _steps  = 0;
    _origin = LR_vector(0., 0., 0.);
    _target = LR_vector(0., 0., 0.);
}

std::tuple<std::vector<int>, std::string> AFMMovingSphere::init(input_file &inp) {

    // 1) Parse common/base stuff and CAPTURE the particle IDs
    std::vector<int> particle_ids;
    std::string base_desc;
    std::tie(particle_ids, base_desc) = BaseForce::init(inp);   // <-- keep these

    // ----------------- Minimal backstop for -1/all or a single id -----------------
    // Only if BaseForce::init() didn’t populate anything:
    if (particle_ids.empty()) {
        number pnum = 0.0;
        // Accept either key name
        int found = getInputNumber(&inp, "particle", &pnum, 0);
        if (found != KEY_FOUND) {
            found = getInputNumber(&inp, "particles", &pnum, 0);
        }
        if (found == KEY_FOUND) {
            if ((int)pnum == -1) {
                // Expand to all particles
                int N = 0;

                // Preferred: ConfigInfo exposes N() as a method
                N = (int)CONFIG_INFO->N();

                if (N <= 0) {
                    throw oxDNAException("AFMMovingSphere: could not resolve particle list (N unknown).");
                }

                particle_ids.resize(N);
                std::iota(particle_ids.begin(), particle_ids.end(), 0);
            } else {
                particle_ids.push_back((int)pnum);
            }
        }
    }
    // ------------------------------------------------------------------------------

    // 2) Your existing parameter parsing...
    getInputNumber(&inp, "stiff", &_stiff, 1);
    getInputNumber(&inp, "r0",   &_r0,   0);
    getInputNumber(&inp, "rate", &_rate, 0);
    getInputNumber(&inp, "r_ext",&_r_ext, 0);

    // Optional direct center (legacy name: 'center')
    std::string str_center;
    if (getInputString(&inp, "center", str_center, 0) == KEY_FOUND) {
        double tmpf[3];
        if (std::sscanf(str_center.c_str(), "%lf,%lf,%lf", &tmpf[0], &tmpf[1], &tmpf[2]) != 3)
            throw oxDNAException("AFMMovingSphere: could not parse center '%s'.", str_center.c_str());
        _center = LR_vector((number)tmpf[0], (number)tmpf[1], (number)tmpf[2]);
    } else {
        _center = LR_vector(0.,0.,0.);
    }

    // ---- Motion definition: origin/target/steps ----
    bool have_origin = false;
    std::string str_origin;
    if (getInputString(&inp, "origin", str_origin, 0) == KEY_FOUND) {
        double tmpf[3];
        if (std::sscanf(str_origin.c_str(), "%lf,%lf,%lf", &tmpf[0], &tmpf[1], &tmpf[2]) != 3)
            throw oxDNAException("AFMMovingSphere: could not parse origin '%s'.", str_origin.c_str());
        _origin = LR_vector((number)tmpf[0], (number)tmpf[1], (number)tmpf[2]);
        have_origin = true;
    }
    std::string str_refpos;
    if (!have_origin && getInputString(&inp, "ref_position", str_refpos, 0) == KEY_FOUND) {
        double tmpf[3];
        if (std::sscanf(str_refpos.c_str(), "%lf,%lf,%lf", &tmpf[0], &tmpf[1], &tmpf[2]) != 3)
            throw oxDNAException("AFMMovingSphere: could not parse ref_position '%s'.", str_refpos.c_str());
        _origin = LR_vector((number)tmpf[0], (number)tmpf[1], (number)tmpf[2]);
        have_origin = true;
    }
    if (!have_origin) _origin = _center;

    std::string str_target;
    if (getInputString(&inp, "target", str_target, 0) == KEY_FOUND) {
        double tmpf[3];
        if (std::sscanf(str_target.c_str(), "%lf,%lf,%lf", &tmpf[0], &tmpf[1], &tmpf[2]) != 3)
            throw oxDNAException("AFMMovingSphere: could not parse target '%s'.", str_target.c_str());
        _target = LR_vector((number)tmpf[0], (number)tmpf[1], (number)tmpf[2]);
    } else {
        _target = _origin;
    }

    // move_steps / steps
    number steps_tmp = 0.0;
    if (getInputNumber(&inp, "move_steps", &steps_tmp, 0) == KEY_FOUND) {
        _steps = (llint)steps_tmp;
    } else if (getInputNumber(&inp, "steps", &steps_tmp, 0) == KEY_FOUND) {
        _steps = (llint)steps_tmp;
    } else {
        _steps = 0;
    }

    // AFM controller
    _dir = LR_vector(0.,0.,1.);
    std::string str_dir;
    if (getInputString(&inp, "dir", str_dir, 0) == KEY_FOUND) {
        double d[3];
        if (std::sscanf(str_dir.c_str(), "%lf,%lf,%lf", &d[0], &d[1], &d[2]) != 3)
            throw oxDNAException("AFMMovingSphere: could not parse dir '%s'.", str_dir.c_str());
        _dir = LR_vector((number)d[0], (number)d[1], (number)d[2]);
        number n = _dir.module();
        if (n == 0.) throw oxDNAException("AFMMovingSphere: dir cannot be zero.");
        _dir /= n;
    }

    getInputNumber(&inp, "F_set",    &_F_set,    0);
    getInputNumber(&inp, "Kp",       &_Kp,       0);
    getInputNumber(&inp, "max_step", &_max_step, 0);
    getInputNumber(&inp, "eps",      &_eps,      0);
    number logn = (number)_log_every;
    if (getInputNumber(&inp, "log_every", &logn, 0) == KEY_FOUND) {
        if (logn < 0) logn = 0;
        _log_every = (llint)logn;
    }

    // per-step state
    _center_curr = _origin;
    _last_step   = -1;
    _F_sum       = LR_vector(0.,0.,0.);
    _min_gap     = 1e30;
    _min_id      = -1;

    // ---- AFM scan (optional) ----
    std::string s_origin_xy, s_size_xy, s_pixels, s_csv;
    double tmp[3];

    // origin (x0,y0)
    if (getInputString(&inp, "scan_origin_xy", s_origin_xy, 0) == KEY_FOUND) {
        if (std::sscanf(s_origin_xy.c_str(), "%lf,%lf", &tmp[0], &tmp[1]) != 2)
            throw oxDNAException("AFMMovingSphere: bad scan_origin_xy '%s'.", s_origin_xy.c_str());
        _scan_origin_xy = LR_vector((number)tmp[0], (number)tmp[1], 0.);
    }

    // size (Lx,Ly)
    if (getInputString(&inp, "scan_size_xy", s_size_xy, 0) == KEY_FOUND) {
        if (std::sscanf(s_size_xy.c_str(), "%lf,%lf", &tmp[0], &tmp[1]) != 2)
            throw oxDNAException("AFMMovingSphere: bad scan_size_xy '%s'.", s_size_xy.c_str());
        _scan_size_xy = LR_vector((number)tmp[0], (number)tmp[1], 0.);
    }

    // pixels Nx,Ny
    int Nx=0, Ny=0;
    if (getInputString(&inp, "scan_pixels", s_pixels, 0) == KEY_FOUND) {
        if (std::sscanf(s_pixels.c_str(), "%d,%d", &Nx, &Ny) != 2 || Nx < 1 || Ny < 1)
            throw oxDNAException("AFMMovingSphere: bad scan_pixels '%s'.", s_pixels.c_str());
        _Nx = Nx; _Ny = Ny;
    }

    // serpentine (0/1)
    number serp = 1.0;
    if (getInputNumber(&inp, "scan_serpentine", &serp, 0) == KEY_FOUND) {
        _serp = (serp != 0.0);
    }

    // phases
    number n;
    if (getInputNumber(&inp, "pixel_step_guard", &n, 0) == KEY_FOUND) _guard  = (llint) std::max(0.0, n);
    if (getInputNumber(&inp, "pixel_settle_steps", &n, 0) == KEY_FOUND) _settle = (llint) std::max(0.0, n);
    if (getInputNumber(&inp, "pixel_sample_steps", &n, 0) == KEY_FOUND) _sample = (llint) std::max(1.0, n);

    // CSV path
    if (getInputString(&inp, "save_csv", s_csv, 0) == KEY_FOUND) _csv_path = s_csv;

    // enable scan if configured
    _scan_enabled = (_Nx > 0 && _Ny > 0 && _scan_size_xy.x >= 0.0 && _scan_size_xy.y >= 0.0);

    // init raster state
    _ix = 0; _iy = 0; _phase_ctr = 0;
    _sum_z = _sum_Fn = 0.0; _n_sample = 0; _min_gap_pix = 1e30;

    // position tip at the first pixel's lateral location; keep current z
    if (_scan_enabled) {
        open_csv_if_needed();
        auto xy = pixel_xy(_ix, _iy);
        set_lateral(xy.first, xy.second);
    }

    if (_Kp > 0.0 && _max_step <= 0.0) {
        // optional warning
    }
    if (_r0 < 0.0) {
        throw oxDNAException("AFMMovingSphere: r0 must be >= 0.");
    }

    // Final checks
    if (particle_ids.empty()) {
        throw oxDNAException("AFMMovingSphere: no particles selected. "
                             "Did you set 'particle = -1' or a valid list?");
    }
    return std::make_tuple(particle_ids, std::string("AFMMovingSphere"));
}


// Helper: linear-interpolated center for this step
LR_vector AFMMovingSphere::center_for_step(llint step) const {
    if (_steps <= 0) return _origin; // static at origin if steps not set
    number t = (number)step / (number)_steps;
    if (t < 0.) t = 0.;
    if (t > 1.) t = 1.;
    return _origin + (_target - _origin) * t;
}

std::pair<number,number> AFMMovingSphere::pixel_xy(int ix, int iy) const {
    // Map indices to coordinates covering the FULL rectangle [x0, x0+Lx] and [y0, y0+Ly]
    number x0 = _scan_origin_xy.x, y0 = _scan_origin_xy.y;
    number Lx = _scan_size_xy.x,  Ly = _scan_size_xy.y;

    // serpentine in x (so successive rows reverse direction)
    int jx = ix;
    if (_serp && (iy % 2 == 1)) jx = (_Nx - 1) - ix;

    number dx = (_Nx > 1) ? (Lx / (number)(_Nx - 1)) : 0.0;
    number dy = (_Ny > 1) ? (Ly / (number)(_Ny - 1)) : 0.0;

    number cx = x0 + jx * dx;
    number cy = y0 + iy * dy;
    return {cx, cy};
}

void AFMMovingSphere::set_lateral(number cx, number cy) {
    // Keep current z; only change x/y of the tip center
    _center_curr.x = cx;
    _center_curr.y = cy;
}

void AFMMovingSphere::advance_pixel_indices() {
    bool new_row = false;

    _ix++;
    if (_ix >= _Nx) {
        _ix = 0;
        _iy++;
        new_row = true;
    }

    _phase_ctr = 0;
    _sum_z = _sum_Fn = 0.0;
    _n_sample = 0;
    _min_gap_pix = 1e30;

    if (_iy < _Ny) {
        auto xy = pixel_xy(_ix, _iy);
        set_lateral(xy.first, xy.second);

        if (new_row) {
            // optional "lift" or "reset" so each new line starts from a sane height
            // idea 1: move the tip back up some safe offset along +dir
            _center_curr += _dir * 5.0; // <-- tunable safety lift in your length units
        }
    } else {
        _scan_enabled = false;
        if (_csv) { std::fflush(_csv); std::fclose(_csv); _csv = nullptr; }
    }
}


void AFMMovingSphere::open_csv_if_needed() {
    if (_csv || _csv_path.empty()) return;
    _csv = std::fopen(_csv_path.c_str(), "w");
    if (_csv) {
        std::fprintf(_csv, "ix,iy,x,y,z_avg,Fn_avg,min_gap\n");
        std::fflush(_csv);
    }
}

void AFMMovingSphere::write_pixel_row(number z_avg, number Fn_avg) {
    if (!_csv) return;
    auto xy = pixel_xy(_ix, _iy);
    std::fprintf(_csv, "%d,%d,%.9g,%.9g,%.9g,%.9g,%.9g\n",
                 _ix, _iy, (double)xy.first, (double)xy.second,
                 (double)z_avg, (double)Fn_avg, (double)_min_gap_pix);
    std::fflush(_csv);  // <-- this part is CRITICAL
}


LR_vector AFMMovingSphere::value(llint step, LR_vector &pos) {
    if (step != _last_step) {
        if (_scan_enabled && _last_step >= 0) {
            _phase_ctr++;  // move phase forward at the start of the new MD step
        }
        // ----- finalize previous step (if any) -----
        if (_last_step >= 0) {
            // Compute finalized normal force BEFORE moving the tip
            number F_normal_prev = _F_sum * _dir;
            number z_prev = _center_curr * _dir;  // projection along approach direction

            // If scanning and we're inside the SAMPLE window for the current pixel, accumulate
            if (_scan_enabled) {
                llint g = _guard, s = _settle, m = _sample;
                llint k = _phase_ctr; // 0..(g+s+m-1)
                if (k >= (g + s) && k < (g + s + m)) {
                    _sum_z  += z_prev;
                    _sum_Fn += F_normal_prev;
                    _n_sample++;
                }
            }

            // ---- Constant-force controller with contact safety ----
            if (_Kp > 0.0) {
                // proposed motion along the AFM approach direction
                number dz = _Kp * (_F_set - F_normal_prev);

                // clamp dz magnitude
                if (_max_step > 0.0) {
                    if (dz >  _max_step) dz =  _max_step;
                    if (dz < -_max_step) dz = -_max_step;
                }

                // contact detection: _min_gap < 0 means we're overlapping the sample
                bool in_contact = (_min_gap <= 0.0);

                // Interpret positive dz as "move further along +_dir".
                // We forbid *further penetration* once contact exists.
                if (in_contact && dz > 0.0) {
                    dz = 0.0;
                }

                _center_curr += _dir * dz;
            }

            // Optional periodic logging
            if (_log_every > 0 && (_last_step % _log_every) == 0) {
                // You can print z_prev, F_normal_prev, _min_gap here if desired
            }

            // ---- Advance scan phase and finish pixel if needed ----
            if (_scan_enabled) {
                _phase_ctr++;
                llint total = _guard + _settle + _sample;
                if (_phase_ctr >= total) {
                    // finalize pixel
                    number z_avg  = (_n_sample > 0) ? (_sum_z  / (number)_n_sample) : (_center_curr * _dir);
                    number Fn_avg = (_n_sample > 0) ? (_sum_Fn / (number)_n_sample) : (F_normal_prev);
                    write_pixel_row(z_avg, Fn_avg);
                    // std::cout << "[AFM PIXEL DONE] ix=" << _ix
                    //         << " iy=" << _iy
                    //         << " z_avg=" << z_avg
                    //         << " Fn_avg=" << Fn_avg
                    //         << " min_gap=" << _min_gap_pix
                    //         << std::endl;
                    advance_pixel_indices();  // moves lateral to next pixel and resets per-pixel accumulators
                }
            }
        }

        // ----- reset accumulators for THIS step -----
        _F_sum     = LR_vector(0.,0.,0.);
        _min_gap   = 1e30;
        _min_id    = -1;
        _last_step = step;

        // If controller disabled, keep your scripted motion
        if (_Kp == 0.0) {
            _center_curr = center_for_step(step);
        }

        std::cout << "[AFM STEP] step=" << step
                  << " center=("
                  << _center_curr.x << ","
                  << _center_curr.y << ","
                  << _center_curr.z << ")"
                  << " min_gap=" << _min_gap
                  << std::endl;
    }

    // ---- Compute per-particle repulsion ----
    LR_vector c = _center_curr;
    LR_vector dist = CONFIG_INFO->box->min_image(c, pos);
    number mdist = dist.module();
    number radius = _r0 + _rate * (number) step;
    // number sigma = 0.5;

    number overlap = radius - mdist; // positive if particle is INSIDE the sphere

    LR_vector fi(0.,0.,0.);

    // Case A: particle inside sphere -> strong push outward
    if (overlap > 0.0) {
        // Hooke-like contact: F = k * overlap in the outward normal
        number k_contact = _stiff; // reuse _stiff but now it means contact stiffness
        LR_vector n_hat = dist / (mdist + 1e-12); // unit normal (center -> particle)
        fi = n_hat * (k_contact * overlap);

    // Case B: otherwise, if you still want a "soft halo" before contact:
    } else if (mdist < _r_ext) {
        number gap = -overlap; // = mdist - radius
        number sigma = 0.5;
        number fmag = _stiff * exp( - (gap*gap) / (2*sigma*sigma) );
        LR_vector n_hat = dist / (mdist + 1e-12);
        fi = n_hat * fmag;
    }

    // Accumulate reaction force on tip
    _F_sum -= fi;

    // Track minimum gap (contact detection)
    number gap = mdist - radius;
    if (gap < _min_gap) _min_gap = gap;
    if (gap < _min_gap_pix) _min_gap_pix = gap;   // <-- add this

    return fi;
}


number AFMMovingSphere::potential(llint step, LR_vector &pos) {
    // Use the same effective center as in value():
    LR_vector c = (_Kp == 0.0) ? center_for_step(step) : _center_curr;

    // Relative position and distance
    LR_vector dist = CONFIG_INFO->box->min_image(c, pos);
    number mdist = dist.module();

    // Sphere geometry
    number radius = _r0 + _rate * (number) step;
    number overlap = radius - mdist; // >0 means particle is inside the tip
    number sigma = 0.5;

    // Outside the cutoff? no interaction
    if (mdist >= _r_ext) {
        return 0.0;
    }

    // Case A: inside the tip (contact region) -> quadratic penalty
    if (overlap > 0.0) {
        number k_contact = _stiff; // same stiffness you used for fi in value()
        number U = 0.5 * k_contact * overlap * overlap;
        return U;
    }

    // Case B: just outside the tip surface -> soft Gaussian halo
    // gap = mdist - radius = -overlap >= 0
    number gap = -overlap;
    number U = _stiff * std::exp( - (gap * gap) / (2.0 * sigma * sigma) );
    return U;
}
