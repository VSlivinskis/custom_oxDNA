/**
 * @file    custom_force.h
 * @date    22/oct/2025
 * @author  Victor
 */

#ifndef AFMMovingSphere_H_
#define AFMMovingSphere_H_

#include <tuple>
#include <vector>
#include <string>

#include "BaseForce.h"

class AFMMovingSphere : public BaseForce {
public:
    AFMMovingSphere();
    ~AFMMovingSphere() override = default;

    // oxDNA hooks
    std::tuple<std::vector<int>, std::string> init(input_file &inp) override;
    LR_vector value(llint step, LR_vector &pos) override;
    number    potential(llint step, LR_vector &pos) override;

private:
    // Interpolated center at a given MD step
    LR_vector center_for_step(llint step) const;

    // Parameters
    number    _stiff{};          // harmonic stiffness
    number    _r0{};             // initial radius
    number    _rate{};           // growth rate (per step)
    number    _r_ext{1e10};      // outer cutoff

    // Positions
    LR_vector _center{0., 0., 0.};  // optional static center if used
    LR_vector _origin{0., 0., 0.};  // start of motion
    LR_vector _target{0., 0., 0.};  // end of motion
    llint     _steps{0};            // steps to go origin -> target

    // NEW: AFM control
    LR_vector _dir{0.,0.,1.};      // approach direction (normalized)
    number    _F_set{0.};          // desired normal force (setpoint, in your force units)
    number    _Kp{0.0};            // proportional gain (nm per force unit)
    number    _max_step{0.0};      // max tip motion per MD step (nm)
    number    _eps{1e-3};          // contact tolerance (nm)
    llint     _log_every{1000};    // debug cadence
    number _min_gap{1e30};   // smallest (distance - radius) seen this step
    int    _min_id{-1};      // id of the particle closest to sphere

    // step-to-step bookkeeping
    llint     _last_step{-1};
    LR_vector _center_curr{0.,0.,0.};
    LR_vector _F_sum{0.,0.,0.};    // accumulated tip force from previous value() calls

    // ----- AFM SCAN (new) -----
    bool     _scan_enabled{false};
    LR_vector _scan_origin_xy{0.,0.,0.};   // lower-left corner (x0,y0)
    LR_vector _scan_size_xy{0.,0.,0.};     // (Lx,Ly)
    int      _Nx{0}, _Ny{0};               // pixels in x,y
    bool     _serp{true};                  // serpentine lines

    llint    _guard{200};                  // steps to wait after lateral move
    llint    _settle{2000};                // steps to let controller settle
    llint    _sample{500};                 // steps to average z / Fn
    llint    _phase_ctr{0};                // steps elapsed within current pixel

    // per-pixel accumulators
    number   _sum_z{0.0}, _sum_Fn{0.0};
    llint    _n_sample{0};
    number   _min_gap_pix{1e30};           // best (closest) gap across pixel

    // raster indices
    int      _ix{0}, _iy{0};               // current pixel
    FILE*    _csv{nullptr};
    std::string _csv_path;

    // helpers
    std::pair<number,number> pixel_xy(int ix, int iy) const;
    void set_lateral(number cx, number cy);
    void advance_pixel_indices();
    void open_csv_if_needed();
    void write_pixel_row(number z_avg, number Fn_avg);

};

#endif // AFMMovingSphere_H_

