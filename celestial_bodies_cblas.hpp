#pragma once

#include <iostream>
#include <matplotlibcpp.h>
#include <vector>
#include <math.h>
#include <string>
#include <thread>
#include <chrono>
#include <cblas.h>

#define G 6.67430e-11   // Gravitational constant (m^3 kg^-1 s^-2)
#define c 299792458     // Speed of light (m/s)
#define N 3             // number of deminsions

namespace plt = matplotlibcpp;

// Circular orbit: e = 0
// Elliptic orbit: 0 < e < 1
// Parabolic trajectory: e = 1
// Hyperbolic trajectory: e > 1

// e = 1 - 2 / (ra/rp + 1)
// e = (ra - rp) / (ra + rp)
// ra = aphelion
// rp = perihelion

enum cel_type {
    PLANET, 
    MOON, 
    STAR,
    ASTERIOD, 
    COMET, 
    DWARF_PLANET, 
    GALAXY, 
    BLACK_HOLE, 
    BLAZAR, 
    QUASER, 
    EXOPLANET, 
    NEUTRON_STAR,
    PULSAR, 
};


class celestial {
    public: 
        static constexpr int rec_decmi = 100;    // how many samples to decimate before recording
        static constexpr double ep = 1000;
 
        cel_type type;
        std::string name;
        std::string color;

        double m; 
        double drag;

        double pos[N];

        bool record;
        int inc;
        std::vector<double> history_x_coords;
        std::vector<double> history_y_coords;
        std::vector<double> history_z_coords;


        celestial(std::string ic_name, double ic_m, double ic_x, double ic_y, double ic_z, double ic_dx_dt, double ic_dy_dt, double ic_dz_dt);
        celestial(std::string ic_name, double ic_m, double coords_ics[N], double vel_ics[N]);
        celestial(std::string ic_name, double ic_M, double ic_m, double ic_r);

        void accum_accel(celestial *cj); 
        void step_position(std::vector<celestial *> c_);
        void step_position(std::vector<celestial *> c_, double *origin_x, double *origin_y, double *origin_z); 

    private: 
        // static constexpr double ep = 100;
        double vel[N];
        double accel[N];
};

celestial::celestial(std::string ic_name, double ic_m, double ic_x, double ic_y, double ic_z, double ic_dx_dt, double ic_dy_dt, double ic_dz_dt) {
    this->name = ic_name;
    this->color = "def";

    this->m = ic_m;
    this->drag = 0;

    this->pos[0] = ic_x;
    this->pos[1] = ic_y;
    this->pos[2] = ic_z;

    this->vel[0] = ic_dx_dt;
    this->vel[1] = ic_dy_dt;
    this->vel[2] = ic_dz_dt;


    this->inc = 0;
    this->record = true;
}

celestial::celestial(std::string ic_name, double ic_M, double ic_m, double ic_r) {
    this->name = ic_name;
    this->color = "def";

    this->m = ic_m;
    this->drag = 0;

    this->pos[0] = ic_r;
    this->pos[1] = 0;
    this->pos[2] = 0;

    this->vel[0] = 0;
    this->vel[1] = sqrt(G*ic_M/ic_r);
    this->vel[2] = 0;

    this->inc = 0;
    this->record = true;
}

celestial::celestial(std::string ic_name, double ic_m, double coords_ics[N], double vel_ics[N]) {
    this->name = ic_name;
    this->color = "def";

    this->m = ic_m;
    this->drag = 0;

    cblas_dcopy(N, coords_ics, 1, this->pos, 1);
    cblas_dcopy(N, vel_ics, 1, this->vel, 1);

    this->inc = 0;
    this->record = true;
}

void celestial::accum_accel(celestial *cj) {
    double diff[N];     // diff = this - cj
    cblas_dcopy(N, this->pos, 1, diff, 1);
    cblas_daxpy(N, -1, cj->pos, 1, diff, 1);

    double r_sq = cblas_ddot(N, diff, 1, diff, 1);
    double r_pow_3 = pow(r_sq, 1.5);

    // if on top, sigularity, don't compute
    if (r_pow_3 == 0) {
        return;
    }

    // accel -= G*m*diff/r_pow_3; 
    cblas_daxpy(N, -G*cj->m/r_pow_3, diff, 1, this->accel, 1);
}

void celestial::step_position(std::vector<celestial *> c_) {
    // reset accelerations
    cblas_dscal(N, 0.0, accel, 1);

    // sum accelerations
    for (celestial *ci : c_) {
        this->accum_accel(ci);
    }

    // step velocities
    cblas_daxpy(N, ep, this->accel, 1, this->vel, 1);
    
    // step positions
    cblas_daxpy(N, ep, this->vel, 1, this->pos, 1);

    if (this->record && (this->inc++ % rec_decmi == 0)) {
        this->history_x_coords.push_back(this->pos[0]);
        this->history_y_coords.push_back(this->pos[1]);
        this->history_z_coords.push_back(this->pos[2]);
    }
}