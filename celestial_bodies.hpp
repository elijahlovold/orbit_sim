#pragma once

#include <iostream>
#include <matplotlibcpp.h>
#include <vector>
#include <math.h>
#include <string>
// #include <cblas.h>

#define G 6.67430e-11   // Gravitational constant (m^3 kg^-1 s^-2)
#define c 299792458     // Speed of light (m/s)

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
        static constexpr int rec_decmi = 10;    // how many samples to decimate before recording
        static constexpr double ep = 1000;

        cel_type type;
        std::string name;
        std::string color;

        double m; 
        
        double x;
        double y;
        double z;

        double dx_dt;
        double dy_dt;
        double dz_dt;

        double dx2_2dt;
        double dy2_2dt;
        double dz2_2dt;


        int inc;
        std::vector<double> history_x_coords;
        std::vector<double> history_y_coords;
        std::vector<double> history_z_coords;

        celestial(std::string ic_name, double ic_m, double ic_x, double ic_y, double ic_z, double ic_dx_dt, double ic_dy_dt, double ic_dz_dt);

        void accum_accel(celestial *cj); 
        void step_position(std::vector<celestial *> c_);

};


// double celestial::compute_grav_attract_pp(double r, double m0, double m1) {
//     return G*m0*m1/(r*r);
// }

celestial::celestial(std::string ic_name, double ic_m, double ic_x, double ic_y, double ic_z, double ic_dx_dt, double ic_dy_dt, double ic_dz_dt) {
    this->name = ic_name;
    this->color = "def";

    this->m = ic_m;

    this->x = ic_x;
    this->y = ic_y;
    this->z = ic_z;

    this->dx_dt = ic_dx_dt;
    this->dy_dt = ic_dy_dt;
    this->dz_dt = ic_dz_dt;

    this->inc = 0;
}

void celestial::accum_accel(celestial *cj) {
    double x_diff = this->x - cj->x;
    double y_diff = this->y - cj->y;
    double z_diff = this->z - cj->z;

    double r_sq = x_diff*x_diff + y_diff*y_diff + z_diff*z_diff;
    double r_pow_3 = pow(r_sq, 1.5);

    // if on top, sigularity, don't compute
    if (r_pow_3 == 0) {
        return;
    }

    this->dx2_2dt -= G*cj->m*x_diff/r_pow_3; 
    this->dy2_2dt -= G*cj->m*y_diff/r_pow_3; 
    this->dz2_2dt -= G*cj->m*z_diff/r_pow_3; 
}

void celestial::step_position(std::vector<celestial *> c_) {
    // reset accelerations
    this->dx2_2dt = 0; 
    this->dy2_2dt = 0; 
    this->dz2_2dt = 0; 

    // sum accelerations
    for (celestial *ci : c_) {
        this->accum_accel(ci);
    }

    // step velocities
    this->dx_dt += this->dx2_2dt*ep;
    this->dy_dt += this->dy2_2dt*ep;
    this->dz_dt += this->dz2_2dt*ep;
    
    // step positions
    this->x += this->dx_dt*ep;
    this->y += this->dy_dt*ep;
    this->z += this->dz_dt*ep;

    if (this->inc++ % rec_decmi == 0) {
        this->history_x_coords.push_back(this->x);
        this->history_y_coords.push_back(this->y);
        this->history_z_coords.push_back(this->z);
    }
}