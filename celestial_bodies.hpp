#pragma once

#include <iostream>
#include <matplotlibcpp.h>
#include <vector>
#include <math.h>
#include <string>
#include <thread>
#include <chrono>

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
        static constexpr int rec_decmi = 100;    // how many samples to decimate before recording
        static constexpr double ep = 1000;
 
        cel_type type;
        std::string name;
        std::string color;

        double m; 
        double drag;

        double x;
        double y;
        double z;

        bool record;
        int inc;
        std::vector<double> history_x_coords;
        std::vector<double> history_y_coords;
        std::vector<double> history_z_coords;

        celestial(std::string ic_name, double ic_m, double ic_x, double ic_y, double ic_z, double ic_dx_dt, double ic_dy_dt, double ic_dz_dt);
        celestial(std::string ic_name, double ic_M, double ic_m, double ic_r);

        void accum_accel(celestial *cj); 
        void step_position(std::vector<celestial *> c_);
        void step_position(std::vector<celestial *> c_, double *origin_x, double *origin_y, double *origin_z); 
        void set_xyz(double set_x, double set_y, double set_z); 
        void set_dxyzdt(double dx_dt, double dy_dt, double dz_dt);

    private: 
        // static constexpr double ep = 100;

        double dx_dt;
        double dy_dt;
        double dz_dt;

        double dx2_2dt;
        double dy2_2dt;
        double dz2_2dt;

};

celestial::celestial(std::string ic_name, double ic_M, double ic_m, double ic_r) {
    this->name = ic_name;
    this->color = "def";

    this->m = ic_m;
    this->drag = 0;

    this->set_xyz(ic_r, 0, 0);
    this->set_dxyzdt(0, sqrt(G*ic_M/ic_r), 0);

    this->inc = 0;
    this->record = true;
}

celestial::celestial(std::string ic_name, double ic_m, double ic_x, double ic_y, double ic_z, double ic_dx_dt, double ic_dy_dt, double ic_dz_dt) {
    this->name = ic_name;
    this->color = "def";

    this->m = ic_m;
    this->drag = 0;

    this->set_xyz(ic_x, ic_y, ic_z);
    this->set_dxyzdt(ic_dx_dt, ic_dy_dt, ic_dz_dt);

    this->inc = 0;
    this->record = true;
}

void celestial::set_xyz(double set_x, double set_y, double set_z) {
    this->x = set_x;
    this->y = set_y;
    this->z = set_z;
}

void celestial::set_dxyzdt(double dx_dt, double dy_dt, double dz_dt) {
    this->dx_dt = dx_dt;
    this->dy_dt = dy_dt;
    this->dz_dt = dz_dt;
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

void celestial::step_position(std::vector<celestial *> c_, double *origin_x, double *origin_y, double *origin_z) {
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

    if (this->record && (this->inc++ % rec_decmi == 0)) {
        this->history_x_coords.push_back(this->x - *origin_x);
        this->history_y_coords.push_back(this->y - *origin_y);
        this->history_z_coords.push_back(this->z - *origin_z);
    }
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

    if (this->record && (this->inc++ % rec_decmi == 0)) {
        this->history_x_coords.push_back(this->x);
        this->history_y_coords.push_back(this->y);
        this->history_z_coords.push_back(this->z);
    }
}