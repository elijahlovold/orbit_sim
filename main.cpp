
#include "celestial_bodies.hpp"




int main() {
    celestial sun("sun", 1.9885e30, 0, 0, 0, 0, 0, 0); 

    celestial earth("earth", 5.972e24, 150e9, 0, 0, 0, 29.8e3, 0); 
    celestial moon("moon", 7.342e22, 150e9 + 380e6, 0, 0, 0, 29.8e3 + 1.022e3, 0); 

    celestial jupiter("jupiter", 1.898e27, 778e9, 0, 0, 0, 13.07e3, 0); 

    std::vector<celestial *> objs = {&earth, &moon, &sun, &jupiter}; 

    // celestial earth(5.972e24, 0, 0, 0, 0, 0, 0); 
    // celestial moon(7.342e22, 380e6, 0, 0, 0, 1.022e3, 0); 

    // std::vector<celestial *> objs = {&earth, &moon}; 
    
    for (int i = 0; i < 100000; i++) {
        for (celestial *obj : objs) {
            obj->step_position(objs);
        }
    }

    // Create a plot
    plt::figure();
    for (celestial *obj : objs) {
        plt::plot(obj->history_x_coords, obj->history_y_coords, {{"label", obj->name}});
    }
    plt::title("Sample Plot");
    plt::xlabel("X-axis");
    plt::ylabel("Y-axis");
    plt::legend();

    // Display the plot
    plt::show();

    return 0;

}