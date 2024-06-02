
#include "celestial_bodies.hpp"
// #include "celestial_bodies_cblas.hpp"




int main(int argc, char *argv[]) {

    int iterations;

    if (argc > 1) {
        try {
            // Convert the first argument to an integer
            iterations = std::stoi(argv[1]);
            std::cout << "Number of iterations: " << iterations << std::endl;
        } catch (const std::invalid_argument &e) {
            std::cerr << "Invalid argument: " << argv[1] << " is not a valid number." << std::endl;
            return 1; // Return an error code
        } catch (const std::out_of_range &e) {
            std::cerr << "Argument out of range: " << argv[1] << std::endl;
            return 1; // Return an error code
        }
        std::cout << argv[1] << std::endl;
    } else {
        iterations = 100000;
    }

    auto time_s = std::chrono::high_resolution_clock::now();

    double *origin_x;
    double *origin_y;
    double *origin_z;

    celestial sun("sun", 1.9885e30, 0, 0, 0, 0, 0, 0); 

    celestial mercury("mercury", sun.m, 3.301e23, 57.91e9); 
    celestial venus("venus", sun.m, 4.867e24, 108.21e9); 
    celestial earth("earth", sun.m, 5.972e24, 150e9); 
    celestial moon("moon", 7.342e22, 150e9 + 380e6, 0, 0, 0, 29.8e3 + 1.02e3, 0); 
    celestial mars("mars", sun.m, 6.417e23, 227.94e9); 
    celestial jupiter("jupiter", sun.m, 1.898e27, 778e9); 
    celestial saturn("saturn", sun.m, 5.683e26, 1433.53e9); 
    celestial uranus("uranus", sun.m, 8.681e25, 2870.97e9); 
    celestial neptune("neptune", sun.m, 1.024e26, 4500e9); 

    // origin_x = &earth.x;
    // origin_y = &earth.y;
    // origin_z = &earth.z;

    std::vector<celestial *> objs = {&sun, &mercury, &venus, &earth, &moon, &mars, &jupiter, &saturn, &uranus, &neptune}; 
    // std::vector<celestial *> objs = {&earth, &moon}; 
    
    for (int i = 0; i < iterations; i++) {
        for (celestial *obj : objs) {
            obj->step_position(objs);
            // obj->step_position(objs, origin_x, origin_y, origin_z);
        }
    }

    auto time_e = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration = (time_e - time_s);

    std::cout << "decmi: " << celestial::rec_decmi << std::endl;
    std::cout << "time step: " << celestial::ep << std::endl;
    std::cout << "time: " << duration.count() << " seconds" << std::endl;

    // Create a plot
    plt::figure();
    for (celestial *obj : objs) {
        if (obj->record) {
            plt::plot(obj->history_x_coords, obj->history_y_coords, {{"label", obj->name}});
        }
    }
    plt::title("Sample Plot");
    plt::xlabel("X-axis");
    plt::ylabel("Y-axis");
    plt::legend();

    // Display the plot
    plt::show();

    return 0;

}


// std
// decmi: 100
// time step: 1000
// time: 1.02233 seconds

// cblas
// decmi: 100
// time step: 1000
// time: 1.2567 seconds