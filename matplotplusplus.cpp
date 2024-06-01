#include <matplot/matplot.h>
#include <vector>

int main() {
    using namespace matplot;

    // Sample data
    std::vector<double> x = {1, 2, 3, 4, 5};
    std::vector<double> y = {1, 4, 9, 16, 25};

    // Create a plot
    plot(x, y);
    title("2D Plot");
    xlabel("X-axis");
    ylabel("Y-axis");

    // Show the plot
    show();

    return 0;
}