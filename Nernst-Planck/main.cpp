#include "model.hpp"

int main() {
    int NX = 100;
    int NY = 30;
    double tau_h = 1.0;
    double tau_g = 1.0;
    double F = 1.0;
    double kappa = 1.0e-1;
    double Deff = 1.0e-2;
    double RT = 3.0;
    int steps = 10000;

    Solver solver(NX, NY, tau_h, tau_g, F, kappa, Deff, RT);
    for (int i{}; i < steps; i++) {
        solver.update_moments();
        solver.do_step_h();
    }

    for (int i{}; i < steps; i++) {
        solver.update_moments();
        solver.do_step_g();
        solver.do_step_h();
    }
    solver.save("result_phi.csv", "result_C.csv");
}