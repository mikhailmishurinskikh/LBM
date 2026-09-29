#include <fstream>
#include <string>

#define STAY 0
#define RIGHT 1
#define UP 2
#define LEFT 3
#define DOWN 4
#define UP_RIGHT 5
#define UP_LEFT 6
#define DOWN_LEFT 7
#define DOWN_RIGHT 8

inline void swap(double*& ptr1, double*& ptr2) {
    double* tmp = ptr1;
    ptr1 = ptr2;
    ptr2 = tmp;
}

class Solver {
    /*
    D2Q9:
    Directions numeration:
    6  2  5
    3  0  1
    7  4  8

    */
    int NX {100};
    int NY {50};

    double* f_new;
    double* f_old;

    double* rho;
    double* ux;
    double* uy;

    double dt;
    double tau;
    double omega;
    double time;

    double Fx;
    double Fy;

public:

    Solver(int nx, int ny, double fx, double p_tau) : NX(nx), NY(ny), Fx(fx), tau(p_tau)
    {
        dt = 1.0;
        omega = dt / tau;
        time = 0;
        f_old = new double[NX*NY*9];
        f_new = new double[NX*NY*9];
        rho = new double[NX*NY];
        ux = new double[NX*NY];
        uy = new double[NX*NY];
        Fy = 0;

        init();
    }

    ~Solver() {
        delete[] f_old;
        delete[] f_new;
        delete[] rho;
        delete[] ux;
        delete[] uy;
    }

private:

    inline int idx_vec(int x, int y, int n) const {
        return (x * NY + y) * 9 + n;
    }

    inline int idx_scal(int x, int y) const {
        return x * NY + y;
    }    

    void init() {
        for (int nx{0}; nx < NX; nx++) {
            for (int ny{0}; ny < NY; ny++) {
                f_old[idx_vec(nx, ny, 0)] = (4.0/9.0);

                f_old[idx_vec(nx, ny, 1)] = (1.0/9.0);
                f_old[idx_vec(nx, ny, 2)] = (1.0/9.0);
                f_old[idx_vec(nx, ny, 3)] = (1.0/9.0);
                f_old[idx_vec(nx, ny, 4)] = (1.0/9.0);

                f_old[idx_vec(nx, ny, 5)] = (1.0/36.0);
                f_old[idx_vec(nx, ny, 6)] = (1.0/36.0);
                f_old[idx_vec(nx, ny, 7)] = (1.0/36.0);
                f_old[idx_vec(nx, ny, 8)] = (1.0/36.0);
            }
        }
    }

    void calc_rho() {
        for (int nx{0}; nx < NX; nx++) {
            for (int ny{0}; ny < NY; ny++) {
                rho[idx_scal(nx, ny)] = f_old[idx_vec(nx, ny, STAY)]
                            + f_old[idx_vec(nx, ny, RIGHT)]
                            + f_old[idx_vec(nx, ny, UP)]
                            + f_old[idx_vec(nx, ny, LEFT)]
                            + f_old[idx_vec(nx, ny, DOWN)]
                            + f_old[idx_vec(nx, ny, UP_RIGHT)]
                            + f_old[idx_vec(nx, ny, UP_LEFT)]
                            + f_old[idx_vec(nx, ny, DOWN_LEFT)]
                            + f_old[idx_vec(nx, ny, DOWN_RIGHT)];
            }
        }
    }

    void calc_u() {
        for (int nx{0}; nx < NX; nx++) {
            for (int ny{0}; ny < NY; ny++) {
                double fx = Fx * dt / 2 / rho[idx_scal(nx, ny)];
                double fy = Fy * dt / 2 / rho[idx_scal(nx, ny)];
                ux[idx_scal(nx, ny)] = (f_old[idx_vec(nx, ny, RIGHT)] 
                            - f_old[idx_vec(nx, ny, LEFT)] 
                            + f_old[idx_vec(nx, ny, UP_RIGHT)] 
                            - f_old[idx_vec(nx, ny, UP_LEFT)] 
                            - f_old[idx_vec(nx, ny, DOWN_LEFT)] 
                            + f_old[idx_vec(nx, ny, DOWN_RIGHT)]) / rho[idx_scal(nx, ny)] + fx;

                uy[idx_scal(nx, ny)] = (f_old[idx_vec(nx, ny, UP)] 
                            - f_old[idx_vec(nx, ny, DOWN)] 
                            + f_old[idx_vec(nx, ny, UP_RIGHT)] 
                            + f_old[idx_vec(nx, ny, UP_LEFT)] 
                            - f_old[idx_vec(nx, ny, DOWN_LEFT)] 
                            - f_old[idx_vec(nx, ny, DOWN_RIGHT)]) / rho[idx_scal(nx, ny)] + fy;
            }
        }
    }

    void moment_update() {
        calc_rho();
        calc_u();
    }

    void collision() {
        for (int nx{0}; nx < NX; nx++) {
            for (int ny{0}; ny < NY; ny++) {
                double u2 = ux[idx_scal(nx, ny)] * ux[idx_scal(nx, ny)] + uy[idx_scal(nx, ny)] * uy[idx_scal(nx, ny)];

                double f_eq0 = (4.0/9.0) * rho[idx_scal(nx, ny)] * (1 - (3.0/2.0) * u2);

                double f_eq1 = (1.0/9.0) * rho[idx_scal(nx, ny)] * (1 + 3.0 * ux[idx_scal(nx, ny)] + (9.0/2.0) * ux[idx_scal(nx, ny)] * ux[idx_scal(nx, ny)] - (3.0/2.0) * u2);
                double f_eq2 = (1.0/9.0) * rho[idx_scal(nx, ny)] * (1 + 3.0 * uy[idx_scal(nx, ny)] + (9.0/2.0) * uy[idx_scal(nx, ny)] * uy[idx_scal(nx, ny)] - (3.0/2.0) * u2);
                double f_eq3 = (1.0/9.0) * rho[idx_scal(nx, ny)] * (1 - 3.0 * ux[idx_scal(nx, ny)] + (9.0/2.0) * ux[idx_scal(nx, ny)] * ux[idx_scal(nx, ny)] - (3.0/2.0) * u2);
                double f_eq4 = (1.0/9.0) * rho[idx_scal(nx, ny)] * (1 - 3.0 * uy[idx_scal(nx, ny)] + (9.0/2.0) * uy[idx_scal(nx, ny)] * uy[idx_scal(nx, ny)] - (3.0/2.0) * u2);

                double f_eq5 = (1.0/36.0) * rho[idx_scal(nx, ny)] * (1 + 3.0 * (ux[idx_scal(nx, ny)] + uy[idx_scal(nx, ny)]) + (9.0/2.0) * (ux[idx_scal(nx, ny)] + uy[idx_scal(nx, ny)]) * (ux[idx_scal(nx, ny)] + uy[idx_scal(nx, ny)]) - (3.0/2.0) * u2);
                double f_eq6 = (1.0/36.0) * rho[idx_scal(nx, ny)] * (1 + 3.0 * (-ux[idx_scal(nx, ny)] + uy[idx_scal(nx, ny)]) + (9.0/2.0) * (-ux[idx_scal(nx, ny)] + uy[idx_scal(nx, ny)]) * (-ux[idx_scal(nx, ny)] + uy[idx_scal(nx, ny)]) - (3.0/2.0) * u2);
                double f_eq7 = (1.0/36.0) * rho[idx_scal(nx, ny)] * (1 + 3.0 * (-ux[idx_scal(nx, ny)] - uy[idx_scal(nx, ny)]) + (9.0/2.0) * (-ux[idx_scal(nx, ny)] - uy[idx_scal(nx, ny)]) * (-ux[idx_scal(nx, ny)] - uy[idx_scal(nx, ny)]) - (3.0/2.0) * u2);
                double f_eq8 = (1.0/36.0) * rho[idx_scal(nx, ny)] * (1 + 3.0 * (ux[idx_scal(nx, ny)] - uy[idx_scal(nx, ny)]) + (9.0/2.0) * (ux[idx_scal(nx, ny)] - uy[idx_scal(nx, ny)]) * (ux[idx_scal(nx, ny)] - uy[idx_scal(nx, ny)]) - (3.0/2.0) * u2);

                double S0 = 0.0;
                double S1 = (1.0/9.0) * 3.0 * Fx;
                double S2 = (1.0/9.0) * 3.0 * Fy;
                double S3 = (1.0/9.0) * 3.0 * (-Fx);
                double S4 = (1.0/9.0) * 3.0 * (-Fy);
                double S5 = (1.0/36.0) * 3.0 * (Fx + Fy);
                double S6 = (1.0/36.0) * 3.0 * (-Fx + Fy);
                double S7 = (1.0/36.0) * 3.0 * (-Fx - Fy);
                double S8 = (1.0/36.0) * 3.0 * (Fx - Fy);

                f_old[idx_vec(nx, ny, STAY)] = f_old[idx_vec(nx, ny, STAY)] * (1 - omega) + f_eq0 * omega + S0 * (1 - 0.5*omega);
                f_old[idx_vec(nx, ny, RIGHT)] = f_old[idx_vec(nx, ny, RIGHT)] * (1 - omega) + f_eq1 * omega + S1 * (1 - 0.5*omega);
                f_old[idx_vec(nx, ny, UP)] = f_old[idx_vec(nx, ny, UP)] * (1 - omega) + f_eq2 * omega + S2 * (1 - 0.5*omega);
                f_old[idx_vec(nx, ny, LEFT)] = f_old[idx_vec(nx, ny, LEFT)] * (1 - omega) + f_eq3 * omega + S3 * (1 - 0.5*omega);
                f_old[idx_vec(nx, ny, DOWN)] = f_old[idx_vec(nx, ny, DOWN)] * (1 - omega) + f_eq4 * omega + S4 * (1 - 0.5*omega);
                f_old[idx_vec(nx, ny, UP_RIGHT)] = f_old[idx_vec(nx, ny, UP_RIGHT)] * (1 - omega) + f_eq5 * omega + S5 * (1 - 0.5*omega);
                f_old[idx_vec(nx, ny, UP_LEFT)] = f_old[idx_vec(nx, ny, UP_LEFT)] * (1 - omega) + f_eq6 * omega + S6 * (1 - 0.5*omega);
                f_old[idx_vec(nx, ny, DOWN_LEFT)] = f_old[idx_vec(nx, ny, DOWN_LEFT)] * (1 - omega) + f_eq7 * omega + S7 * (1 - 0.5*omega);
                f_old[idx_vec(nx, ny, DOWN_RIGHT)] = f_old[idx_vec(nx, ny, DOWN_RIGHT)] * (1 - omega) + f_eq8 * omega + S8 * (1 - 0.5*omega);

            }
        }
    }

    void propagation() {
        for (int nx{1}; nx < NX-1; nx++) {
            for (int ny{1}; ny < NY-1; ny++) {
                f_new[idx_vec(nx, ny, STAY)] = f_old[idx_vec(nx, ny, STAY)];
                
                f_new[idx_vec(nx, ny, RIGHT)] = f_old[idx_vec(nx-1, ny, RIGHT)];
                f_new[idx_vec(nx, ny, UP)] = f_old[idx_vec(nx, ny+1, UP)];
                f_new[idx_vec(nx, ny, LEFT)] = f_old[idx_vec(nx+1, ny, LEFT)];
                f_new[idx_vec(nx, ny, DOWN)] = f_old[idx_vec(nx, ny-1, DOWN)];
                
                f_new[idx_vec(nx, ny, UP_RIGHT)] = f_old[idx_vec(nx-1, ny+1, UP_RIGHT)];
                f_new[idx_vec(nx, ny, UP_LEFT)] = f_old[idx_vec(nx+1, ny+1, UP_LEFT)];
                f_new[idx_vec(nx, ny, DOWN_RIGHT)] = f_old[idx_vec(nx-1, ny-1, DOWN_RIGHT)];
                f_new[idx_vec(nx, ny, DOWN_LEFT)] = f_old[idx_vec(nx+1, ny-1, DOWN_LEFT)];
            }
        }

        for (int ny = 1; ny < NY - 1; ny++) {
            // periodic boundary conditions
            f_new[idx_vec(0, ny, STAY)] = f_old[idx_vec(0, ny, STAY)];
            f_new[idx_vec(0, ny, RIGHT)] = f_old[idx_vec(NX - 1, ny, RIGHT)];
            f_new[idx_vec(0, ny, UP)] = f_old[idx_vec(0, ny + 1, UP)];
            f_new[idx_vec(0, ny, LEFT)] = f_old[idx_vec(1, ny, LEFT)];
            f_new[idx_vec(0, ny, DOWN)] = f_old[idx_vec(0, ny - 1, DOWN)];
            f_new[idx_vec(0, ny, UP_RIGHT)] = f_old[idx_vec(NX - 1, ny + 1, UP_RIGHT)];
            f_new[idx_vec(0, ny, UP_LEFT)] = f_old[idx_vec(1, ny + 1, UP_LEFT)];
            f_new[idx_vec(0, ny, DOWN_LEFT)] = f_old[idx_vec(1, ny - 1, DOWN_LEFT)];
            f_new[idx_vec(0, ny, DOWN_RIGHT)] = f_old[idx_vec(NX - 1, ny - 1, DOWN_RIGHT)];

            f_new[idx_vec(NX - 1, ny, STAY)] = f_old[idx_vec(NX - 1, ny, STAY)];
            f_new[idx_vec(NX - 1, ny, RIGHT)] = f_old[idx_vec(NX - 2, ny, RIGHT)];
            f_new[idx_vec(NX - 1, ny, UP)] = f_old[idx_vec(NX - 1, ny + 1, UP)];
            f_new[idx_vec(NX - 1, ny, LEFT)] = f_old[idx_vec(0, ny, LEFT)];
            f_new[idx_vec(NX - 1, ny, DOWN)] = f_old[idx_vec(NX - 1, ny - 1, DOWN)];
            f_new[idx_vec(NX - 1, ny, UP_RIGHT)] = f_old[idx_vec(NX-2, ny + 1, UP_RIGHT)];
            f_new[idx_vec(NX - 1, ny, UP_LEFT)] = f_old[idx_vec(0, ny + 1, UP_LEFT)];
            f_new[idx_vec(NX - 1, ny, DOWN_LEFT)] = f_old[idx_vec(0, ny - 1, DOWN_LEFT)];
            f_new[idx_vec(NX - 1, ny, DOWN_RIGHT)] = f_old[idx_vec(NX-2, ny - 1, DOWN_RIGHT)];
        }

        for (int nx = 1; nx < NX - 1; nx++) {
            // no-slip conditions
            // upper wall
            f_new[idx_vec(nx, 0, STAY)] = f_old[idx_vec(nx, 0, STAY)];
            f_new[idx_vec(nx, 0, RIGHT)] = f_old[idx_vec(nx-1, 0, RIGHT)];
            f_new[idx_vec(nx, 0, UP)] = f_old[idx_vec(nx, 1, UP)];
            f_new[idx_vec(nx, 0, LEFT)] = f_old[idx_vec(nx+1, 0, LEFT)];
            f_new[idx_vec(nx, 0, DOWN)] = f_old[idx_vec(nx, 0, UP)];
            f_new[idx_vec(nx, 0, UP_RIGHT)] = f_old[idx_vec(nx-1, 1, UP_RIGHT)];
            f_new[idx_vec(nx, 0, UP_LEFT)] = f_old[idx_vec(nx+1, 1, UP_LEFT)];
            f_new[idx_vec(nx, 0, DOWN_LEFT)] = f_old[idx_vec(nx, 0, UP_RIGHT)];
            f_new[idx_vec(nx, 0, DOWN_RIGHT)] = f_old[idx_vec(nx, 0, UP_LEFT)];

            // lower wall
            f_new[idx_vec(nx, NY - 1, STAY)] = f_old[idx_vec(nx, NY-1, STAY)];
            f_new[idx_vec(nx, NY - 1, RIGHT)] = f_old[idx_vec(nx-1, NY-1, RIGHT)];
            f_new[idx_vec(nx, NY - 1, UP)] = f_old[idx_vec(nx, NY-1, DOWN)];
            f_new[idx_vec(nx, NY - 1, LEFT)] = f_old[idx_vec(nx+1, NY-1, LEFT)];
            f_new[idx_vec(nx, NY - 1, DOWN)] = f_old[idx_vec(nx, NY-2, DOWN)];
            f_new[idx_vec(nx, NY - 1, UP_RIGHT)] = f_old[idx_vec(nx, NY-1, DOWN_LEFT)];
            f_new[idx_vec(nx, NY - 1, UP_LEFT)] = f_old[idx_vec(nx, NY-1, DOWN_RIGHT)];
            f_new[idx_vec(nx, NY - 1, DOWN_LEFT)] = f_old[idx_vec(nx+1, NY-2, DOWN_LEFT)];
            f_new[idx_vec(nx, NY - 1, DOWN_RIGHT)] = f_old[idx_vec(nx-1, NY-2, DOWN_RIGHT)];
        }

        // corners
        // UP LEFT
        f_new[idx_vec(0, 0, STAY)] = f_old[idx_vec(0, 0, STAY)];
        f_new[idx_vec(0, 0, RIGHT)] = f_old[idx_vec(NX-1, 0, RIGHT)];
        f_new[idx_vec(0, 0, UP)] = f_old[idx_vec(0, 1, UP)];
        f_new[idx_vec(0, 0, LEFT)] = f_old[idx_vec(1, 0, LEFT)];
        f_new[idx_vec(0, 0, DOWN)] = f_old[idx_vec(0, 0, UP)];
        f_new[idx_vec(0, 0, UP_RIGHT)] = f_old[idx_vec(NX-1, 1, UP_RIGHT)];
        f_new[idx_vec(0, 0, UP_LEFT)] = f_old[idx_vec(1, 1, UP_LEFT)];
        f_new[idx_vec(0, 0, DOWN_LEFT)] = f_old[idx_vec(0, 0, UP_RIGHT)];
        f_new[idx_vec(0, 0, DOWN_RIGHT)] = f_old[idx_vec(0, 0, UP_LEFT)];

        // UP RIGHT
        f_new[idx_vec(NX-1, 0, STAY)] = f_old[idx_vec(NX-1, 0, STAY)];
        f_new[idx_vec(NX-1, 0, RIGHT)] = f_old[idx_vec(NX-2, 0, RIGHT)];
        f_new[idx_vec(NX-1, 0, UP)] = f_old[idx_vec(NX-1, 1, UP)];
        f_new[idx_vec(NX-1, 0, LEFT)] = f_old[idx_vec(0, 0, LEFT)];
        f_new[idx_vec(NX-1, 0, DOWN)] = f_old[idx_vec(NX-1, 0, UP)];
        f_new[idx_vec(NX-1, 0, UP_RIGHT)] = f_old[idx_vec(NX-2, 1, UP_RIGHT)];
        f_new[idx_vec(NX-1, 0, UP_LEFT)] = f_old[idx_vec(0, 1, UP_LEFT)];
        f_new[idx_vec(NX-1, 0, DOWN_LEFT)] = f_old[idx_vec(NX-1, 0, UP_RIGHT)];
        f_new[idx_vec(NX-1, 0, DOWN_RIGHT)] = f_old[idx_vec(NX-1, 0, UP_LEFT)];

        // DOWN LEFT
        f_new[idx_vec(0, NY-1, STAY)] = f_old[idx_vec(0, NY-1, STAY)];
        f_new[idx_vec(0, NY-1, RIGHT)] = f_old[idx_vec(NX-1, NY-1, RIGHT)];
        f_new[idx_vec(0, NY-1, UP)] = f_old[idx_vec(0, NY-1, DOWN)];
        f_new[idx_vec(0, NY-1, LEFT)] = f_old[idx_vec(1, NY-1, LEFT)];
        f_new[idx_vec(0, NY-1, DOWN)] = f_old[idx_vec(0, NY-2, DOWN)];
        f_new[idx_vec(0, NY-1, UP_RIGHT)] = f_old[idx_vec(0, NY-1, DOWN_LEFT)];
        f_new[idx_vec(0, NY-1, UP_LEFT)] = f_old[idx_vec(0, NY-1, DOWN_RIGHT)];
        f_new[idx_vec(0, NY-1, DOWN_LEFT)] = f_old[idx_vec(1, NY-2, DOWN_LEFT)];
        f_new[idx_vec(0, NY-1, DOWN_RIGHT)] = f_old[idx_vec(NX-1, NY-2, DOWN_RIGHT)];

        // DOWN RIGHT
        f_new[idx_vec(NX-1, NY-1, STAY)] = f_old[idx_vec(NX-1, NY-1, STAY)];
        f_new[idx_vec(NX-1, NY-1, RIGHT)] = f_old[idx_vec(NX-2, NY-1, RIGHT)];
        f_new[idx_vec(NX-1, NY-1, UP)] = f_old[idx_vec(NX-1, NY-1, DOWN)];
        f_new[idx_vec(NX-1, NY-1, LEFT)] = f_old[idx_vec(0, NY-1, LEFT)];
        f_new[idx_vec(NX-1, NY-1, DOWN)] = f_old[idx_vec(NX-1, NY-2, DOWN)];
        f_new[idx_vec(NX-1, NY-1, UP_RIGHT)] = f_old[idx_vec(NX-1, NY-1, DOWN_LEFT)];
        f_new[idx_vec(NX-1, NY-1, UP_LEFT)] = f_old[idx_vec(NX-1, NY-1, DOWN_RIGHT)];
        f_new[idx_vec(NX-1, NY-1, DOWN_LEFT)] = f_old[idx_vec(0, NY-2, DOWN_LEFT)];
        f_new[idx_vec(NX-1, NY-1, DOWN_RIGHT)] = f_old[idx_vec(NX-2, NY-2, DOWN_RIGHT)];

        swap(f_new, f_old);
    }

public:

    void do_step() {
        moment_update();
        collision();
        propagation();
        time += dt;
    }

    void save(const char* file_path) {
        std::ofstream file(file_path);

        for (int ny{}; ny < NY; ny++) {
            for (int nx{}; nx < NX-1; nx++) {
                file << ux[idx_scal(nx, ny)] << " ";
            }
            file << ux[idx_scal(NX-1, ny)] << "\n";
        }
        file.close();
    }
};


int main(int argc, char* argv[]) {
    int NX = std::stoi(argv[1]);
    int NY = std::stoi(argv[2]);
    double FX = std::stod(argv[3]);
    double tau = std::stod(argv[4]);
    int steps = std::stoi(argv[5]);

    Solver solver(NX, NY, FX, tau);
    for (int i{}; i < steps; i++) {
        solver.do_step();
    }
    solver.save("result.csv");
}