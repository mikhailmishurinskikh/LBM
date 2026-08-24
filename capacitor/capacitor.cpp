#include <fstream>

#define STAY 0
#define RIGHT 1
#define UP 2
#define LEFT 3
#define DOWN 4

#define W0 1.0/3.0
#define W14 1.0/6.0

inline void swap(double*& ptr1, double*& ptr2) {
    double* tmp = ptr1;
    ptr1 = ptr2;
    ptr2 = tmp;
}

class Solver {
    /*
    D2Q5:
    Directions numeration:
       2  
    3  0  1
       4  
    
    Capacitor: left - anode, right - catode
    */
    int NX;
    int NY;

    // Electric potential
    double* h_new;
    double* h_old;

    double* C;
    double* phi;

    double dt;

    // omega = dt / tau
    double omega_h;
    double time;

    double left_phi{1.0};
    double right_phi{0.0};

public:

    Solver(int nx, int ny, double tau_h) : NX(nx), NY(ny)
    {
        dt = 1.0;
        omega_h = dt / tau_h;
        time = 0;
        h_old = new double[NX*NY*5];
        h_new = new double[NX*NY*5];
        phi = new double[NX*NY];

        init();
    }

    ~Solver() {
        delete[] h_old;
        delete[] h_new;
        delete[] C;
        delete[] phi;
    }

private:

    inline int idx_vec(int x, int y, int n) const {
        return (x * NY + y) * 5 + n;
    }

    inline int idx_scal(int x, int y) const {
        return x * NY + y;
    }    

    void init() {
        for (int nx{0}; nx < NX; nx++) {
            for (int ny{0}; ny < NY; ny++) {
                h_old[idx_vec(nx, ny, STAY)] = W0;
                h_old[idx_vec(nx, ny, RIGHT)] = W14;
                h_old[idx_vec(nx, ny, UP)] = W14;
                h_old[idx_vec(nx, ny, LEFT)] = W14;
                h_old[idx_vec(nx, ny, DOWN)] = W14;
            }
        }
    }

    void calc_phi() {
        for (int nx{0}; nx < NX; nx++) {
            for (int ny{0}; ny < NY; ny++) {
                phi[idx_scal(nx, ny)] = (
                    h_old[idx_vec(nx, ny, STAY)] +
                    h_old[idx_vec(nx, ny, RIGHT)] +
                    h_old[idx_vec(nx, ny, UP)] +
                    h_old[idx_vec(nx, ny, LEFT)] +
                    h_old[idx_vec(nx, ny, DOWN)]
                );
            }
        }
    }

    void moment_update() {
        calc_phi();
    }

    void collision() {
        for (int nx{0}; nx < NX; nx++) {
            for (int ny{0}; ny < NY; ny++) {
                double h_eq0 = W0 * phi[idx_scal(nx, ny)];
                double h_eq14 = W14 * phi[idx_scal(nx, ny)];

                h_old[idx_vec(nx, ny, STAY)] = h_old[idx_vec(nx, ny, STAY)] * (1 - omega_h) + h_eq0 * omega_h;
                h_old[idx_vec(nx, ny, RIGHT)] = h_old[idx_vec(nx, ny, RIGHT)] * (1 - omega_h) + h_eq14 * omega_h;
                h_old[idx_vec(nx, ny, UP)] = h_old[idx_vec(nx, ny, UP)] * (1 - omega_h) + h_eq14 * omega_h;
                h_old[idx_vec(nx, ny, LEFT)] = h_old[idx_vec(nx, ny, LEFT)] * (1 - omega_h) + h_eq14 * omega_h;
                h_old[idx_vec(nx, ny, DOWN)] = h_old[idx_vec(nx, ny, DOWN)] * (1 - omega_h) + h_eq14 * omega_h;
            }
        }
    }

    void propagation() {
        for (int nx{1}; nx < NX-1; nx++) {
            for (int ny{1}; ny < NY-1; ny++) {
                h_new[idx_vec(nx, ny, STAY)] = h_old[idx_vec(nx, ny, STAY)];
                
                h_new[idx_vec(nx, ny, RIGHT)] = h_old[idx_vec(nx-1, ny, RIGHT)];
                h_new[idx_vec(nx, ny, UP)] = h_old[idx_vec(nx, ny+1, UP)];
                h_new[idx_vec(nx, ny, LEFT)] = h_old[idx_vec(nx+1, ny, LEFT)];
                h_new[idx_vec(nx, ny, DOWN)] = h_old[idx_vec(nx, ny-1, DOWN)];
            }
        }

        for (int ny = 1; ny < NY - 1; ny++) {
            // Dirichlet boundary conditions
            // left electrode (anode)
            h_new[idx_vec(0, ny, STAY)] = h_old[idx_vec(0, ny, STAY)];
            h_new[idx_vec(0, ny, RIGHT)] = -h_old[idx_vec(0, ny, LEFT)] + 2 * W14 * left_phi;
            h_new[idx_vec(0, ny, UP)] = h_old[idx_vec(0, ny + 1, UP)];
            h_new[idx_vec(0, ny, LEFT)] = h_old[idx_vec(1, ny, LEFT)];
            h_new[idx_vec(0, ny, DOWN)] = h_old[idx_vec(0, ny - 1, DOWN)];

            // right electrode (cathode)
            h_new[idx_vec(NX - 1, ny, STAY)] = h_old[idx_vec(NX - 1, ny, STAY)];
            h_new[idx_vec(NX - 1, ny, RIGHT)] = h_old[idx_vec(NX - 2, ny, RIGHT)];
            h_new[idx_vec(NX - 1, ny, UP)] = h_old[idx_vec(NX - 1, ny + 1, UP)];
            h_new[idx_vec(NX - 1, ny, LEFT)] = -h_old[idx_vec(NX - 1, ny, RIGHT)] + 2 * W14 * right_phi;
            h_new[idx_vec(NX - 1, ny, DOWN)] = h_old[idx_vec(NX - 1, ny - 1, DOWN)];
        }

        for (int nx = 1; nx < NX - 1; nx++) {
            // bounce-back conditions
            // upper wall
            h_new[idx_vec(nx, 0, STAY)] = h_old[idx_vec(nx, 0, STAY)];
            h_new[idx_vec(nx, 0, RIGHT)] = h_old[idx_vec(nx-1, 0, RIGHT)];
            h_new[idx_vec(nx, 0, UP)] = h_old[idx_vec(nx, 1, UP)];
            h_new[idx_vec(nx, 0, LEFT)] = h_old[idx_vec(nx+1, 0, LEFT)];
            h_new[idx_vec(nx, 0, DOWN)] = h_old[idx_vec(nx, 0, UP)];

            // lower wall
            h_new[idx_vec(nx, NY - 1, STAY)] = h_old[idx_vec(nx, NY-1, STAY)];
            h_new[idx_vec(nx, NY - 1, RIGHT)] = h_old[idx_vec(nx-1, NY-1, RIGHT)];
            h_new[idx_vec(nx, NY - 1, UP)] = h_old[idx_vec(nx, NY-1, DOWN)];
            h_new[idx_vec(nx, NY - 1, LEFT)] = h_old[idx_vec(nx+1, NY-1, LEFT)];
            h_new[idx_vec(nx, NY - 1, DOWN)] = h_old[idx_vec(nx, NY-2, DOWN)];
        }

        // corners
        // UP LEFT
        h_new[idx_vec(0, 0, STAY)] = h_old[idx_vec(0, 0, STAY)];
        h_new[idx_vec(0, 0, RIGHT)] = -h_old[idx_vec(0, 0, LEFT)] + 2 * W14 * left_phi;
        h_new[idx_vec(0, 0, UP)] = h_old[idx_vec(0, 1, UP)];
        h_new[idx_vec(0, 0, LEFT)] = h_old[idx_vec(1, 0, LEFT)];
        h_new[idx_vec(0, 0, DOWN)] = h_old[idx_vec(0, 0, UP)];

        // UP RIGHT
        h_new[idx_vec(NX-1, 0, STAY)] = h_old[idx_vec(NX-1, 0, STAY)];
        h_new[idx_vec(NX-1, 0, RIGHT)] = h_old[idx_vec(NX-2, 0, RIGHT)];
        h_new[idx_vec(NX-1, 0, UP)] = h_old[idx_vec(NX-1, 1, UP)];
        h_new[idx_vec(NX-1, 0, LEFT)] = -h_old[idx_vec(NX-1, 0, RIGHT)] + 2 * W14 * right_phi;
        h_new[idx_vec(NX-1, 0, DOWN)] = h_old[idx_vec(NX-1, 0, UP)];

        // DOWN LEFT
        h_new[idx_vec(0, NY-1, STAY)] = h_old[idx_vec(0, NY-1, STAY)];
        h_new[idx_vec(0, NY-1, RIGHT)] = -h_old[idx_vec(0, NY-1, LEFT)] + 2 * W14 * left_phi;
        h_new[idx_vec(0, NY-1, UP)] = h_old[idx_vec(0, NY-1, DOWN)];
        h_new[idx_vec(0, NY-1, LEFT)] = h_old[idx_vec(1, NY-1, LEFT)];
        h_new[idx_vec(0, NY-1, DOWN)] = h_old[idx_vec(0, NY-2, DOWN)];

        // DOWN RIGHT
        h_new[idx_vec(NX-1, NY-1, STAY)] = h_old[idx_vec(NX-1, NY-1, STAY)];
        h_new[idx_vec(NX-1, NY-1, RIGHT)] = h_old[idx_vec(NX-2, NY-1, RIGHT)];
        h_new[idx_vec(NX-1, NY-1, UP)] = h_old[idx_vec(NX-1, NY-1, DOWN)];
        h_new[idx_vec(NX-1, NY-1, LEFT)] = -h_old[idx_vec(NX-1, NY-1, RIGHT)] + 2 * W14 * right_phi;
        h_new[idx_vec(NX-1, NY-1, DOWN)] = h_old[idx_vec(NX-1, NY-2, DOWN)];

        swap(h_new, h_old);
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
                file << phi[idx_scal(nx, ny)] << " ";
            }
            file << phi[idx_scal(NX-1, ny)] << "\n";
        }
        file.close();
    }
};


int main(int argc, char* argv[]) {
    int NX = 50;
    int NY = 30;
    double tau = 1;
    int steps = 10000;

    Solver solver(NX, NY, tau);
    for (int i{}; i < steps; i++) {
        solver.do_step();
    }
    solver.save("result.csv");
}