#include "model.hpp"
#include <fstream>

/*
    D2Q5:
    Directions numeration:
       2  
    3  0  1
       4  
*/
#define Qg 5
#define W0g (1.0 / 3.0)
#define W14g (1.0 / 6.0)
#define STAYg 0
#define RIGHTg 1
#define UPg 2
#define LEFTg 3
#define DOWNg 4

/*
    D2Q4:
    Directions numeration:
       1  
    2     0
       3  
*/
#define Qh 4
#define Wh (1.0 / 4.0)
#define RIGHTh 0
#define UPh 1
#define LEFTh 2
#define DOWNh 3



Solver::Solver(int nx, int ny, double tau_h, double tau_g, double p_A, double p_kappa, double p_D)
    : NX(nx), NY(ny), A(p_A), kappa(p_kappa), D(p_D) {
    dt = 1.0;
    omega_h = dt / tau_h;
    omega_g = dt / tau_g;
    time = 0;
    h_old = new double[NX*NY*Qh];
    h_new = new double[NX*NY*Qh];
    g_old = new double[NX*NY*Qg];
    g_new = new double[NX*NY*Qg];
    phi = new double[NX*NY];
    C = new double[NX*NY];

    init();
}

Solver::~Solver() {
    delete[] h_old;
    delete[] h_new;
    delete[] g_old;
    delete[] g_new;
    delete[] phi;
    delete[] C;
}


inline int Solver::idx_vec_g(int x, int y, int n) const {
    return (x * NY + y) * Qg + n;
}

inline int Solver::idx_vec_h(int x, int y, int n) const {
    return (x * NY + y) * Qh + n;
}

inline int Solver::idx_scal(int x, int y) const {
    return x * NY + y;
}

void Solver::do_step() {
    moment_update();
    collision();
    propagation();
    time += dt;
}

void Solver::save(const char* file_path)
{
    std::ofstream file(file_path);

    for (int ny{}; ny < NY; ny++) {
        for (int nx{}; nx < NX-1; nx++) {
            file << phi[idx_scal(nx, ny)] << " ";
        }
        file << phi[idx_scal(NX-1, ny)] << "\n";
    }
    file.close();
}

inline void Solver::init()
{
    for (int nx{0}; nx < NX; nx++) {
        for (int ny{0}; ny < NY; ny++) {
            h_old[idx_vec_h(nx, ny, RIGHTh)] = Wh;
            h_old[idx_vec_h(nx, ny, UPh)] = Wh;
            h_old[idx_vec_h(nx, ny, LEFTh)] = Wh;
            h_old[idx_vec_h(nx, ny, DOWNh)] = Wh;

            g_old[idx_vec_g(nx, ny, STAYg)] = W0g;
            g_old[idx_vec_g(nx, ny, RIGHTg)] = W14g;
            g_old[idx_vec_g(nx, ny, UPg)] = W14g;
            g_old[idx_vec_g(nx, ny, LEFTg)] = W14g;
            g_old[idx_vec_g(nx, ny, DOWNg)] = W14g;
        }
    }
}

inline void Solver::moment_update()
{
    calc_phi();
    calc_C();
}

inline void Solver::calc_phi() {
    for (int nx{0}; nx < NX; nx++) {
        for (int ny{0}; ny < NY; ny++) {
            phi[idx_scal(nx, ny)] = (
                h_old[idx_vec_h(nx, ny, RIGHTh)] +
                h_old[idx_vec_h(nx, ny, UPh)] +
                h_old[idx_vec_h(nx, ny, LEFTh)] +
                h_old[idx_vec_h(nx, ny, DOWNh)]
            );
        }
    }
}

inline void Solver::calc_C() {
    for (int nx{0}; nx < NX; nx++) {
        for (int ny{0}; ny < NY; ny++) {
            C[idx_scal(nx, ny)] = (
                g_old[idx_vec_g(nx, ny, STAYg)] +
                g_old[idx_vec_g(nx, ny, RIGHTg)] +
                g_old[idx_vec_g(nx, ny, UPg)] +
                g_old[idx_vec_g(nx, ny, LEFTg)] +
                g_old[idx_vec_g(nx, ny, DOWNg)]
            );
        }
    }
}

inline void Solver::collision() {
    for (int nx{0}; nx < NX; nx++) {
        for (int ny{0}; ny < NY; ny++) {
            double h_eq14 = W14 * phi[idx_scal(nx, ny)];

            h_old[idx_vec(nx, ny, RIGHT)] = h_old[idx_vec(nx, ny, RIGHT)] * (1 - omega_h) + h_eq14 * omega_h;
            h_old[idx_vec(nx, ny, UP)] = h_old[idx_vec(nx, ny, UP)] * (1 - omega_h) + h_eq14 * omega_h;
            h_old[idx_vec(nx, ny, LEFT)] = h_old[idx_vec(nx, ny, LEFT)] * (1 - omega_h) + h_eq14 * omega_h;
            h_old[idx_vec(nx, ny, DOWN)] = h_old[idx_vec(nx, ny, DOWN)] * (1 - omega_h) + h_eq14 * omega_h;
        }
    }
}