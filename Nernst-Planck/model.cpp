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



Solver::Solver(int nx, int ny, double tau_h, double tau_g, double F, double kappa, double Deff, double RT)
    : NX(nx), NY(ny), F(F), kappa(kappa), Deff(Deff), RT(RT) {
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

void Solver::do_step_g() {
    calc_phi();
    collision_g();
    propagation_g();
    time += dt;
}

void Solver::do_step_h() {
    calc_C();
    collision_h();
    propagation_h();
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
    for (int nx{0}; nx < NX - 1; nx++) {
        for (int ny{0}; ny < NY - 1; ny++) {
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

inline void Solver::collision_h() {
    // internal dots
    for (int nx{1}; nx < NX - 1; nx++) {
        for (int ny{1}; ny < NY - 1; ny++) {
            double laplace_C = (
                C[idx_scal(nx, ny+1)] + C[idx_scal(nx+1, ny)]
                + C[idx_scal(nx, ny-1)] + C[idx_scal(nx-1, ny)]
                - 4 * C[idx_scal(nx, ny)]
            );
            double Q = Wh * F * Deff / (2 * kappa) * (0.5 * dt - tau_h) * laplace_C * dt;

            double h_eq = Wh * phi[idx_scal(nx, ny)];
            
            h_old[idx_vec_h(nx, ny, RIGHTh)]  = (1.0 - omega_h) * h_old[idx_vec_h(nx, ny, RIGHTh)]  + omega_h * h_eq + Q;
            h_old[idx_vec_h(nx, ny, UPh)]     = (1.0 - omega_h) * h_old[idx_vec_h(nx, ny, UPh)]     + omega_h * h_eq + Q;
            h_old[idx_vec_h(nx, ny, LEFTh)]   = (1.0 - omega_h) * h_old[idx_vec_h(nx, ny, LEFTh)]   + omega_h * h_eq + Q;
            h_old[idx_vec_h(nx, ny, DOWNh)]   = (1.0 - omega_h) * h_old[idx_vec_h(nx, ny, DOWNh)]   + omega_h * h_eq + Q;
        }
    }

    // upper wall
    for (int nx{1}; nx < NX - 1; nx++) {
        double laplace_C = (
            2.0 * C[idx_scal(nx, 1)] + C[idx_scal(nx+1, 0)]
            + C[idx_scal(nx-1, 0)]
            - 4 * C[idx_scal(nx, 0)]
        );
        double Q = Wh * F * Deff / (2 * kappa) * (0.5 * dt - tau_h) * laplace_C * dt;

        double h_eq = Wh * phi[idx_scal(nx, 0)];
        h_old[idx_vec_h(nx, 0, RIGHTh)]  = (1.0 - omega_h) * h_old[idx_vec_h(nx, 0, RIGHTh)]  + omega_h * h_eq + Q;
        h_old[idx_vec_h(nx, 0, UPh)]     = (1.0 - omega_h) * h_old[idx_vec_h(nx, 0, UPh)]     + omega_h * h_eq + Q;
        h_old[idx_vec_h(nx, 0, LEFTh)]   = (1.0 - omega_h) * h_old[idx_vec_h(nx, 0, LEFTh)]   + omega_h * h_eq + Q;
        h_old[idx_vec_h(nx, 0, DOWNh)]   = (1.0 - omega_h) * h_old[idx_vec_h(nx, 0, DOWNh)]   + omega_h * h_eq + Q;
    }

    // bottom wall
    for (int nx{1}; nx < NX - 1; nx++) {
        double laplace_C = (
            2.0 * C[idx_scal(nx, NY-2)] + C[idx_scal(nx+1, NY-1)]
            + C[idx_scal(nx-1, NY-1)]
            - 4.0 * C[idx_scal(nx, NY-1)]
        );
        
        double Q = Wh * F * Deff / (2.0 * kappa) * (0.5 * dt - tau_h) * laplace_C * dt;
        double h_eq = Wh * phi[idx_scal(nx, NY-1)];
        
        h_old[idx_vec_h(nx, NY-1, RIGHTh)] = (1.0 - omega_h) * h_old[idx_vec_h(nx, NY-1, RIGHTh)]  + omega_h * h_eq + Q;
        h_old[idx_vec_h(nx, NY-1, UPh)]    = (1.0 - omega_h) * h_old[idx_vec_h(nx, NY-1, UPh)]     + omega_h * h_eq + Q;
        h_old[idx_vec_h(nx, NY-1, LEFTh)]  = (1.0 - omega_h) * h_old[idx_vec_h(nx, NY-1, LEFTh)]   + omega_h * h_eq + Q;
        h_old[idx_vec_h(nx, NY-1, DOWNh)]  = (1.0 - omega_h) * h_old[idx_vec_h(nx, NY-1, DOWNh)]   + omega_h * h_eq + Q;
    }

    // left wall
    for (int ny{1}; ny < NY - 1; ny++) {
        double laplace_C = (
            2.0 * C[idx_scal(1, ny)] + C[idx_scal(0, ny+1)]
            + C[idx_scal(0, ny-1)]
            - 4.0 * C[idx_scal(0, ny)]
        );
        
        double Q = Wh * F * Deff / (2.0 * kappa) * (0.5 * dt - tau_h) * laplace_C * dt;
        double h_eq = Wh * phi[idx_scal(0, ny)];
        
        h_old[idx_vec_h(0, ny, RIGHTh)]  = (1.0 - omega_h) * h_old[idx_vec_h(0, ny, RIGHTh)]  + omega_h * h_eq + Q;
        h_old[idx_vec_h(0, ny, UPh)]     = (1.0 - omega_h) * h_old[idx_vec_h(0, ny, UPh)]     + omega_h * h_eq + Q;
        h_old[idx_vec_h(0, ny, LEFTh)]   = (1.0 - omega_h) * h_old[idx_vec_h(0, ny, LEFTh)]   + omega_h * h_eq + Q;
        h_old[idx_vec_h(0, ny, DOWNh)]   = (1.0 - omega_h) * h_old[idx_vec_h(0, ny, DOWNh)]   + omega_h * h_eq + Q;
    }

    // right wall
    for (int ny{1}; ny < NY - 1; ny++) {
        double laplace_C = (
            2.0 * C[idx_scal(NX-2, ny)] + C[idx_scal(NX-1, ny+1)]
            + C[idx_scal(NX-1, ny-1)]
            - 4.0 * C[idx_scal(NX-1, ny)]
        );
        
        double Q = Wh * F * Deff / (2.0 * kappa) * (0.5 * dt - tau_h) * laplace_C * dt;
        double h_eq = Wh * phi[idx_scal(NX-1, ny)];
        
        h_old[idx_vec_h(NX-1, ny, RIGHTh)]  = (1.0 - omega_h) * h_old[idx_vec_h(NX-1, ny, RIGHTh)]  + omega_h * h_eq + Q;
        h_old[idx_vec_h(NX-1, ny, UPh)]     = (1.0 - omega_h) * h_old[idx_vec_h(NX-1, ny, UPh)]     + omega_h * h_eq + Q;
        h_old[idx_vec_h(NX-1, ny, LEFTh)]   = (1.0 - omega_h) * h_old[idx_vec_h(NX-1, ny, LEFTh)]   + omega_h * h_eq + Q;
        h_old[idx_vec_h(NX-1, ny, DOWNh)]   = (1.0 - omega_h) * h_old[idx_vec_h(NX-1, ny, DOWNh)]   + omega_h * h_eq + Q;
    }

    // corners
    // UP LEFT
    {
        double laplace_C = (
            2.0 * C[idx_scal(1, 0)] + 2.0 * C[idx_scal(0, 1)]
            - 4.0 * C[idx_scal(0, 0)]
        );
        
        double Q = Wh * F * Deff / (2.0 * kappa) * (0.5 * dt - tau_h) * laplace_C * dt;
        double h_eq = Wh * phi[idx_scal(0, 0)];
        
        h_old[idx_vec_h(0, 0, RIGHTh)]  = (1.0 - omega_h) * h_old[idx_vec_h(0, 0, RIGHTh)]  + omega_h * h_eq + Q;
        h_old[idx_vec_h(0, 0, UPh)]     = (1.0 - omega_h) * h_old[idx_vec_h(0, 0, UPh)]     + omega_h * h_eq + Q;
        h_old[idx_vec_h(0, 0, LEFTh)]   = (1.0 - omega_h) * h_old[idx_vec_h(0, 0, LEFTh)]   + omega_h * h_eq + Q;
        h_old[idx_vec_h(0, 0, DOWNh)]   = (1.0 - omega_h) * h_old[idx_vec_h(0, 0, DOWNh)]   + omega_h * h_eq + Q;
    }

    // UP RIGHT
    {
        double laplace_C = (
            2.0 * C[idx_scal(NX-2, 0)] + 2.0 * C[idx_scal(NX-1, 1)]
            - 4.0 * C[idx_scal(NX-1, 0)]
        );
        
        double Q = Wh * F * Deff / (2.0 * kappa) * (0.5 * dt - tau_h) * laplace_C * dt;
        double h_eq = Wh * phi[idx_scal(NX-1, 0)];
        
        h_old[idx_vec_h(NX-1, 0, RIGHTh)]  = (1.0 - omega_h) * h_old[idx_vec_h(NX-1, 0, RIGHTh)]  + omega_h * h_eq + Q;
        h_old[idx_vec_h(NX-1, 0, UPh)]     = (1.0 - omega_h) * h_old[idx_vec_h(NX-1, 0, UPh)]     + omega_h * h_eq + Q;
        h_old[idx_vec_h(NX-1, 0, LEFTh)]   = (1.0 - omega_h) * h_old[idx_vec_h(NX-1, 0, LEFTh)]   + omega_h * h_eq + Q;
        h_old[idx_vec_h(NX-1, 0, DOWNh)]   = (1.0 - omega_h) * h_old[idx_vec_h(NX-1, 0, DOWNh)]   + omega_h * h_eq + Q;
    }

    // DOWN LEFT
    {
        double laplace_C = (
            2.0 * C[idx_scal(1, NY-1)] + 2.0 * C[idx_scal(0, NY-2)]
            - 4.0 * C[idx_scal(0, NY-1)]
        );
        
        double Q = Wh * F * Deff / (2.0 * kappa) * (0.5 * dt - tau_h) * laplace_C * dt;
        double h_eq = Wh * phi[idx_scal(0, NY-1)];
        
        h_old[idx_vec_h(0, NY-1, RIGHTh)]  = (1.0 - omega_h) * h_old[idx_vec_h(0, NY-1, RIGHTh)]  + omega_h * h_eq + Q;
        h_old[idx_vec_h(0, NY-1, UPh)]     = (1.0 - omega_h) * h_old[idx_vec_h(0, NY-1, UPh)]     + omega_h * h_eq + Q;
        h_old[idx_vec_h(0, NY-1, LEFTh)]   = (1.0 - omega_h) * h_old[idx_vec_h(0, NY-1, LEFTh)]   + omega_h * h_eq + Q;
        h_old[idx_vec_h(0, NY-1, DOWNh)]   = (1.0 - omega_h) * h_old[idx_vec_h(0, NY-1, DOWNh)]   + omega_h * h_eq + Q;
    }

    // DOWN RIGHT
    {
        double laplace_C = (
            2.0 * C[idx_scal(NX-2, NY-1)] + 2.0 * C[idx_scal(NX-1, NY-2)]
            - 4.0 * C[idx_scal(NX-1, NY-1)]
        );
        
        double Q = Wh * F * Deff / (2.0 * kappa) * (0.5 * dt - tau_h) * laplace_C * dt;
        double h_eq = Wh * phi[idx_scal(NX-1, NY-1)];
        
        h_old[idx_vec_h(NX-1, NY-1, RIGHTh)]  = (1.0 - omega_h) * h_old[idx_vec_h(NX-1, NY-1, RIGHTh)]  + omega_h * h_eq + Q;
        h_old[idx_vec_h(NX-1, NY-1, UPh)]     = (1.0 - omega_h) * h_old[idx_vec_h(NX-1, NY-1, UPh)]     + omega_h * h_eq + Q;
        h_old[idx_vec_h(NX-1, NY-1, LEFTh)]   = (1.0 - omega_h) * h_old[idx_vec_h(NX-1, NY-1, LEFTh)]   + omega_h * h_eq + Q;
        h_old[idx_vec_h(NX-1, NY-1, DOWNh)]   = (1.0 - omega_h) * h_old[idx_vec_h(NX-1, NY-1, DOWNh)]   + omega_h * h_eq + Q;
    }
}

inline void Solver::propagation_h() {
    for (int nx{1}; nx < NX-1; nx++) {
        for (int ny{1}; ny < NY-1; ny++) {            
            h_new[idx_vec_h(nx, ny, RIGHTh)] = h_old[idx_vec_h(nx-1, ny, RIGHTh)];
            h_new[idx_vec_h(nx, ny, UPh)] = h_old[idx_vec_h(nx, ny+1, UPh)];
            h_new[idx_vec_h(nx, ny, LEFTh)] = h_old[idx_vec_h(nx+1, ny, LEFTh)];
            h_new[idx_vec_h(nx, ny, DOWNh)] = h_old[idx_vec_h(nx, ny-1, DOWNh)];
        }
    }

    double known_sum_h, inamuro_phi;
    for (int ny = 1; ny < NY - 1; ny++) {
        // Dirichlet boundary conditions
        // left electrode (anode)
        h_new[idx_vec_h(0, ny, UPh)] = h_old[idx_vec_h(0, ny + 1, UPh)];
        h_new[idx_vec_h(0, ny, LEFTh)] = h_old[idx_vec_h(1, ny, LEFTh)];
        h_new[idx_vec_h(0, ny, DOWNh)] = h_old[idx_vec_h(0, ny - 1, DOWNh)];

        known_sum_h = h_new[idx_vec_h(0, ny, UPh)] 
                + h_new[idx_vec_h(0, ny, LEFTh)] 
                + h_new[idx_vec_h(0, ny, DOWNh)];

        inamuro_phi = (left_phi - known_sum_h) / Wh;
        h_new[idx_vec_h(0, ny, RIGHTh)] = Wh * inamuro_phi;

        // right electrode (cathode)
        h_new[idx_vec_h(NX - 1, ny, RIGHTh)] = h_old[idx_vec_h(NX - 2, ny, RIGHTh)];
        h_new[idx_vec_h(NX - 1, ny, UPh)] = h_old[idx_vec_h(NX - 1, ny + 1, UPh)];
        h_new[idx_vec_h(NX - 1, ny, DOWNh)] = h_old[idx_vec_h(NX - 1, ny - 1, DOWNh)];

        known_sum_h = h_new[idx_vec_h(NX-1, ny, UPh)] 
                + h_new[idx_vec_h(NX-1, ny, RIGHTh)] 
                + h_new[idx_vec_h(NX-1, ny, DOWNh)];
        inamuro_phi = (right_phi - known_sum_h) / Wh;
        h_new[idx_vec_h(NX - 1, ny, LEFTh)] = Wh * inamuro_phi;
    }

    for (int nx = 1; nx < NX - 1; nx++) {
        // bounce-back conditions (wet-node)
        // upper wall
        h_new[idx_vec_h(nx, 0, RIGHTh)] = h_old[idx_vec_h(nx-1, 0, RIGHTh)];
        h_new[idx_vec_h(nx, 0, UPh)] = 0;
        h_new[idx_vec_h(nx, 0, LEFTh)] = h_old[idx_vec_h(nx+1, 0, LEFTh)];
        h_new[idx_vec_h(nx, 0, DOWNh)] = h_old[idx_vec_h(nx, 1, UPh)];

        // lower wall
        h_new[idx_vec_h(nx, NY - 1, RIGHTh)] = h_old[idx_vec_h(nx-1, NY-1, RIGHTh)];
        h_new[idx_vec_h(nx, NY - 1, DOWNh)] = 0;
        h_new[idx_vec_h(nx, NY - 1, LEFTh)] = h_old[idx_vec_h(nx+1, NY-1, LEFTh)];
        h_new[idx_vec_h(nx, NY - 1, UPh)] = h_old[idx_vec_h(nx, NY-2, DOWNh)];
    }

    // corners
    // UP LEFT
    h_new[idx_vec_h(0, 0, UPh)] = h_old[idx_vec_h(0, 1, UPh)];
    h_new[idx_vec_h(0, 0, LEFTh)] = h_old[idx_vec_h(1, 0, LEFTh)];
    h_new[idx_vec_h(0, 0, DOWNh)] = h_old[idx_vec_h(0, 1, UPh)];

    known_sum_h = h_new[idx_vec_h(0, 0, UPh)] 
                + h_new[idx_vec_h(0, 0, LEFTh)]
                + h_new[idx_vec_h(0, 0, DOWNh)];
    inamuro_phi = (left_phi - known_sum_h) / Wh;
    h_new[idx_vec_h(0, 0, RIGHTh)] = Wh * inamuro_phi;

    // UP RIGHT
    h_new[idx_vec_h(NX-1, 0, RIGHTh)] = h_old[idx_vec_h(NX-2, 0, RIGHTh)];
    h_new[idx_vec_h(NX-1, 0, UPh)] = h_old[idx_vec_h(NX-1, 1, UPh)];
    h_new[idx_vec_h(NX-1, 0, DOWNh)] = h_old[idx_vec_h(NX-1, 1, UPh)];

    known_sum_h = h_new[idx_vec_h(NX-1, 0, UPh)] 
                + h_new[idx_vec_h(NX-1, 0, RIGHTh)]
                + h_new[idx_vec_h(NX-1, 0, DOWNh)];
    inamuro_phi = (right_phi - known_sum_h) / Wh;
    h_new[idx_vec_h(NX-1, 0, LEFTh)] = Wh * inamuro_phi;

    // DOWN LEFT
    h_new[idx_vec_h(0, NY-1, UPh)] = h_old[idx_vec_h(0, NY-2, DOWNh)];
    h_new[idx_vec_h(0, NY-1, LEFTh)] = h_old[idx_vec_h(1, NY-1, LEFTh)];
    h_new[idx_vec_h(0, NY-1, DOWNh)] = h_old[idx_vec_h(0, NY-2, DOWNh)];

    known_sum_h = h_new[idx_vec_h(0, NY-1, UPh)] 
                + h_new[idx_vec_h(0, NY-1, LEFTh)]
                + h_new[idx_vec_h(0, NY-1, DOWNh)];
    inamuro_phi = (left_phi - known_sum_h) / Wh;
    h_new[idx_vec_h(0, NY-1, RIGHTh)] = Wh * inamuro_phi;

    // DOWN RIGHT
    h_new[idx_vec_h(NX-1, NY-1, RIGHTh)] = h_old[idx_vec_h(NX-2, NY-1, RIGHTh)];
    h_new[idx_vec_h(NX-1, NY-1, UPh)] = h_old[idx_vec_h(NX-1, NY-2, DOWNh)];
    h_new[idx_vec_h(NX-1, NY-1, DOWNh)] = h_old[idx_vec_h(NX-1, NY-2, DOWNh)];

    known_sum_h = h_new[idx_vec_h(NX-1, NY-1, UPh)] 
                + h_new[idx_vec_h(NX-1, NY-1, RIGHTh)]
                + h_new[idx_vec_h(NX-1, NY-1, DOWNh)];
    inamuro_phi = (right_phi - known_sum_h) / Wh;
    h_new[idx_vec_h(NX-1, NY-1, LEFTh)] = Wh * inamuro_phi;

    swap(h_new, h_old);
}