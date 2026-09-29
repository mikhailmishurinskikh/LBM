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

inline void swap(double*& ptr1, double*& ptr2)
{
    double* tmp = ptr1;
    ptr1 = ptr2;
    ptr2 = tmp;
}

Solver::Solver(int nx, int ny, double tau_h, double tau_g, double F, double kappa, double Deff, double RT)
    : NX(nx), NY(ny), F(F), kappa(kappa), Deff(Deff), RT(RT), tau_h(tau_h), tau_g(tau_g) {
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

void Solver::update_moments() {
    calc_C();
    calc_phi();
}

void Solver::do_step_g() {
    collision_g();
    propagation_g();
    time += dt;
}

void Solver::do_step_h() {
    collision_h();
    propagation_h();
}

void Solver::save(const char* phi_path, const char* C_path)
{
    std::ofstream file_phi(phi_path);
    for (int ny{}; ny < NY; ny++) {
        for (int nx{}; nx < NX-1; nx++) {
            file_phi << phi[idx_scal(nx, ny)] << " ";
        }
        file_phi << phi[idx_scal(NX-1, ny)] << "\n";
    }
    file_phi.close();

    std::ofstream file_C(C_path);
    for (int ny{}; ny < NY; ny++) {
        for (int nx{}; nx < NX-1; nx++) {
            file_C << C[idx_scal(nx, ny)] << " ";
        }
        file_C << C[idx_scal(NX-1, ny)] << "\n";
    }
    file_C.close();
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

inline void Solver::collision_g() {
    // internal dots
    for (int nx{1}; nx < NX - 1; nx++) {
        for (int ny{1}; ny < NY - 1; ny++) {
            double grad_phi1 = (
                phi[idx_scal(nx+1, ny)] - phi[idx_scal(nx-1, ny)]
            ) / 2.0;
            double grad_phi2 = (
                phi[idx_scal(nx, ny+1)] - phi[idx_scal(nx, ny-1)]
            ) / 2.0;
            double F1 = -W14g * (1 - 0.5 * omega_g) * F / RT * grad_phi1 * C[idx_scal(nx, ny)] * dt;
            double F2 = -W14g * (1 - 0.5 * omega_g) * F / RT * grad_phi2 * C[idx_scal(nx, ny)] * dt;
            double F3 = -F1;
            double F4 = -F2;

            double g_eq0 = W0g * C[idx_scal(nx, ny)];
            double g_eq14 = W14g * C[idx_scal(nx, ny)];
            
            g_old[idx_vec_g(nx, ny, STAYg)]   = (1.0 - omega_g) * g_old[idx_vec_g(nx, ny, STAYg)]   + omega_g * g_eq0;
            g_old[idx_vec_g(nx, ny, RIGHTg)]  = (1.0 - omega_g) * g_old[idx_vec_g(nx, ny, RIGHTg)]  + omega_g * g_eq14 + F1;
            g_old[idx_vec_g(nx, ny, UPg)]     = (1.0 - omega_g) * g_old[idx_vec_g(nx, ny, UPg)]     + omega_g * g_eq14 + F2;
            g_old[idx_vec_g(nx, ny, LEFTg)]   = (1.0 - omega_g) * g_old[idx_vec_g(nx, ny, LEFTg)]   + omega_g * g_eq14 + F3;
            g_old[idx_vec_g(nx, ny, DOWNg)]   = (1.0 - omega_g) * g_old[idx_vec_g(nx, ny, DOWNg)]   + omega_g * g_eq14 + F4;
        }
    }

    // upper wall
    for (int nx{1}; nx < NX - 1; nx++) {
        double grad_phi1 = (
            phi[idx_scal(nx+1, 0)] - phi[idx_scal(nx-1, 0)]
        ) / 2.0;
        double F1 = -W14g * (1 - 0.5 * omega_g) * F / RT * grad_phi1 * C[idx_scal(nx, 0)] * dt;
        double F3 = -F1;

        double g_eq0 = W0g * C[idx_scal(nx, 0)];
        double g_eq14 = W14g * C[idx_scal(nx, 0)];
        
        g_old[idx_vec_g(nx, 0, STAYg)]    = (1.0 - omega_g) * g_old[idx_vec_g(nx, 0, STAYg)]    + omega_g * g_eq0;
        g_old[idx_vec_g(nx, 0, RIGHTg)]   = (1.0 - omega_g) * g_old[idx_vec_g(nx, 0, RIGHTg)]   + omega_g * g_eq14 + F1;
        g_old[idx_vec_g(nx, 0, UPg)]      = (1.0 - omega_g) * g_old[idx_vec_g(nx, 0, UPg)]      + omega_g * g_eq14;
        g_old[idx_vec_g(nx, 0, LEFTg)]    = (1.0 - omega_g) * g_old[idx_vec_g(nx, 0, LEFTg)]    + omega_g * g_eq14 + F3;
        g_old[idx_vec_g(nx, 0, DOWNg)]    = (1.0 - omega_g) * g_old[idx_vec_g(nx, 0, DOWNg)]    + omega_g * g_eq14;
    }

    // bottom wall
    for (int nx{1}; nx < NX - 1; nx++) {
        double grad_phi1 = (
            phi[idx_scal(nx+1, NY-1)] - phi[idx_scal(nx-1, NY-1)]
        ) / 2.0;
        double F1 = -W14g * (1 - 0.5 * omega_g) * F / RT * grad_phi1 * C[idx_scal(nx, NY-1)] * dt;
        double F3 = -F1;

        double g_eq0 = W0g * C[idx_scal(nx, NY-1)];
        double g_eq14 = W14g * C[idx_scal(nx, NY-1)];
        
        g_old[idx_vec_g(nx, NY-1, STAYg)]  = (1.0 - omega_g) * g_old[idx_vec_g(nx, NY-1, STAYg)]  + omega_g * g_eq0;
        g_old[idx_vec_g(nx, NY-1, RIGHTg)] = (1.0 - omega_g) * g_old[idx_vec_g(nx, NY-1, RIGHTg)] + omega_g * g_eq14 + F1;
        g_old[idx_vec_g(nx, NY-1, UPg)]    = (1.0 - omega_g) * g_old[idx_vec_g(nx, NY-1, UPg)]    + omega_g * g_eq14;
        g_old[idx_vec_g(nx, NY-1, LEFTg)]  = (1.0 - omega_g) * g_old[idx_vec_g(nx, NY-1, LEFTg)]  + omega_g * g_eq14 + F3;
        g_old[idx_vec_g(nx, NY-1, DOWNg)]  = (1.0 - omega_g) * g_old[idx_vec_g(nx, NY-1, DOWNg)]  + omega_g * g_eq14;
    }

    // left wall
    for (int ny{1}; ny < NY - 1; ny++) {
        double grad_phi1 = (
            phi[idx_scal(1, ny)] - phi[idx_scal(0, ny)]
        );
        double grad_phi2 = (
            phi[idx_scal(0, ny+1)] - phi[idx_scal(0, ny-1)]
        ) / 2.0;
        double F1 = -W14g * (1 - 0.5 * omega_g) * F / RT * grad_phi1 * C[idx_scal(0, ny)] * dt;
        double F2 = -W14g * (1 - 0.5 * omega_g) * F / RT * grad_phi2 * C[idx_scal(0, ny)] * dt;
        double F3 = -F1;
        double F4 = -F2;

        double g_eq0 = W0g * C[idx_scal(0, ny)];
        double g_eq14 = W14g * C[idx_scal(0, ny)];
        
        g_old[idx_vec_g(0, ny, STAYg)]   = (1.0 - omega_g) * g_old[idx_vec_g(0, ny, STAYg)]   + omega_g * g_eq0;
        g_old[idx_vec_g(0, ny, RIGHTg)]  = (1.0 - omega_g) * g_old[idx_vec_g(0, ny, RIGHTg)]  + omega_g * g_eq14 + F1;
        g_old[idx_vec_g(0, ny, UPg)]     = (1.0 - omega_g) * g_old[idx_vec_g(0, ny, UPg)]     + omega_g * g_eq14 + F2;
        g_old[idx_vec_g(0, ny, LEFTg)]   = (1.0 - omega_g) * g_old[idx_vec_g(0, ny, LEFTg)]   + omega_g * g_eq14 + F3;
        g_old[idx_vec_g(0, ny, DOWNg)]   = (1.0 - omega_g) * g_old[idx_vec_g(0, ny, DOWNg)]   + omega_g * g_eq14 + F4;
    }

    // right wall
    for (int ny{1}; ny < NY - 1; ny++) {
        double grad_phi1 = (
            phi[idx_scal(NX-1, ny)] - phi[idx_scal(NX-2, ny)]
        );
        double grad_phi2 = (
            phi[idx_scal(NX-1, ny+1)] - phi[idx_scal(NX-1, ny-1)]
        ) / 2.0;
        double F1 = -W14g * (1 - 0.5 * omega_g) * F / RT * grad_phi1 * C[idx_scal(NX-1, ny)] * dt;
        double F2 = -W14g * (1 - 0.5 * omega_g) * F / RT * grad_phi2 * C[idx_scal(NX-1, ny)] * dt;
        double F3 = -F1;
        double F4 = -F2;

        double g_eq0 = W0g * C[idx_scal(NX-1, ny)];
        double g_eq14 = W14g * C[idx_scal(NX-1, ny)];
        
        g_old[idx_vec_g(NX-1, ny, STAYg)]   = (1.0 - omega_g) * g_old[idx_vec_g(NX-1, ny, STAYg)]   + omega_g * g_eq0;
        g_old[idx_vec_g(NX-1, ny, RIGHTg)]  = (1.0 - omega_g) * g_old[idx_vec_g(NX-1, ny, RIGHTg)]  + omega_g * g_eq14 + F1;
        g_old[idx_vec_g(NX-1, ny, UPg)]     = (1.0 - omega_g) * g_old[idx_vec_g(NX-1, ny, UPg)]     + omega_g * g_eq14 + F2;
        g_old[idx_vec_g(NX-1, ny, LEFTg)]   = (1.0 - omega_g) * g_old[idx_vec_g(NX-1, ny, LEFTg)]   + omega_g * g_eq14 + F3;
        g_old[idx_vec_g(NX-1, ny, DOWNg)]   = (1.0 - omega_g) * g_old[idx_vec_g(NX-1, ny, DOWNg)]   + omega_g * g_eq14 + F4;
    }

    // corners
    // UP LEFT
    {
        double grad_phi1 = (
            phi[idx_scal(1, 0)] - phi[idx_scal(0, 0)]
        );
        double F1 = -W14g * (1 - 0.5 * omega_g) * F / RT * grad_phi1 * C[idx_scal(0, 0)] * dt;
        double F3 = -F1;

        double g_eq0 = W0g * C[idx_scal(0, 0)];
        double g_eq14 = W14g * C[idx_scal(0, 0)];

        g_old[idx_vec_g(0, 0, STAYg)]   = (1.0 - omega_g) * g_old[idx_vec_g(0, 0, STAYg)]   + omega_g * g_eq0;
        g_old[idx_vec_g(0, 0, RIGHTg)]  = (1.0 - omega_g) * g_old[idx_vec_g(0, 0, RIGHTg)]  + omega_g * g_eq14 + F1;
        g_old[idx_vec_g(0, 0, UPg)]     = (1.0 - omega_g) * g_old[idx_vec_g(0, 0, UPg)]     + omega_g * g_eq14;
        g_old[idx_vec_g(0, 0, LEFTg)]   = (1.0 - omega_g) * g_old[idx_vec_g(0, 0, LEFTg)]   + omega_g * g_eq14 + F3;
        g_old[idx_vec_g(0, 0, DOWNg)]   = (1.0 - omega_g) * g_old[idx_vec_g(0, 0, DOWNg)]   + omega_g * g_eq14;
    }

    // UP RIGHT
    {
        double grad_phi1 = (
            phi[idx_scal(NX-1, 0)] - phi[idx_scal(NX-2, 0)]
        );
        double F1 = -W14g * (1 - 0.5 * omega_g) * F / RT * grad_phi1 * C[idx_scal(NX-1, 0)] * dt;
        double F3 = -F1;

        double g_eq0 = W0g * C[idx_scal(NX-1, 0)];
        double g_eq14 = W14g * C[idx_scal(NX-1, 0)];

        g_old[idx_vec_g(NX-1, 0, STAYg)]   = (1.0 - omega_g) * g_old[idx_vec_g(NX-1, 0, STAYg)]   + omega_g * g_eq0;
        g_old[idx_vec_g(NX-1, 0, RIGHTg)]  = (1.0 - omega_g) * g_old[idx_vec_g(NX-1, 0, RIGHTg)]  + omega_g * g_eq14 + F1;
        g_old[idx_vec_g(NX-1, 0, UPg)]     = (1.0 - omega_g) * g_old[idx_vec_g(NX-1, 0, UPg)]     + omega_g * g_eq14;
        g_old[idx_vec_g(NX-1, 0, LEFTg)]   = (1.0 - omega_g) * g_old[idx_vec_g(NX-1, 0, LEFTg)]   + omega_g * g_eq14 + F3;
        g_old[idx_vec_g(NX-1, 0, DOWNg)]   = (1.0 - omega_g) * g_old[idx_vec_g(NX-1, 0, DOWNg)]   + omega_g * g_eq14;
    }

    // DOWN LEFT
    {
        double grad_phi1 = (
            phi[idx_scal(1, NY-1)] - phi[idx_scal(0, NY-1)]
        );
        double F1 = -W14g * (1 - 0.5 * omega_g) * F / RT * grad_phi1 * C[idx_scal(0, NY-1)] * dt;
        double F3 = -F1;

        double g_eq0 = W0g * C[idx_scal(0, NY-1)];
        double g_eq14 = W14g * C[idx_scal(0, NY-1)];

        g_old[idx_vec_g(0, NY-1, STAYg)]   = (1.0 - omega_g) * g_old[idx_vec_g(0, NY-1, STAYg)]   + omega_g * g_eq0;
        g_old[idx_vec_g(0, NY-1, RIGHTg)]  = (1.0 - omega_g) * g_old[idx_vec_g(0, NY-1, RIGHTg)]  + omega_g * g_eq14 + F1;
        g_old[idx_vec_g(0, NY-1, UPg)]     = (1.0 - omega_g) * g_old[idx_vec_g(0, NY-1, UPg)]     + omega_g * g_eq14;
        g_old[idx_vec_g(0, NY-1, LEFTg)]   = (1.0 - omega_g) * g_old[idx_vec_g(0, NY-1, LEFTg)]   + omega_g * g_eq14 + F3;
        g_old[idx_vec_g(0, NY-1, DOWNg)]   = (1.0 - omega_g) * g_old[idx_vec_g(0, NY-1, DOWNg)]   + omega_g * g_eq14;
    }

    // DOWN RIGHT
    {
        double grad_phi1 = (
            phi[idx_scal(NX-1, NY-1)] - phi[idx_scal(NX-2, NY-1)]
        );
        double F1 = -W14g * (1 - 0.5 * omega_g) * F / RT * grad_phi1 * C[idx_scal(NX-1, NY-1)] * dt;
        double F3 = -F1;

        double g_eq0 = W0g * C[idx_scal(NX-1, NY-1)];
        double g_eq14 = W14g * C[idx_scal(NX-1, NY-1)];

        g_old[idx_vec_g(NX-1, NY-1, STAYg)]   = (1.0 - omega_g) * g_old[idx_vec_g(NX-1, NY-1, STAYg)]   + omega_g * g_eq0;
        g_old[idx_vec_g(NX-1, NY-1, RIGHTg)]  = (1.0 - omega_g) * g_old[idx_vec_g(NX-1, NY-1, RIGHTg)]  + omega_g * g_eq14 + F1;
        g_old[idx_vec_g(NX-1, NY-1, UPg)]     = (1.0 - omega_g) * g_old[idx_vec_g(NX-1, NY-1, UPg)]     + omega_g * g_eq14;
        g_old[idx_vec_g(NX-1, NY-1, LEFTg)]   = (1.0 - omega_g) * g_old[idx_vec_g(NX-1, NY-1, LEFTg)]   + omega_g * g_eq14 + F3;
        g_old[idx_vec_g(NX-1, NY-1, DOWNg)]   = (1.0 - omega_g) * g_old[idx_vec_g(NX-1, NY-1, DOWNg)]   + omega_g * g_eq14;
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
        h_new[idx_vec_h(nx, 0, UPh)] = h_old[idx_vec_h(nx, 1, UPh)];
        h_new[idx_vec_h(nx, 0, LEFTh)] = h_old[idx_vec_h(nx+1, 0, LEFTh)];
        h_new[idx_vec_h(nx, 0, DOWNh)] = h_old[idx_vec_h(nx, 1, UPh)];

        // lower wall
        h_new[idx_vec_h(nx, NY - 1, RIGHTh)] = h_old[idx_vec_h(nx-1, NY-1, RIGHTh)];
        h_new[idx_vec_h(nx, NY - 1, DOWNh)] = h_old[idx_vec_h(nx, NY-2, DOWNh)];
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

inline void Solver::propagation_g() {
    for (int nx{1}; nx < NX-1; nx++) {
        for (int ny{1}; ny < NY-1; ny++) {    
            g_new[idx_vec_g(nx, ny, STAYg)] = g_old[idx_vec_g(nx, ny, STAYg)];        
            g_new[idx_vec_g(nx, ny, RIGHTg)] = g_old[idx_vec_g(nx-1, ny, RIGHTg)];
            g_new[idx_vec_g(nx, ny, UPg)] = g_old[idx_vec_g(nx, ny+1, UPg)];
            g_new[idx_vec_g(nx, ny, LEFTg)] = g_old[idx_vec_g(nx+1, ny, LEFTg)];
            g_new[idx_vec_g(nx, ny, DOWNg)] = g_old[idx_vec_g(nx, ny-1, DOWNg)];
        }
    }

    // bounce-back conditions (wet-node)
    // upper wall
    for (int nx = 1; nx < NX - 1; nx++) {
        g_new[idx_vec_g(nx, 0, STAYg)] = g_old[idx_vec_g(nx, 0, STAYg)];
        g_new[idx_vec_g(nx, 0, RIGHTg)] = g_old[idx_vec_g(nx-1, 0, RIGHTg)];
        g_new[idx_vec_g(nx, 0, UPg)] = g_old[idx_vec_g(nx, 1, UPg)];
        g_new[idx_vec_g(nx, 0, LEFTg)] = g_old[idx_vec_g(nx+1, 0, LEFTg)];
        g_new[idx_vec_g(nx, 0, DOWNg)] = g_old[idx_vec_g(nx, 1, UPg)];
    }

    // bottom wall
    for (int nx = 1; nx < NX - 1; nx++) {
        g_new[idx_vec_g(nx, NY - 1, STAYg)] = g_old[idx_vec_g(nx, NY - 1, STAYg)];
        g_new[idx_vec_g(nx, NY - 1, RIGHTg)] = g_old[idx_vec_g(nx-1, NY - 1, RIGHTg)];
        g_new[idx_vec_g(nx, NY - 1, DOWNg)] = g_old[idx_vec_g(nx, NY - 2, DOWNg)];
        g_new[idx_vec_g(nx, NY - 1, LEFTg)] = g_old[idx_vec_g(nx+1, NY - 1, LEFTg)];
        g_new[idx_vec_g(nx, NY - 1, UPg)] = g_old[idx_vec_g(nx, NY - 2, DOWNg)];
    }
    
    // left wall
    for (int ny = 1; ny < NY - 1; ny++) {
        g_new[idx_vec_g(0, ny, STAYg)] = g_old[idx_vec_g(0, ny, STAYg)];
        g_new[idx_vec_g(0, ny, RIGHTg)] = g_old[idx_vec_g(1, ny, LEFTg)];
        g_new[idx_vec_g(0, ny, UPg)] = g_old[idx_vec_g(0, ny+1, UPg)];
        g_new[idx_vec_g(0, ny, LEFTg)] = g_old[idx_vec_g(1, ny, LEFTg)];
        g_new[idx_vec_g(0, ny, DOWNg)] = g_old[idx_vec_g(0, ny-1, DOWNg)];
    }

    // right wall
    for (int ny = 1; ny < NY - 1; ny++) {
        g_new[idx_vec_g(NX-1, ny, STAYg)] = g_old[idx_vec_g(NX-1, ny, STAYg)];
        g_new[idx_vec_g(NX-1, ny, RIGHTg)] = g_old[idx_vec_g(NX-2, ny, RIGHTg)];
        g_new[idx_vec_g(NX-1, ny, UPg)] = g_old[idx_vec_g(NX-1, ny+1, UPg)];
        g_new[idx_vec_g(NX-1, ny, LEFTg)] = g_old[idx_vec_g(NX-2, ny, RIGHTg)];
        g_new[idx_vec_g(NX-1, ny, DOWNg)] = g_old[idx_vec_g(NX-1, ny-1, DOWNg)];
    }

    // corners
    // UP LEFT
    g_new[idx_vec_g(0, 0, STAYg)] = g_old[idx_vec_g(0, 0, STAYg)];
    g_new[idx_vec_g(0, 0, RIGHTg)] = g_old[idx_vec_g(1, 0, LEFTg)];
    g_new[idx_vec_g(0, 0, UPg)] = g_old[idx_vec_g(0, 1, UPg)];
    g_new[idx_vec_g(0, 0, LEFTg)] = g_old[idx_vec_g(1, 0, LEFTg)];
    g_new[idx_vec_g(0, 0, DOWNg)] = g_old[idx_vec_g(0, 1, UPg)];

    // UP RIGHT
    g_new[idx_vec_g(NX-1, 0, STAYg)] = g_old[idx_vec_g(NX-1, 0, STAYg)];
    g_new[idx_vec_g(NX-1, 0, RIGHTg)] = g_old[idx_vec_g(NX-2, 0, RIGHTg)];
    g_new[idx_vec_g(NX-1, 0, UPg)] = g_old[idx_vec_g(NX-1, 1, UPg)];
    g_new[idx_vec_g(NX-1, 0, LEFTg)] = g_old[idx_vec_g(NX-2, 0, RIGHTg)];
    g_new[idx_vec_g(NX-1, 0, DOWNg)] = g_old[idx_vec_g(NX-1, 1, UPg)];

    // DOWN LEFT
    g_new[idx_vec_g(0, NY-1, STAYg)] = g_old[idx_vec_g(0, NY-1, STAYg)];
    g_new[idx_vec_g(0, NY-1, RIGHTg)] = g_old[idx_vec_g(1, NY-1, LEFTg)];
    g_new[idx_vec_g(0, NY-1, UPg)] = g_old[idx_vec_g(0, NY-2, DOWNg)];
    g_new[idx_vec_g(0, NY-1, LEFTg)] = g_old[idx_vec_g(1, NY-1, LEFTg)];
    g_new[idx_vec_g(0, NY-1, DOWNg)] = g_old[idx_vec_g(0, NY-2, DOWNg)];

    // DOWN RIGHT
    g_new[idx_vec_g(NX-1, NY-1, STAYg)] = g_old[idx_vec_g(NX-1, NY-1, STAYg)];
    g_new[idx_vec_g(NX-1, NY-1, RIGHTg)] = g_old[idx_vec_g(NX-2, NY-1, RIGHTg)];
    g_new[idx_vec_g(NX-1, NY-1, UPg)] = g_old[idx_vec_g(NX-1, NY-2, DOWNg)];
    g_new[idx_vec_g(NX-1, NY-1, LEFTg)] = g_old[idx_vec_g(NX-2, NY-1, RIGHTg)];
    g_new[idx_vec_g(NX-1, NY-1, DOWNg)] = g_old[idx_vec_g(NX-1, NY-2, DOWNg)];

    swap(g_new, g_old);
}