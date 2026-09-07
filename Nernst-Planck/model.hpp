inline void swap(double*& ptr1, double*& ptr2);

class Solver {
    /*
    Capacitor: left - anode, right - catode
    */
protected:

    // Grid size
    int NX;
    int NY;

    // Electric potential
    double* h_new;
    double* h_old;

    // Concentration of ions
    double* g_new;
    double* g_old;

    // macro parameters
    double* phi;
    double* C;

    // equation parameters
    double A, kappa, D;

    double dt;

    // omega = dt / tau
    double omega_h;
    double omega_g;
    
    double time;

    double left_phi{1.0};
    double right_phi{0.0};

public:
    Solver(int nx, int ny, double tau_h, double tau_g, double A, double kappa, double D);
    ~Solver();
    void do_step();
    void save(const char* file_path);

protected:
    inline int idx_vec_g(int x, int y, int n) const;
    inline int idx_vec_h(int x, int y, int n) const;
    inline int idx_scal(int x, int y) const;
    
    inline void init();
    inline void calc_phi();
    inline void calc_C();
    inline void moment_update();
    inline void collision();
    inline void propagation();
};