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

    double* phi;

    double dt;

    // omega = dt / tau
    double omega_h;
    double time;

    double left_phi{1.0};
    double right_phi{0.0};

public:
    Solver(int nx, int ny, double tau_h);
    ~Solver();
    void do_step();
    void save(const char* file_path);

protected:
    inline int idx_vec(int x, int y, int n) const;
    inline int idx_scal(int x, int y) const;
    
    inline void init();
    inline void calc_phi();
    inline void moment_update();
    inline void collision();
    inline void propagation();
};