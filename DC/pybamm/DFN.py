import pybamm

class DFN(pybamm.BaseModel):
    def __init__(self, name="DFN"):
        super().__init__(name)
        
        # Spatial variables
        self.x = pybamm.SpatialVariable(
            "x", domain=["anode", "separator", "cathode"], coord_sys="cartesian"
        )

        self.r_anode = pybamm.SpatialVariable(
            "r_anode", domain=["anode particle"], auxiliary_domains={"secondary": "anode"},
            coord_sys="spherical polar",
        )

        self.r_cathode = pybamm.SpatialVariable(
            "r_cathode", domain=["cathode particle"], auxiliary_domains={"secondary": "cathode"},
            coord_sys="spherical polar",
        )


        # Parameters
        
        # Butler-Folmer coefficients anode (C)
        Fk_anode = pybamm.Parameter("Fk_anode")
        c_max_anode = pybamm.Parameter("C_max_anode")
        U_anode = lambda sto: pybamm.FunctionParameter("U_anode", {"sto" : sto})
        
        
        # Butler-Folmer coefficients cathode (NMC)
        Fk_cathode = pybamm.Parameter("Fk_cathode")
        c_max_cathode = pybamm.Parameter("C_max_cathode")
        U_cathode = lambda sto: pybamm.FunctionParameter("U_cathode", {"sto" : sto})
        

        # Electrolyte
        D_e = lambda c: pybamm.FunctionParameter("D_e", {"c" : c})
        kappa_e = lambda c: pybamm.FunctionParameter("kappa_e", {"c" : c})
        eps_sep = pybamm.Parameter("eps_sep")
        tau_anode = pybamm.Parameter("tau_anode")
        tau_sep = pybamm.Parameter("tau_sep")
        tau_cathode = pybamm.Parameter("tau_cathode")
        t_plus = pybamm.Parameter("t_plus")

        # Anode
        kappa_s_anode = pybamm.Parameter("kappa_s_anode")
        eps_s_anode = pybamm.Parameter("eps_s_anode")
        self.R_p_anode = pybamm.Parameter("R_p_anode")
        D_s_anode = pybamm.Parameter("D_s_anode")
        a_anode = 3 * eps_s_anode / self.R_p_anode
        eps_anode = 1 - eps_s_anode

        # Cathode
        kappa_s_cathode = pybamm.Parameter("kappa_s_cathode")
        eps_s_cathode = pybamm.Parameter("eps_s_cathode")
        self.R_p_cathode = pybamm.Parameter("R_p_cathode")
        D_s_cathode = pybamm.Parameter("D_s_cathode")
        a_cathode = 3 * eps_s_cathode / self.R_p_cathode
        eps_cathode = 1 - eps_s_cathode

        # General
        T = pybamm.Parameter("T")
        I = pybamm.FunctionParameter("I", {"t": pybamm.t})


        # Variables
        c_e_anode = pybamm.Variable("c_e_anode", domain="anode")
        c_e_sep = pybamm.Variable("c_e_sep", domain="separator")
        c_e_cathode = pybamm.Variable("c_e_cathode", domain="cathode")

        phi_e_anode = pybamm.Variable("phi_e_anode", domain="anode")
        phi_e_sep = pybamm.Variable("phi_e_sep", domain="separator")
        phi_e_cathode = pybamm.Variable("phi_e_cathode", domain="cathode")

        phi_s_anode = pybamm.Variable("phi_s_anode", domain="anode")
        phi_s_cathode = pybamm.Variable("phi_s_cathode", domain="cathode")

        c_s_anode = pybamm.Variable(
            "c_s_anode",
            domain="anode particle",
            auxiliary_domains={"secondary": "anode"},
        )

        c_s_cathode = pybamm.Variable(
            "c_s_cathode",
            domain="cathode particle",
            auxiliary_domains={"secondary": "cathode"},
        )

        c_e = pybamm.concatenation(c_e_anode, c_e_sep, c_e_cathode)
        phi_e = pybamm.concatenation(phi_e_anode, phi_e_sep, phi_e_cathode)
        c_s_surf_anode = pybamm.surf(c_s_anode)
        c_s_surf_cathode = pybamm.surf(c_s_cathode)


        # Transport
        eps = pybamm.concatenation(
            pybamm.PrimaryBroadcast(eps_anode, "anode"),
            pybamm.PrimaryBroadcast(eps_sep, "separator"),
            pybamm.PrimaryBroadcast(eps_cathode, "cathode"),
        )
        eps_tau = pybamm.concatenation(
            pybamm.PrimaryBroadcast(eps_anode / tau_anode, "anode"),
            pybamm.PrimaryBroadcast(eps_sep / tau_sep, "separator"),
            pybamm.PrimaryBroadcast(eps_cathode / tau_cathode, "cathode"),
        )

        i_e = eps_tau * kappa_e(c_e) * (
            -pybamm.grad(phi_e) + 2 * (1 - t_plus) * pybamm.constants.R * T / pybamm.constants.F * pybamm.grad(c_e) / c_e
        )
        i_s_anode = -kappa_s_anode * pybamm.grad(phi_s_anode)
        i_s_cathode = -kappa_s_cathode * pybamm.grad(phi_s_cathode)

        N_e = - eps_tau * D_e(c_e) * pybamm.grad(c_e) + t_plus / pybamm.constants.F * i_e
        N_s_anode = -D_s_anode * pybamm.grad(c_s_anode)
        N_s_cathode = -D_s_cathode * pybamm.grad(c_s_cathode)

        # Reaction anode
        j0_anode = Fk_anode * pybamm.sqrt(c_e_anode * c_s_surf_anode * (c_max_anode - c_s_surf_anode))
        eta_anode = phi_s_anode - phi_e_anode - U_anode(c_s_surf_anode / c_max_anode)
        j_anode = 2 * j0_anode * pybamm.sinh(pybamm.constants.F / (2 * pybamm.constants.R * T) * eta_anode)

        # Reaction cathode
        j0_cathode = Fk_cathode * pybamm.sqrt(c_e_cathode * c_s_surf_cathode * (c_max_cathode - c_s_surf_cathode))
        eta_cathode = phi_s_cathode - phi_e_cathode - U_cathode(c_s_surf_cathode / c_max_cathode)
        j_cathode = 2 * j0_cathode * pybamm.sinh(pybamm.constants.F / (2 * pybamm.constants.R * T) * eta_cathode)

        reaction = pybamm.concatenation(
            a_anode * j_anode,
            pybamm.PrimaryBroadcast(0, "separator"),
            a_cathode * j_cathode,
        )
        voltage = pybamm.boundary_value(phi_s_cathode, "right") - pybamm.boundary_value(phi_s_anode, "left")


        # Equations
        self.rhs[c_e] = (-pybamm.div(N_e) + 1 / pybamm.constants.F * reaction) / eps
        self.rhs[c_s_anode] = -pybamm.div(N_s_anode)
        self.rhs[c_s_cathode] = -pybamm.div(N_s_cathode)

        self.algebraic[phi_e] = pybamm.div(i_e) - reaction
        self.algebraic[phi_s_anode] = pybamm.div(i_s_anode) + a_anode * j_anode
        self.algebraic[phi_s_cathode] = pybamm.div(i_s_cathode) + a_cathode * j_cathode


        # Boundary conditions
        self.boundary_conditions[phi_e] = {
            "left": (pybamm.Scalar(0), "Neumann"),
            "right": (pybamm.Scalar(0), "Neumann"),
        }

        self.boundary_conditions[c_e] = {
            "left":  (pybamm.Scalar(0), "Neumann"),
            "right": (pybamm.Scalar(0), "Neumann"),
        }

        self.boundary_conditions[phi_s_anode] = {
            "left":  (pybamm.Scalar(0), "Dirichlet"),
            "right": (pybamm.Scalar(0), "Neumann"),
        }

        self.boundary_conditions[phi_s_cathode] = {
            "left":  (pybamm.Scalar(0), "Neumann"),
            "right": (-I / kappa_s_cathode, "Neumann"),
        }

        self.boundary_conditions[c_s_anode] = {
            "left":  (pybamm.Scalar(0), "Neumann"),
            "right": (-j_anode / (pybamm.constants.F * D_s_anode), "Neumann"),
        }

        self.boundary_conditions[c_s_cathode] = {
            "left":  (pybamm.Scalar(0), "Neumann"),
            "right": (-j_cathode / (pybamm.constants.F * D_s_cathode), "Neumann"),
        }

        # Initials conditions
        C_e_0 = pybamm.Parameter("C_e_0")
        C_s_0_anode = pybamm.Parameter("C_s_0_anode")
        C_s_0_cathode = pybamm.Parameter("C_s_0_cathode")
        U_anode_0 = U_anode(C_s_0_anode / c_max_anode)
        U_cathode_0 = U_cathode(C_s_0_cathode / c_max_cathode)

        self.initial_conditions = {
            c_e: C_e_0,
            c_s_anode: C_s_0_anode,
            c_s_cathode: C_s_0_cathode,
            phi_e : -U_anode_0,
            phi_s_anode: pybamm.Scalar(0),
            phi_s_cathode: U_cathode_0 - U_anode_0,
        }


        # Variables for save
        self.variables = {
            "c_e": c_e,
            "phi_e": phi_e,
            "phi_s_anode": phi_s_anode,
            "phi_s_cathode": phi_s_cathode,
            "c_s_anode": c_s_anode,
            "c_s_cathode": c_s_cathode,
            "c_s_anode_int": pybamm.r_average(c_s_anode),
            "c_s_cathode_int": pybamm.r_average(c_s_cathode),
            "c_s_surf_anode": c_s_surf_anode,
            "c_s_surf_cathode": c_s_surf_cathode,
            "voltage": voltage
        }


        # Cutoff voltage
        self.events = [
            pybamm.Event(
                "Minimum voltage",
                voltage - pybamm.Parameter("u_min"),
            ),
            
            pybamm.Event(
                "Maximum voltage",
                -voltage + pybamm.Parameter("u_max")
            )
        ]
        
    
    @property
    def default_geometry(self):
        L_anode = pybamm.Parameter("L_anode")
        L_sep = pybamm.Parameter("L_sep")
        L_cathode = pybamm.Parameter("L_cathode")

        return {
            "anode": {
                self.x: {"min": 0, "max": L_anode}
            },
            "separator": {
                self.x: {"min": L_anode, "max": L_anode + L_sep}
            },
            "cathode": {
                self.x: {"min": L_anode + L_sep, "max": L_anode + L_sep + L_cathode}
            },
            "anode particle": {
                self.r_anode: {"min": 0, "max": self.R_p_anode}
            },
            "cathode particle": {
                self.r_cathode: {"min": 0, "max": self.R_p_cathode}
            },
        }
        
    
    @property
    def default_spatial_methods(self):
        return {
            "anode": pybamm.FiniteVolume(),
            "separator": pybamm.FiniteVolume(),
            "cathode": pybamm.FiniteVolume(),
            "anode particle": pybamm.FiniteVolume(),
            "cathode particle": pybamm.FiniteVolume(),
        }
        
        
    @property
    def default_submesh_types(self):
        return {
            "anode": pybamm.Uniform1DSubMesh,
            "separator": pybamm.Uniform1DSubMesh,
            "cathode": pybamm.Uniform1DSubMesh,
            "anode particle": pybamm.Uniform1DSubMesh,
            "cathode particle": pybamm.Uniform1DSubMesh,
        }
        
        
    def run(self, params, x_dots, r_dots, t_eval):
        pts = {
            self.x: x_dots,
            self.r_anode: r_dots,
            self.r_cathode: r_dots,
        }
        
        sim = pybamm.Simulation(
            self,
            parameter_values=params,
            var_pts=pts,
            solver=pybamm.IDAKLUSolver(),
        )
        
        return sim.solve(t_eval)