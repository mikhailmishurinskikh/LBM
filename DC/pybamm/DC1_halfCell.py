import pybamm

class DC1(pybamm.BaseModel):
    def __init__(self, name="DC1"):
        super().__init__(name)

        # Spatial variables
        self.x = pybamm.SpatialVariable(
            "x", domain=["separator", "anode"], coord_sys="cartesian"
        )


        # Parameters
        # Butler-Folmer coefficients
        Fk = pybamm.Parameter("Fk")
        c_max = pybamm.Parameter("C_max")
        s = pybamm.Parameter("<s>")
        U = lambda sto: pybamm.FunctionParameter("U", {"sto" : sto})
        j0_li = pybamm.Parameter("j0_li")

        # Electrolyte
        D_e = pybamm.Parameter("D_e")
        kappa_e = pybamm.Parameter("kappa_e")
        eps_sep = pybamm.Parameter("eps_sep")
        tau_anode = pybamm.Parameter("tau_anode")
        tau_sep = pybamm.Parameter("tau_sep")
        t_plus = pybamm.Parameter("t_plus")

        # Anode
        kappa_s_eq = pybamm.Parameter("kappa_s_eq")
        eps_s = pybamm.Parameter("eps_s")
        a = pybamm.Parameter("a")
        eps_anode = 1 - eps_s

        # General
        T = pybamm.Parameter("T")
        I = pybamm.FunctionParameter("I", {"t": pybamm.t})


        # Variables
        c_e_sep = pybamm.Variable("c_e_sep", domain="separator")
        c_e_anode = pybamm.Variable("c_e_anode", domain="anode")

        phi_e_sep = pybamm.Variable("phi_e_sep", domain="separator")
        phi_e_anode = pybamm.Variable("phi_e_anode", domain="anode")

        phi_s = pybamm.Variable("phi_s", domain="anode")

        c_s_int = pybamm.Variable("c_s_int", domain="anode")
        c_s_surf = pybamm.Variable("c_s_surf", domain="anode")

        c_e = pybamm.concatenation(c_e_sep, c_e_anode)
        phi_e = pybamm.concatenation(phi_e_sep, phi_e_anode)


        # Transport
        eps = pybamm.concatenation(pybamm.PrimaryBroadcast(eps_sep, "separator"), pybamm.PrimaryBroadcast(eps_anode, "anode"))
        eps_tau = pybamm.concatenation(pybamm.PrimaryBroadcast(eps_sep / tau_sep, "separator"), pybamm.PrimaryBroadcast(eps_anode / tau_anode, "anode"))

        i_e = eps_tau * kappa_e * (
            -pybamm.grad(phi_e) + 2 * (1 - t_plus) * pybamm.constants.R * T / pybamm.constants.F * pybamm.grad(c_e) / c_e
        )
        i_s = -kappa_s_eq * pybamm.grad(phi_s)

        N_e = - eps_tau * D_e * pybamm.grad(c_e) + t_plus / pybamm.constants.F * i_e

        # Reaction anode
        j0 = Fk * pybamm.sqrt(c_e_anode * c_s_surf * (c_max - c_s_surf))
        eta = phi_s - phi_e_anode - U(c_s_surf / c_max)
        j_w = 2 * j0 * pybamm.sinh(pybamm.constants.F / (2 * pybamm.constants.R * T) * eta)
        reaction = pybamm.concatenation(pybamm.PrimaryBroadcast(0, "separator"), a * j_w)
        voltage = pybamm.boundary_value(phi_s, "right")


        # Equations
        self.rhs[c_e] = (-pybamm.div(N_e) + 1 / pybamm.constants.F * reaction) / eps
        self.rhs[c_s_int] = -a / (pybamm.constants.F * eps_s) * j_w

        self.algebraic[phi_e] = pybamm.div(i_e) - reaction
        self.algebraic[phi_s] = pybamm.div(i_s) + a * j_w
        self.algebraic[c_s_surf] = c_s_surf - c_s_int - j_w * s


        # Boundary conditions
        self.boundary_conditions[phi_e] = {
            "left": (-2 * pybamm.constants.R * T / pybamm.constants.F * pybamm.arcsinh(I / (2 * j0_li)), "Dirichlet"),
            "right": (pybamm.Scalar(0), "Neumann"),
        }

        self.boundary_conditions[c_e] = {
            "left":  (-(1 - t_plus) * I / (pybamm.constants.F * eps_sep / tau_sep * D_e), "Neumann"),
            "right": (pybamm.Scalar(0), "Neumann"),
        }

        self.boundary_conditions[phi_s] = {
            "left":  (pybamm.Scalar(0), "Neumann"),
            "right": (-I / kappa_s_eq, "Neumann"),
        }

        self.boundary_conditions[c_s_int] = {
            "left":  (pybamm.Scalar(0), "Neumann"),
            "right": (pybamm.Scalar(0), "Neumann"),
        }

        # Initials conditions
        self.initial_conditions = {
            c_e: pybamm.Parameter("C_e_0"),
            c_s_int: pybamm.Parameter("C_s_0"),
            phi_e : 0,
            phi_s : 4,
            c_s_surf : pybamm.Parameter("C_s_0"),
        }


        # Variables for save
        self.variables = {
            "c_e": c_e,
            "phi_e": phi_e,
            "phi_s": phi_s,
            "c_s_int": c_s_int,
            "c_s_surf": c_s_surf,
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
                pybamm.Parameter("u_max") - voltage
            )
        ]


    @property
    def default_geometry(self):
        return {
            "separator": {
                self.x: {"min": -pybamm.Parameter("L_sep"), "max": 0}
            },
            "anode": {
                self.x: {"min": 0, "max": pybamm.Parameter("L_anode")}
            },
        }


    @property
    def default_spatial_methods(self):
        return {
            "separator": pybamm.FiniteVolume(),
            "anode": pybamm.FiniteVolume(),
        }


    @property
    def default_submesh_types(self):
        return {
            "separator": pybamm.Uniform1DSubMesh,
            "anode": pybamm.Uniform1DSubMesh,
        }


    def run(self, params, x_dots, t_eval):
        pts = {
            self.x: x_dots
        }

        sim = pybamm.Simulation(
            self,
            parameter_values=params,
            var_pts=pts,
            solver=pybamm.IDAKLUSolver(),
        )

        return sim.solve(t_eval)