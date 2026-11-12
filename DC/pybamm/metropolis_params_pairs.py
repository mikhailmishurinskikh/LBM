def Metropolis(names, step):
    import pybamm
    import numpy as np
    from DFN import DFN
    
    data_exp = np.loadtxt("exp_data_1C.txt")
    t_eval, v_exp = data_exp[:, 0], data_exp[:, 1]
    
    params = pybamm.ParameterValues({
        # Geometry
        "L_sep": 12e-6,
        "L_anode": 85.2e-6,
        "L_cathode": 75.6e-6,

        # Butler-Folmer anode, negative (C)
        "Fk_anode": "[input]",
        "C_max_anode": 33133,
        "U_anode": lambda x: (
            1.9793 * pybamm.exp(-39.3631 * x) + 0.2482
            - 0.0909 * pybamm.tanh(29.8538 * (x - 0.1234))
            - 0.04478 * pybamm.tanh(14.9159 * (x - 0.2769))
            - 0.0205 * pybamm.tanh(30.4444 * (x - 0.6103))
        ),

        # Butler-Folmer cathode, positive (NMC)
        "Fk_cathode": "[input]",
        "C_max_cathode": 63104,
        "U_cathode": lambda x: (
            -0.8090 * x + 4.4875
            - 0.0428 * pybamm.tanh(18.5138 * (x - 0.5542))
            - 17.7326 * pybamm.tanh(15.7890 * (x - 0.3117))
            + 17.5842 * pybamm.tanh(15.9308 * (x - 0.3120))
        ),

        # Electrolyte
        "D_e": lambda c: 8.794e-11 * (c/1000) ** 2 - 3.972e-10 * (c/1000) + 4.862e-10,
        "kappa_e": lambda c: 0.1297 * (c/1000) ** 3 - 2.51 * (c/1000) ** 1.5 + 3.329 * (c/1000),
        "t_plus": 0.2594,
        "eps_sep": 0.47,
        "tau_sep": 0.47 / 0.47**1.5,
        "tau_anode": 0.25 / 0.25**1.5,
        "tau_cathode": 0.335 / 0.335**1.5,

        # Solid anode
        "eps_s_anode": 0.75,
        "kappa_s_anode": 215,
        "D_s_anode": "[input]",
        "R_p_anode": 5.86e-6,

        # Solid cathode
        "eps_s_cathode": 0.665,
        "kappa_s_cathode": 0.18,
        "D_s_cathode": "[input]",
        "R_p_cathode": 5.22e-6,

        # General
        "T": 298.15,

        # Initial
        "C_e_0": 1000,
        "C_s_0_anode": 0.9014 * 33133,
        "C_s_0_cathode": 0.27 * 63104,
        
        # Cutoff
        "u_min" : 2.5,
        "u_max" : 4.4,
        
        "I": lambda t: 5 / 0.1027
    })
    
    inputs = {
        "Fk_anode": 6.48e-7,
        "Fk_cathode": 3.42e-6,
        "D_s_anode": 3.3e-14,
        "D_s_cathode": 4e-15,
    }

    model = DFN()
    model.init(params, 30, 30)

    # Parameters
    n_steps = 1000
    a1 = 1 / 20
    a2 = 1 / 0.02
    t_rmse = 3000.0
    accept_th = 2.53


    rng = np.random.default_rng(0)
    t_end_exp = t_eval[-1]
    t_run = np.linspace(0, 1.3 * t_end_exp, 1000)
    mask_rmse = t_eval <= t_rmse


    def evaluate(theta):
        run_inputs = dict(inputs, **{k: 10 ** t for k, t in zip(names, theta)})
        try:
            sol = model.run(inputs=run_inputs, t_eval=t_run, sens=False)
        except Exception as e:
            return None
        t_end = sol.t[-1]
        time_error = abs(t_end - t_end_exp)
        if t_end < t_rmse:
            return None
        v = np.interp(t_eval[mask_rmse], sol.t, sol["voltage"].entries)
        rmse = np.sqrt(np.mean((v - v_exp[mask_rmse]) ** 2))
        return a1 * time_error + a2 * rmse


    theta = np.log10([inputs[k] for k in names])
    state = evaluate(theta)
    
    theta0 = np.log10([inputs[k] for k in names])

    chain = [theta.copy()]
    history = [state]
    n_acc = 1

    for i in range(n_steps):
        prop = theta + step * rng.standard_normal(len(names))
        out = evaluate(prop)

        if out is not None and out < accept_th:
            theta, state = prop, out
            n_acc += 1

        chain.append(prop.copy())
        history.append(out)
        
        print(f"{names[0]}-{names[1]} Step {i} of {n_steps}, accepted percentage: {n_acc / (i + 2) * 100:.2f}%")
        
    chain = np.array(chain)
    history = np.array(history)
    history[history == None] = np.nan

    np.savetxt(f"data_metropolis/{names[0]}-{names[1]}.data", np.column_stack((chain, history, history < accept_th)), fmt="%.2e",
        header=f"log10 {names[0]}, log10 {names[1]}, a1 * RMSE + a2 * time error, accepted\n"
        "model : DFN, Chen2020\n"
        f"a1 = {a1}, a2 = {a2}, accept_th = {accept_th}\n"
        f"theta0 = {theta0}\n"
        f"accepted percentage = {n_acc / n_steps * 100:.2f}%"
    )
    
    
if __name__ == "__main__":
    from multiprocessing import Pool
    import numpy as np

    args = [
        (("Fk_cathode", "D_s_cathode"), np.array([1, 1])),
        (("Fk_anode", "D_s_anode"), np.array([0.5, 0.5])),
        (("Fk_cathode", "Fk_anode"), np.array([1, 1])),
        (("D_s_cathode", "D_s_anode"), np.array([0.5, 0.5]))
    ]

    with Pool(processes=4) as pool:
        pool.starmap(Metropolis, args)