
class PlantBed():
    """Internal to the plant component. Defines a single bed.
    """

    def register_equations(self, m, planting_day, harvest_day, t, T, I, N,
                           **kwargs):
        """Registers the equations for a single grow bed into m.

        Parameters
        ----------
        m : GEKKO
            The gekko model.
        planting_day : number (days)
            The day on which the bed is planted.
        harvest_day : number (days)
            The day on which the bed is harvested.
        t : gekko variable (days)
            The current simulation time.
        T : gekko variable (deg C)
            Temperature.
        I : gekko variable (J / m^2 / day)
            Solar radiation.
        N : gekko variable (mM)
            Nitrogen available.
        """

        # ----------
        # Parameters
        # ----------

        kswitch = kwargs.get('kswitch', 100)

        # Nitrogen uptake
        K = m.Param(value=kwargs.get('K', 0.02))
        Jmax_0 = m.Param(value=kwargs.get('Jmax_0', 0.0374*24))
        alpha = m.Param(value=kwargs.get('alpha', .151))
        plant_mumax = m.Param(value=kwargs.get('plant_mumax', 0.01 * 24))
        cq10_mu = m.Param(value=kwargs.get('cq10_mu', 1.6))
        YG = m.Param(value=kwargs.get('YG', 0.8))
        theta = m.Param(value=kwargs.get('theta', 0.68))
        # k = m.Param(value=kwargs.get('plant_k', 0.9))
        xi = m.Param(value=kwargs.get('xi', 14e-6))
        sigma = m.Param(value=kwargs.get('plant_sigma', 7.2 * 24))
        beta = m.Param(value=kwargs.get('plant_beta', 0.36 * 24))
        CCO2 = m.Param(value=kwargs.get('CCO2', 0.8))

        wG0 = kwargs.get('wG0', 1/3)
        wS0 = kwargs.get('wS0', 2/3)

        # ------------------------
        # Internal State Variables
        # ------------------------

        Nup = m.Var(value=0)
        wG = m.Var(value=wG0)
        wS = m.Var(value=wS0)
        w = m.Var(value=0)
        dNup = m.Var(value=0)

        # ---------------------------
        # Equations and Intermediates
        # ---------------------------

        cN = N

        on_harvest = m.switch(1, 0, t, harvest_day, kswitch)
        on = m.switch(0, on_harvest, t, planting_day, kswitch)

        # Dry weight
        m.Equation(w == wG + wS)

        # Dry matter pool
        dWG = m.Intermediate(
            plant_mumax * (wS / (wS + wG)) * wG * cq10_mu ** ((T - 20) / 10)
        )
        # m.Equation(wG.dt() == m.max(dWG, 0, k=kswitch))
        m.Equation(wG.dt() == dWG * on)
        # Structural pool
        LAIexp = m.switch(0, 1, w, 27.5, k=2)
        fNup = m.switch(0, 1, Nup / w, 0.035, k=150)
        A = m.switch(0, 0.02 * 2.4, w, 1.2, k=2)
        Pg = m.Intermediate(
            A * (1 - LAIexp) * fNup *
            (xi * I * (sigma * CCO2 - beta)) / (xi * I + sigma * CCO2)  # noqa
        )
        dWS = m.Intermediate(theta * Pg - (1 / YG) * dWG)
        m.Equation(wS.dt() == dWS * on)
        # Nitrogen Uptake
        Jmax = m.Intermediate(Jmax_0 * m.exp(-alpha * w))
        m.Equation(dNup == m.Intermediate(Jmax * cN * w / (cN + K)) * on)
        m.Equation(Nup.dt() == dNup)

        return w, dNup
