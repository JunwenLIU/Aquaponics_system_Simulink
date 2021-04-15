from .CoreComponent import CoreComponent


class HydroPlantComponent(CoreComponent):
    """Registers the plant component of the hydroponics model. Meant to run
    with the HydroNitrogenComponent.

    Model taken from [3](see Aquaponics.py for sources).

    Parameters
    ----------
    m : GEKKO
        The gekko model.

    Connecting Variables
    --------------------
    """

    def __init__(self, m, aqua):
        CoreComponent.__init__(self, m, aqua)

    def register_connecting(self, **kwargs):
        """Registers all connecting variables associated with this component to
        the gekko model.

        kwargs
        ------
        """
        # Aliases
        a = self.aqua
        CV = self.m.CV
        FV = self.m.FV
        MV = self.m.MV
        SV = self.m.Var

        # Initial Conditions
        T0 = kwargs.get('T0', 25)
        I0 = kwargs.get('I0', 5e6)
        Iwindow = kwargs.get('Iwindow', 1e5)
        Twindow = kwargs.get('Twindow', 5)
        N0 = kwargs.get('N0', 0)

        wG0 = kwargs.get('wG0', 1/3)
        wS0 = kwargs.get('wS0', 2/3)

        K0 = kwargs.get('K', 0.02)
        Jmax_00 = kwargs.get('Jmax_0', 0.0374 * 24)
        alpha0 = kwargs.get('alpha', .151)

        # Connecting Variables
        a.register_connecting(
            'T', MV, T0, lb=(T0 - Twindow), ub=(T0 + Twindow)
        )
        a.register_connecting(
            'I', MV, I0, lb=(I0 - Iwindow), ub=(I0 + Iwindow)
        )
        a.register_connecting('cN', MV, N0)

        a.register_connecting('w', CV, 0)
        a.register_connecting('wG', SV, wG0)
        a.register_connecting('wS', SV, wS0)
        a.register_connecting('Nup', SV, 0, lb=0)
        a.register_connecting('dNup', SV, 0)

        a.register_connecting('K', FV, K0)
        a.register_connecting('Jmax_0', FV, Jmax_00)
        a.register_connecting('alpha', FV, alpha0)

    def register_equations(self, **kwargs):
        """Registers all equations and intermediates associated with this
        component to the gekko model.

        kwargs
        ------
        """
        # -------
        # Aliases
        # -------

        m = self.m
        a = self.aqua

        # ----------
        # Parameters
        # ----------

        # kswitch = kwargs.get('kswitch', 100)
        plant_mumax = m.Param(value=kwargs.get('plant_mumax', 0.01 * 24))
        cq10_mu = m.Param(value=kwargs.get('cq10_mu', 1.6))
        YG = m.Param(value=kwargs.get('YG', 0.8))
        theta = m.Param(value=kwargs.get('theta', 0.68))
        # k = m.Param(value=kwargs.get('plant_k', 0.9))
        xi = m.Param(value=kwargs.get('xi', 14e-6))
        sigma = m.Param(value=kwargs.get('plant_sigma', 7.2 * 24))
        beta = m.Param(value=kwargs.get('plant_beta', 0.36 * 24))
        CCO2 = m.Param(value=kwargs.get('CCO2', 0.8))

        # --------------------
        # Connecting Variables
        # --------------------

        T = a.T
        I = a.I  # noqa
        cN = a.cN

        w = a.w
        wG = a.wG
        wS = a.wS
        Nup = a.Nup
        dNup = a.dNup

        K = a.K
        Jmax_0 = a.Jmax_0
        alpha = a.alpha

        # ---------------------------
        # Equations and Intermediates
        # ---------------------------

        # Dry weight
        m.Equation(w == wG + wS)

        # Dry matter pool
        dWG = m.Intermediate(
            plant_mumax * (wS / (wS + wG)) * wG * cq10_mu ** ((T - 20) / 10)
        )
        # m.Equation(wG.dt() == m.max(dWG, 0, k=kswitch))
        m.Equation(wG.dt() == dWG)
        # Structural pool
        LAIexp = m.switch(0, 1, w, 27.5, k=2)
        fNup = m.switch(0, 1, Nup / w, 0.035, k=150)
        A = m.switch(0, 0.02 * 2.4, w, 1.2, k=2)
        Pg = m.Intermediate(
            A * (1 - LAIexp) * fNup *
            (xi * I * (sigma * CCO2 - beta)) / (xi * I + sigma * CCO2)  # noqa
        )
        dWS = m.Intermediate(theta * Pg - (1 / YG) * dWG)
        m.Equation(wS.dt() == dWS)
        # Nitrogen Uptake
        Jmax = m.Intermediate(Jmax_0 * m.exp(-alpha * w))
        m.Equation(dNup == m.Intermediate(Jmax * cN * w / (cN + K)))
        m.Equation(Nup.dt() == dNup)
