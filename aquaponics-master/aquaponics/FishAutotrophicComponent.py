from .CoreComponent import CoreComponent


class FishAutotrophicComponent(CoreComponent):
    """Registers the autotrophic nutrient sector component of the fish model.

    Model taken from [2] (see Aquaponics.py for sources).

    Parameters
    ----------
    m : GEKKO
        The gekko model.

    Connecting Variables
    --------------------
    T : m.MV (deg C)
        Input. Receives from environment.
        The temperature over time.
    I0 : m.MV (10^6 cal / m^2 / day)
        Input. Receives from environment.
        The light reaching the surface of the water.
    INC : m.MV (mg N / l)
        Input. Receives from nitrogen connection.
        Total inorganic nitrogen concentration.
    IPC : m.MV (mg P / l)
        Input.
        Total inorganic phosphorus concentration.
    HFe : m.MV (kcal / pond)
        Input. Receives from heterotrophic component.
        Quantity of food nutrients in terms of energy.
    AFC : m.MV (kcal / day / pond)
        Input. Receives from fish anabolism.
        Autotrophic food loss rate due to tilapia grazing.
    AFe : m.SV (kcal / pond)
        Output.
        The autotrophic food quantity in terms of energy.
    AFp : m.SV (g protein / pond)
        Output.
        The autotrophic food quantity in terms of protein.
    AFM : m.SV (kcal / day / pond)
        Output.
        Autotrophic food entering heterotrophic food pool.
    AFR : m.SV (kcal / day / pond)
        Output.
        Autotrophic food loss due to phytoplankton respiration.
    AFR : m.SV (kcal / day / pond)
        Output.
        Autotrophic food loss due to phytoplankton growth.
    """

    def __init__(self, m, aqua):
        CoreComponent.__init__(self, m, aqua)

    def register_connecting(self, **kwargs):
        """Registers all connecting variables associated with this component to
        the gekko model.

        kwargs
        ------
        AFe_0 : float, default = 0 (kcal / pond)
        AFp_0 : float, default = 0 (g protein / pond)
        AFC_0 : float, default = 0 (kcal / day / pond)
        HFe_0 : float, default = 0 (kcal / pond)
        INC_0 : float, default = 1.63 (mg N / l)
        IPC_0 : float, default = 16.86 (mg P / l)
        T_0 : float, default = 20
        I0_0 : float, default = 2.5 / 4.184 (10^6 cal / m^2 / day)
        """
        # Aliases
        a = self.aqua
        MV = self.m.MV
        SV = self.m.SV

        # Initial Conditions
        T_0 = kwargs.get('T_0', 20)
        I0_0 = kwargs.get('I0_0', 2.5 / 4.184)

        AFe_0 = kwargs.get('AFe_0', 0)
        AFp_0 = kwargs.get('AFp_0', 0)
        AFC_0 = kwargs.get('AFC_0', 0)
        HFe_0 = kwargs.get('HFe_0', 0)

        INC_0 = kwargs.get('INC_0', 1.63)
        IPC_0 = kwargs.get('IPC_0', 16.86)

        # Connecting Variables
        a.register_connecting('T', MV, T_0)
        a.register_connecting('I0', MV, I0_0)
        a.register_connecting('AFC', MV, AFC_0)
        a.register_connecting('HFe', MV, HFe_0, lb=0)
        a.register_connecting('INC', MV, INC_0)
        a.register_connecting('IPC', MV, IPC_0)

        a.register_connecting('AFe', SV, AFe_0, lb=0)
        a.register_connecting('AFp', SV, AFp_0, lb=0)
        a.register_connecting('AFM', SV, 0)
        a.register_connecting('AFR', SV, 0)
        a.register_connecting('AFG', SV, 0)

    def register_equations(self, **kwargs):
        """Registers all equations and intermediates associated with this
        component to the gekko model.

        kwargs
        ------
        kswitch : int > 0
            The quality of the switch/max/min functions. Larger is better.
        (coming soon, see paper for model parameters)
        """
        # -------
        # Aliases
        # -------

        m = self.m
        a = self.aqua

        # ----------
        # Parameters
        # ----------

        kswitch = kwargs.get('kswitch', 100)

        # Tilapia params from paper
        kp1 = m.Param(value=kwargs.get('kp1', 0.14))
        mumax = m.Param(value=kwargs.get('mumax', 1.6))
        km1 = m.Param(value=kwargs.get('km1', 0.6))
        kr = m.Param(value=kwargs.get('kr', 0.1))
        hN = m.Param(value=kwargs.get('hN', 0.2))
        hP = m.Param(value=kwargs.get('hP', 0.02))
        _a = m.Param(value=kwargs.get('a', 0.000017))
        b = m.Param(value=kwargs.get('b', 0.000015))
        kT1 = m.Param(value=kwargs.get('kT1', 0.004))
        kT2 = m.Param(value=kwargs.get('kT2', 0.008))

        Ir = m.Param(value=kwargs.get('Ir', 6.547))
        Topta = m.Param(value=kwargs.get('Topta', 30))

        # --------------------
        # Connecting Variables
        # --------------------

        T = a.T
        I0 = a.I0
        AFC = a.AFC
        HFe = a.HFe
        INC = a.INC
        IPC = a.IPC

        AFe = a.AFe
        AFp = a.AFp
        AFM = a.AFM
        AFR = a.AFR
        AFG = a.AFG

        # ---------------------------
        # Equations and Intermediates
        # ---------------------------

        # Elementary nutrient limitation
        fNP = m.min(INC / (INC + hN), IPC / (IPC + hP), k=kswitch)
        # Light limitation
        fI = I0 * m.exp(-(_a * AFe + b * HFe)) / Ir
        # Phytoplankton growth vs. temperature relationship
        f1T = m.switch(
            m.exp(-kT1 * (T - Topta) ** 2),
            m.exp(-kT2 * (Topta - T) ** 2),
            T, Topta, k=kswitch
        )

        # Autotrophic food growth due to phytoplankton growth
        #   (kcal / day / pond)
        m.Equation(AFG == mumax * AFe * fNP * fI * f1T)
        # Autotrophic food loss rate due to phytoplankton respiration
        #   (kcal / day / pond)
        m.Equation(AFR == AFe * kr)
        # Autotrophic food entering heterotrophic food pool (kcal / day / pond)
        m.Equation(AFM == AFe * km1)

        # Autotrophic food quantity dynamics
        m.Equation(AFe.dt() == AFG - AFC - AFR - AFM)
        m.Equation(AFp.dt() == AFe.dt() * kp1)
