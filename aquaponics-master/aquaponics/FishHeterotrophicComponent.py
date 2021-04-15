from .CoreComponent import CoreComponent


class FishHeterotrophicComponent(CoreComponent):
    """Registers the heterotrophic nutrient sector component of the fish model.

    Model taken from [2] (see Aquaponics.py for sources).

    Parameters
    ----------
    m : GEKKO
        The gekko model.

    Connecting Variables
    --------------------
    DO : MV
        Input. Receives from environment.
        Dissolved oxygen.
    DOC : MV (mg DO / l)
        Input. Receives from environment.
        Dissolved oxygen concentration.
    AFM : MV (kcal / day / pond)
        Input. Receives from autotrophic food pool.
        Autotrophic food entering heterotrophic food pool.
    FW : MV (kcal / day / pond)
        Input. Receives from fish anabolism.
        Fish fecal wastes.
    HFC : MV (kcal / day / pond)
        Input. Receives from fish anabolism.
        Heterotrophic food loss from fish grazing.
    HFe : SV (kcal / pond)
        Ouput.
        Quantity of food nutrients in terms of energy.
    HFp : SV (kcal / pond)
        Output.
        Quantity of food nutrients in terms of protein.
    f1DO : SV
        Output.
        Decomposition of heterotrophic particles.
    HFS : SV
        Output.
        Heterotrophic food loss rate due to sedimentation.
    HFD : SV
        Output.
        Heterotrophic food loss rate due to decomposition.
    """

    def __init__(self, m, aqua):
        CoreComponent.__init__(self, m, aqua)

    def register_connecting(self, **kwargs):
        """Registers all connecting variables associated with this component to
        the gekko model.

        kwargs
        ------
        DO_0 : float, default = 5
        DOC_0 : float, default = 5 (mg DO / l)
        HFe_0 : float, default = 0 (kcal / pond)
        HFp_0 : float, default = 0 (kcal / pond)
        AFM_0 : float, default = 0 (kcal / day / pond)
        FW_0 : float, default = 0 (kcal / day / pond)
        HFC_0 : float, default = 0 (kcal / day / pond)
        """
        # Aliases
        a = self.aqua
        MV = self.m.MV
        SV = self.m.SV

        # Initial Conditions
        DO_0 = kwargs.get('DO', 5)
        DOC_0 = kwargs.get('DOC', 5)
        HFe_0 = kwargs.get('HFe_0', 0)
        HFp_0 = kwargs.get('HFp_0', 0)
        AFM_0 = kwargs.get('AFM_0', 0)
        FW_0 = kwargs.get('FW_0', 0)
        HFC_0 = kwargs.get('HFC_0', 0)

        # Connecting Variables
        a.register_connecting('DO', MV, DO_0)
        a.register_connecting('DOC', MV, DOC_0)
        a.register_connecting('AFM', MV, AFM_0)
        a.register_connecting('FW', MV, FW_0)
        a.register_connecting('HFC', MV, HFC_0)
        a.register_connecting('HFe', SV, HFe_0, lb=0)
        a.register_connecting('HFp', SV, HFp_0, lb=0)
        a.register_connecting('f1DO', SV, 0)
        a.register_connecting('HFS', SV, 0)
        a.register_connecting('HFD', SV, 0)

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
        ks = m.Param(value=kwargs.get('ks', 0.14))
        kd = m.Param(value=kwargs.get('kd', 0.12))
        kp2 = m.Param(value=kwargs.get('kp2', 0.12))
        kDO = m.Param(value=kwargs.get('kDO', 0.14))

        # --------------------
        # Connecting Variables
        # --------------------

        DO = a.DO
        DOC = a.DOC
        AFM = a.AFM
        FW = a.FW
        HFC = a.HFC
        HFe = a.HFe
        HFp = a.HFp
        f1DO = a.f1DO
        HFS = a.HFS
        HFD = a.HFD

        # ---------------------------
        # Equations and Intermediates
        # ---------------------------

        # Decomposition of heterotrophic particles
        m.Equation(
            f1DO == m.switch(m.exp(-kDO * (4 - DOC) ** 2), 1, DO, 4, k=kswitch)
        )
        # Heterotrophic food loss rate due to sedimentation
        m.Equation(HFS == HFe * ks)
        # Heterotrophic food loss rate due to decomposition
        m.Equation(HFD == HFe * kd * f1DO)
        # Heterotrophic food quantity dynamics
        m.Equation(HFe.dt() == AFM + FW - HFC - HFS - HFD)
        m.Equation(HFp.dt() == HFe.dt() * kp2)
