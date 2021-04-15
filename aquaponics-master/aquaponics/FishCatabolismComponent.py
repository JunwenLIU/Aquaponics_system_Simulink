from .CoreComponent import CoreComponent


class FishCatabolismComponent(CoreComponent):
    """Registers the fish catabolism component of the fish model.

    Model taken from [2] (see Aquaponics.py for sources).

    Parameters
    ----------
    m : GEKKO
        The gekko model.

    Connecting Variables
    --------------------
    T : MV (deg C)
        Receives from environment.
        The temperature over time.
    DO : MV
        Input. Receives from environment.
        Dissolved oxygen.
    DOC : MV (mg DO / l)
        Input. Receives from environment.
        Dissolved oxygen concentration.
    FB : MV (kcal / pond)
        Input. Receives from fish growth.
        Total fish biomass.
    FP : MV (fish / pond)
        Input. Receives from fish growth.
        Total fish population.
    TFC : m.MV (kcal / day / pond)
        Input. Receives from fish anabolism.
        Total food consumption.
    FC : float (kcal / day / pond)
        Output.
        Fish catabolism.
    """

    def __init__(self, m, aqua):
        CoreComponent.__init__(self, m, aqua)

    def register_connecting(self, **kwargs):
        """Registers all connecting variables associated with this component to
        the gekko model.

        kwargs
        ------
        T_0 : float, default = 20
        DO_0 : float, default = 5
        DOC_0 : float, default = 5 (mg DO / l)
        FBi : float, default = 15 * 15
        FP_0 : float, default = 10
        FB_0 : float, default = FBi * FP_0
        FC_0 : float, default = 0
        """
        # Aliases
        a = self.aqua
        MV = self.m.MV
        SV = self.m.SV

        # Initial Conditions
        T_0 = kwargs.get('T_0', 20)
        DO_0 = kwargs.get('DO', 5)
        DOC_0 = kwargs.get('DOC', 5)
        FBi_val = kwargs.get('FBi', 15 * 15)
        FP_0 = kwargs.get('FP_0', 10)
        FB_0 = kwargs.get('FB_0', FBi_val * FP_0)
        TFC_0 = kwargs.get('TFC_0', 0)
        FC_0 = kwargs.get('FC_0', 0)

        # Connecting Variables
        a.register_connecting('T', MV, T_0)
        a.register_connecting('DO', MV, DO_0)
        a.register_connecting('FB', MV, FB_0)
        a.register_connecting('FP', MV, FP_0)
        a.register_connecting('DOC', MV, DOC_0)
        a.register_connecting('TFC', MV, TFC_0)

        a.register_connecting('FC', SV, FC_0)

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
        FBs = m.Param(value=kwargs.get('FBs', 20))
        kfeed = m.Param(value=kwargs.get('kfeed', 0.31))
        kfast = m.Param(value=kwargs.get('kfast', 0.005))
        n = m.Param(value=kwargs.get('n', -0.12))
        c = m.Param(value=kwargs.get('c', 0.59))
        d = m.Param(value=kwargs.get('d', 0.027))
        DOcrit = m.Param(value=kwargs.get('DOcrit', 1.0))
        kDOf = m.Param(value=kwargs.get('kDOf', 2.5))

        # --------------------
        # Connecting Variables
        # --------------------

        T = a.T
        DO = a.DO
        DOC = a.DOC
        FC = a.FC
        TFC = a.TFC
        FB = a.FB
        FP = a.FP

        # ---------------------------
        # Equations and Intermediates
        # ---------------------------

        # Mean fish biomass
        FBm = m.Intermediate(FB / FP)
        # Effects of fish size on fasting catabolism
        fsn = (FBm / FBs) ** n
        # Effect of dissolved oxygen on fasting catabolism
        f2DO = m.switch(
            m.exp(-kDOf * (DOcrit - DOC) ** 2),
            1,
            DO, DOcrit, kswitch
        )
        # Effect of temperature on fasting catabolism
        f3T = m.Intermediate(c + d * T)
        # Catabolism
        m.Equation(FC == TFC * kfeed + FB * kfast * fsn * f2DO * f3T)
