from .CoreComponent import CoreComponent


class FishGrowthComponent(CoreComponent):
    """Registers the fish growth component of the fish model.

    Model taken from [2], mortality taken from [4] (see Aquaponics.py for
    sources).

    Parameters
    ----------
    m : GEKKO
        The gekko model.

    Connecting Variables
    --------------------
    FPs : MV (fish / day / pond)
        Input, user-supplied.
        Fish stocking number.
    FA : MV (kcal / day / pond)
        Input. Receives from fish anabolism.
        Fish anabolism.
    FC : MV (kcal / day / pond)
        Input. Receives from fish catabolism.
        Fish catabolism.
    INC : MV (mg N / l)
        Input. Receives from nitrogen connection.
        Total inorganic nitrogen concentration.
    FB : SV (kcal / pond)
        Output.
        Total fish biomass.
    FP : SV (fish / pond)
        Output.
        Total fish population.
    """

    def __init__(self, m, aqua):
        CoreComponent.__init__(self, m, aqua)

    def register_connecting(self, **kwargs):
        """Registers all connecting variables associated with this component to
        the gekko model.

        kwargs
        ------
        FBi : float, default = 15 * 15
            Average fish biomass at stocking.
        FP_0 : int, default = 10
        FB_0 : float, default = FBi * FP0
        FPs_0 : int, default=0
        FA_0 : float, default=0
        FC_0 : float, default=0
        INC_0 : float, default=1.63
        """
        # Aliases
        a = self.aqua
        MV = self.m.MV
        SV = self.m.Var

        # Initial Conditions
        FBi_val = kwargs.get('FBi', 15 * 15)
        FP_0 = kwargs.get('FP_0', 10)
        FB_0 = kwargs.get('FB_0', FBi_val * FP_0)
        FPs_0 = kwargs.get('FPs_0', 0)
        FA_0 = kwargs.get('FA_0', 0)
        FC_0 = kwargs.get('FC_0', 0)
        INC_0 = kwargs.get('INC_0', 1.63)

        # Connecting Variables
        a.register_connecting('FB', SV, FB_0, lb=0)
        a.register_connecting('FP', SV, FP_0, lb=0)

        a.register_connecting('FPs', MV, FPs_0)
        a.register_connecting('FA', MV, FA_0)
        a.register_connecting('FC', MV, FC_0)
        a.register_connecting('INC', MV, INC_0)

        a.register_connecting('FPdt', SV, value=0)

    def register_equations(self, **kwargs):
        """Registers all equations and intermediates associated with this
        component to the gekko model.

        kwargs
        ------
        FBi : float, default = 15 * 15
            Average fish biomass at stocking.
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

        # Tilapia params from paper
        FBi_val = kwargs.get('FBi', 15 * 15)
        FBi = m.Param(value=FBi_val)
        kNH = m.Param(value=kwargs.get('kNH', 0.365))

        # --------------------
        # Connecting Variables
        # --------------------

        FB = a.FB
        FP = a.FP

        FPs = a.FPs
        FA = a.FA
        FC = a.FC
        INC = a.INC

        FPdt = a.FPdt

        # ---------------------------
        # Equations and Intermediates
        # ---------------------------

        # Fraction of unincorporated ammonia to INC
        UIA = kNH * INC
        # Computation of fish mortality
        km2 = 1 / (1 + m.exp(-6 * (UIA - 1.25)))
        # m.aqua['km2'] = m.Intermediate(km2)
        # Fish biomass and population dynamics
        FBm = FB / FP
        m.Equation(FB.dt() == FPs * FBi + FA - FC - FP * km2 * FBm)
        m.Equation(FPdt == FP.dt())
        m.Equation(FP.dt() == FPs - FP * km2)
