from .CoreComponent import CoreComponent


class BacteriaComponent(CoreComponent):
    """Registers the bacteria component.

    Model taken from [1] (see Aquaponics.py for sources).

    Parameters
    ----------
    m : GEKKO
        The gekko model.

    Connecting Variables
    --------------------
    T : MV (deg C)
        Receives from environment.
        The temperature over time.
    NH3_exc : MV (mg/l)
        Receives from fish (through nitrogen connection).
        The *total* (integral over time) concentration of ammonia produced by
        the fish.
    NH3 : SV (mg/l)
        Output.
        The total concentration of ammonia.
    NO2 : SV (mg/l)
        Output.
        The total concentration of nitrite.
    NO3 : SV (mg/l)
        Output.
        The total concentration of nitrate.
    NO2_up : MV (mg/l)
        Receives from lettuce.
        The *total* (integral over time) amount of nitrite uptake used by the
        plants.
    NO3_up : MV (mg/l)
        Receives from lettuce.
        The *total* (integral over time) amount of nitrate uptake used by the
        plants.
    Cm : Var (mg/l)
        Potentially measurable output.
        The population concentration of nitrosomonas.
    Cb : Var (mg/l)
        Potentially measurable output.
        The population concentration of nitrobacter.
    """

    def __init__(self, m, aqua):
        CoreComponent.__init__(self, m, aqua)

    def register_connecting(self, **kwargs):
        """Registers all connecting variables associated with this component to
        the gekko model.

        kwargs
        ------
        T_0 : float or array, default = 20 (deg C)
            The temperature over time, if a float, the temperature is constant
            across time.
        Cm_0 : float, default = 0.0025 (mg/l)
            The initial population concentration of nitrosomonas.
        Cb_0 : float, default = 0.0005 (mg/l)
            The initial population concentration of nitrobacter.
            specify across time.
        NH3_0 : float, default = 0 (mg/l)
            The initial concentration of ammonia.
        NO2_0 : float, default = 0 (mg/l)
            The initial concentration of nitrite.
        NO3_0 : float, default = 0 (mg/l)
            The initial concentration of nitrate.
            to specify across time.
        """
        # Aliases
        a = self.aqua
        MV = self.m.MV
        SV = self.m.SV

        # Initial Conditions
        T_0 = kwargs.get('T0', 0)
        Cm_0 = kwargs.get('Cm_0', 0.0025)
        Cb_0 = kwargs.get('Cb_0', 0.0005)
        NH3_0 = kwargs.get('NH3_0', 0)
        NO2_0 = kwargs.get('NO2_0', 0)
        NO3_0 = kwargs.get('NO3_0', 0)

        # Connecting Variables
        a.register_connecting('T', MV, T_0)
        a.register_connecting('NH3_exc', MV, 0)
        a.register_connecting('NH3', SV, NH3_0, lb=0)
        a.register_connecting('NO2', SV, NO2_0, lb=0)
        a.register_connecting('NO3', SV, NO3_0, lb=0)
        a.register_connecting('NO2_up', MV, 0)
        a.register_connecting('NO3_up', MV, 0)
        a.register_connecting('Cm', SV, Cm_0, lb=0)
        a.register_connecting('Cb', SV, Cb_0, lb=0)

    def register_equations(self, **kwargs):
        """Registers all equations and intermediates associated with this
        component to the gekko model.

        kwargs
        ------
        Cm_0 : float, default = 0.0025 (mg/l)
            The initial population concentration of nitrosomonas.
        Cb_0 : float, default = 0.0005 (mg/l)
            The initial population concentration of nitrobacter.
            specify across time.
        NH3_0 : float, default = 0 (mg/l)
            The initial concentration of ammonia.
        NO2_0 : float, default = 0 (mg/l)
            The initial concentration of nitrite.
        NO3_0 : float, default = 0 (mg/l)
            The initial concentration of nitrate.
            to specify across time.
        Em : float, default = 0.05 (unitless)
            Production rate of dry mass of nitrosomonas from ammonia.
        Eb : float, default = 0.02 (unitless)
            Production rate of dry mass of nitrobacter from nitrite.
        fm : float in [0,1], default = 0.99
            The ratio of the mass of nitrite formed to that of ammonia
            oxidized.
        fn : float in [0, 1], default = 0.99
            The ratio of the mass of nitrate formed to that of nitrite
            oxidized.
        """
        # -------
        # Aliases
        # -------

        m = self.m
        a = self.aqua

        # ----------
        # Parameters
        # ----------

        Em = m.Param(value=kwargs.get('Em', 0.05))
        Eb = m.Param(value=kwargs.get('Eb', 0.02))
        fm = m.Param(value=kwargs.get('fm', 0.99))
        fn = m.Param(value=kwargs.get('fn', 0.99))

        # Need initial conditions for equations too
        Cm_0 = kwargs.get('Cm_0', 0.0025)
        Cb_0 = kwargs.get('Cb_0', 0.0005)
        NH3_0 = kwargs.get('NH3_0', 0)
        NO2_0 = kwargs.get('NO2_0', 0)
        NO3_0 = kwargs.get('NO3_0', 0)

        # --------------------
        # Connecting Variables
        # --------------------

        T = a.T
        NH3_exc = a.NH3_exc
        NH3 = a.NH3
        NO2 = a.NO2
        NO3 = a.NO3
        NO2_up = a.NO2_up
        NO3_up = a.NO3_up
        Cm = a.Cm
        Cb = a.Cb

        # ---------------------------
        # Equations and Intermediates
        # ---------------------------

        # Growth Constants
        km = 10 ** (0.0413 * T - 0.944)
        kb = 10 ** (0.0255 * T - 0.492)
        # Saturation constants
        NH3_m = 10 ** (0.051 * T - 1.158)
        NO2_m = 10 ** (0.063 * T - 1.149)
        # Growth of Nitrosomonas
        m.Equation(Cm.dt() == km * Cm * NH3 / (NH3 + NH3_m))
        # Growth of Nitrobacter
        m.Equation(Cb.dt() == kb * Cb * NO2 / (NO2 + NO2_m))
        # Oxidation of ammonia
        m.Equation(Cm - Cm_0 == Em * (NH3_0 + NH3_exc - NH3))
        # m.Equation(NH3.dt() == -Cm.dt() / Em + NH3_exc.dt())
        # # Oxidation of nitrite
        m.Equation(
            Cb - Cb_0 == Eb * (NO2_0 - NO2_up + fm * (NH3_0 + NH3_exc - NH3) -
                               NO2)
        )
        # m.Equation(NO2.dt() == -Cb.dt() / Eb - NO2_up.dt() - fm * NH3.dt())
        # # Balance of nitrate
        # m.Equation(NO3 == NO3_0 + fm * (NH3_0 + NH3_exc - NH3) +
        #            fn * (NO2_0 - NO2_up - NO2) - NO3_up)
        m.Equation(NO3 == NO3_0 - NO3_up + fn * (
            NO2_0 - NO2_up + fm * (NH3_0 + NH3_exc - NH3) - NO2)
        )
