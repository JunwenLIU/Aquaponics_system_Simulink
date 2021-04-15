from .CoreComponent import CoreComponent


class HydroNitrogenComponent(CoreComponent):
    """Registers the nitrogen component of the hydroponics model. Meant to run
    with the HydroPlantComponent.

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
        a = self.aqua
        SV = self.m.SV
        MV = self.m.MV

        N0 = kwargs.get('N0', 0)
        dNadd = kwargs.get('dNadd', 0)

        a.register_connecting('cN', SV, N0, lb=0)
        a.register_connecting('dNup', MV, 0)
        a.register_connecting('dNadd', MV, dNadd, lb=0)
        a.register_connecting('Nadd', SV, value=0)

    def register_equations(self, **kwargs):
        """Registers all equations and intermediates associated with this
        component to the gekko model.

        kwargs
        ------
        """
        m = self.m
        a = self.aqua

        cN = a.cN
        dNup = a.dNup
        dNadd = a.dNadd
        Nadd = a.Nadd

        m.Equation(cN.dt() == -dNup + dNadd)
        m.Equation(Nadd.dt() == dNadd)
