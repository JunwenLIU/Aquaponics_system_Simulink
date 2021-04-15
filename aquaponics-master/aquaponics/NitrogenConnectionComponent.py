from math import pi
from .CoreComponent import CoreComponent
from .FishNutrientComponent import FishNutrientComponent as nutrient
from .BacteriaComponent import BacteriaComponent as bacteria
from .PlantComponent import PlantComponent as plant


class NitrogenConnectionComponent(CoreComponent):
    """Connects the nitrogen used in each of the three models.

    Connecting Variables
    --------------------
    dw : FV (cm)
        User input.
        The fish tank depth.
    rw : FV (cm)
        User input.
        The fish tank radius.
    INC : SV (mg N / l)
        Output.
        Total inorganic nitrogen concentration.
    IPC : SV (mg P / l)
        Output.
        Total inorganic phosphorus concentration.

    Connecting Variables if Fish Nutrient Component Exists
    ------------------------------------------------------

    Connecting Variables if Bacteria Component Exists
    -------------------------------------------------
    NH3_exc : SV (mg/l)
        Output.
        The *total* (integral over time) concentration of ammonia produced by
        the fish.
    NH3 : MV (mg/l)
        Input. Receives from bacteria.
        The total concentration of ammonia.
    NO2 : MV (mg/l)
        Input. Receives from bacteria.
        The total concentration of nitrite.
    NO3 : MV (mg/l)
        Input. Receives from bacteria.
        The total concentration of nitrate.
    """

    def __init__(self, m, aqua):
        CoreComponent.__init__(self, m, aqua)

    def register_connecting(self, **kwargs):
        """Registers all connecting variables associated with this component to
        the gekko model.

        kwargs
        ------
        INC_0 : float, default = 0
        IPC_0 : float, default = 0
        """
        # Aliases
        a = self.aqua
        # MV = self.m.MV
        SV = self.m.SV
        FV = self.m.FV

        # Initial Conditions
        dw = kwargs.get('dw', 40)
        rw = kwargs.get('dw', 120)
        INC_0 = kwargs.get('INC_0', 0)
        IPC_0 = kwargs.get('IPC_0', 0)
        NH3_exc_0 = kwargs.get('NH3_exc_0', 0)
        NO3_up_0 = kwargs.get('NO3_up_0', 0)
        N0 = kwargs.get('N0', 0)

        # Tank dimensions
        a.register_connecting('dw', FV, dw)
        a.register_connecting('rw', FV, rw)

        # Nutrient variables
        if nutrient in a.registered_components:
            a.register_connecting('INC', SV, INC_0, lb=0)
            a.register_connecting('IPC', SV, IPC_0, lb=0)

        # Upgrade NH3_exc to SV if both nutrient and bacteria are registered
        if bacteria in a.registered_components and \
                nutrient in a.registered_components:
            a.register_connecting('NH3_exc', SV, NH3_exc_0)

        if plant in a.registered_components and \
                bacteria in a.registered_components:
            a.register_connecting('N', SV, N0)
            a.register_connecting('NO3_up', SV, NO3_up_0)

    def register_equations(self, **kwargs):
        """Registers all equations and intermediates associated with this
        component to the gekko model.

        kwargs
        ------
        (coming soon, see paper for model parameters)
        """
        # Aliases
        a = self.aqua
        m = self.m

        # Common connecting variables
        dw = a.dw
        rw = a.rw

        # Tank Volume (l / pond)
        vw = m.Intermediate(pi * rw ** 2 * dw * 0.001)

        # Equations
        if bacteria in a.registered_components and \
                nutrient in a.registered_components:
            INC = a.INC
            IPC = a.IPC
            TIN = a.TIN
            TIP = a.TIP
            NH3 = a.NH3
            NH3_exc = a.NH3_exc
            NO2 = a.NO2
            NO3 = a.NO3

            # Change in nitrogen (mg N / l / day)
            dINC = m.Intermediate(TIN.dt() * 1000 / vw)
            m.Equation(NH3_exc.dt() == dINC)
            m.Equation(INC == NH3 + NO2 + NO3)

            # Change in phosphorus (mg P / l / day)
            dIPC = m.Intermediate(TIP.dt() * 1000 / vw)
            m.Equation(IPC.dt() == dIPC)

        elif nutrient in a.registered_components:
            # Nutrient connecting variables
            INC = a.INC
            IPC = a.IPC
            TIN = a.TIN
            TIP = a.TIP

            # Change in nitrogen (mg N / l / day)
            dINC = m.Intermediate(TIN.dt() * 1000 / vw)
            m.Equation(INC.dt() == dINC)

            # Change in phosphorus (mg P / l / day)
            dIPC = m.Intermediate(TIP.dt() * 1000 / vw)
            m.Equation(IPC.dt() == dIPC)

        if plant in a.registered_components and \
                bacteria in a.registered_components:
            mm = 62.0049  # molar mass NO3 (g / mol)
            NO3 = a.NO3
            NO3_up = a.NO3_up
            N = a.N
            dNup = a.dNup

            # Conversion of NO3 to N
            m.Equation(N == NO3 / mm)

            # Conversion of dNup (mmol / day) to NO3_up (mg / l / day)
            dNO3_up = m.Intermediate(dNup * mm / vw)
            m.Equation(NO3_up.dt() == dNO3_up)
