from .CoreComponent import CoreComponent
from .PlantBed import PlantBed


class PlantComponent(CoreComponent):
    """Registers the Plant component.

    Model taken from [3] (see Aquaponics.py for sources).

    Parameters
    ----------
    m : GEKKO
        The gekko model.

    Available kwargs
    ----------------
    T_0 : float or array, default = 20 (deg C)
        The temperature over time, if a float, the temperature is constant
        across time.
    Nup_0 : float, default = ? (mmol)
        The initial amount of nitrogen in the shoot
    J_max_0 : float, default = 0.0423 (mmol/(h g of DM))
        Maximum nitrogen uptake rate per unit of shoot dry matter.
    w : float, default = ? (g)
        Plant shoot dry weight.
    K : float, default = 0.0359 (mM)
        The semi-saturation constant for nitrogen uptake.
    alpha : float, default = 0.154 (1/g of DM)
        Adjustment to J_max for decreased nitrogen demand over plant life.
    mu_max : float, default = ? (hr-1)
        Saturation growth rate at 20 deg C.
    C_Q10 : float, default = ? ()
        Q10 factor for growth.
    theta : float, default = ? ()
        Factor to convert CO_2 to dry matter.
    Y_G : float, default = ? ()
        Conversion (CO_2 to dry matter) efficiency
    A : float, default = ? (m^2)
        The ground cover area per plant.
    k_c : float, default = ?
        The extinction Coefficient
    L_AI : float, default = ? ()
        Leaf area index
    coeff : variable, default = 0
        coefficient betwee 0 to 1 controlled by nitrogen concentration in shoot
    psi : float, default = ? (g(CO_2)/J)
        leaf light use efficiency
    Beta : float, default = ? (g(CO_2)/(m2 hr))
        CO_2 compensation point to account for photrespiration
    I : float, default = ? (J/(m2 hr))
        Incident photosynthetically active radiation
    C_CO2 : float, default = ? (g(CO_2)/m3)
        CO_2 concentration in the air
    sigma : float, default = ? (m/hr)
        Leaf conductance to CO_2 diffusion.
    w_g : SV (g)
        Potential measurable variable. Dry weight of the structural pool.
    w_s : SV (g)
        Potential measurable variable. Dry weight of the non-structural pool.
    P_g : SV (g_CO2/h)
        Gross canopy photosynthesis rate.

    Connecting Gekko Variables Registered
    -------------------------------------
    T : MV (deg C)
        Input. Receives from environment.
        The temperature over time.
    I : MV (J / m^2 / day)
        Input. Receives from environment.
        Solar radiation over time.
    N : MV (mM)
        Input. Receives from nitrogen connection.
        Available nitrogen for uptake.
    ppb : FV (int, unitless)
        User input.
        Plants per bed.
    w : SV (g)
        Output.
        Total plant dry weight across all plants and beds.
    y : SV (kg)
        Output.
        Total crop yield across all plants and beds.
    dNup : SV (mmol / day / bed)
        Output.
        The change in nitrogen uptake.
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
        """
        # Aliases
        a = self.aqua
        MV = self.m.MV
        SV = self.m.Var
        FV = self.m.FV

        # Initial Conditions
        T0 = kwargs.get('T0', 25)
        I0 = kwargs.get('I0', 5e6)
        N0 = kwargs.get('N0', 0)
        ppb0 = kwargs.get('ppb0', 1)
        # beds = kwargs.get('beds', [(0, 30)])
        # wG0 = kwargs.get('wG0', 1/3)
        # wS0 = kwargs.get('wS0', 2/3)
        # w0 = ppb0 * len(beds) * (wG0 + wS0)
        w0 = 0  # weird gekko bug, instability if not set to 0

        # Register connecting
        a.register_connecting('T', MV, T0)
        a.register_connecting('I', MV, I0)
        a.register_connecting('N', MV, N0)
        a.register_connecting('ppb', FV, ppb0)
        # a.ppb = self.m.Param(value=ppb0)

        a.register_connecting('w', SV, w0, lb=0)
        a.register_connecting('dNup', SV, 0, lb=0)

    def register_equations(self, **kwargs):
        """Registers all equations and intermediates associated with this
        component to the gekko model.

        kwargs
        ------
        beds : list of 2-tuple
            List of (planting_day, harvest_day) defining each growing bed.
        NH3_0 : float, default = 0 (mg/l)
            The initial concentration of ammonia.
        NO2_0 : float, default = 0 (mg/l)
            The initial concentration of nitrite.
        NO3_0 : float, default = 0 (mg/l)
            The initial concentration of nitrate.
            to specify across time.
        Nup_0 : float, default = ? (mmol)
            The initial amount of nitrogen in the shoot
        J_max_0 : float, default = ? (mmol/(h g of DM))
            Maximum nitrogen uptake rate per unit of shoot dry matter.
        w : float, default = ? (g)
            Plant shoot dry weight.
        K : float, default = ? (mM)
            The semi-saturation constant for nitrogen uptake.
        alpha : float, default = ? (1/g of DM)
            Adjustment to J_max for decreased nitrogen demand over plant life.
        """
        # -------
        # Aliases
        # -------

        m = self.m
        a = self.aqua

        # ----------
        # Parameters
        # ----------

        # Growing bed definition
        beds = kwargs.get('beds', [(0, 30)])

        # --------------------
        # Connecting Variables
        # --------------------

        T = a.T
        I = a.I   # noqa
        N = a.N
        ppb = a.ppb
        w = a.w
        dNup = a.dNup

        # ---------------------------
        # Equations and Intermediates
        # ---------------------------

        time = m.SV(value=0)
        m.timevar = time
        m.Equation(time.dt() == 1)

        bed_models = [
            (PlantBed(), plant_day, harvest_day)
            for plant_day, harvest_day in beds
        ]
        bed_vars = [
            bed.register_equations(
                m, plant_day, harvest_day, time, T, I, N, **kwargs
            )
            for bed, plant_day, harvest_day in bed_models
        ]

        m.Equation(w == ppb * sum([var[0] for var in bed_vars]))
        m.Equation(dNup == ppb * sum([var[1] for var in bed_vars]))
