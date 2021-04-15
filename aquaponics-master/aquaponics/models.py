"""OBSOLETE!!!

Collection of utility functions for building components of an aquaponics model.

Each component takes a GEKKO model `m` as an imput, a long with some key word
arguments (typically used to set initial conditions), and returns that model
with the parameters, variables, equations, etc. for that component added.
Variables are registered internal to the model and can be accessed through

        m.aqua['VARIABLE NAME']

Each component is designed to be independent from other components. If a
model is registered with a single component, then it can run that single
component as a stand-alone model. If a model is registered with multiple
components, then the variables which are specified as connecting variables
in the components' documentations will be used to build the interconnection of
components, allowing the combined system to be modeled.

A model can register components in any order. However, initial conditions
for any connecting variable must be specified in the first component that
uses it, or on the SV if one of the components creates the variable as an SV;
otherwise, it will be ignored.

Note that components DO NOT specify the model's time or options, and do not set
the FSTATUS, STATUS, etc. of variables. This should be done outside in order
to allow the model to be a simulator, estimator, controller, etc.

Also note that if the fish module is to be used, then the solver should be set
to 1 in order to get integer population counts. However, solver=3 can be used
for a relaxed (floating point) population count.

Example Usage
-------------

        (coming soon)


Sources
-------
[1] Knowles, G., Downing, A. L., & Barrett, M. J. (1965). Determination of
    kinetic constants for nitrifying bacteria in mixed culture, with the aid of
    an electronic computer. Microbiology, 38(2), 263-278.
[2] Li, L., & Yakupitiyage, A. (2003). A model for food nutrient dynamics of
    semi-intensive pond fish culture. Aquacultural Engineering, 27(1), 9-38.
[3] Zhang, K., Burns, I. G., & Turner, M. K. (2008). Derivation of a dynamic
    model of the kinetics of nitrogen uptake throughout the growth of lettuce:
    calibration and validation. Journal of plant nutrition, 31(8), 1440-1460.
[4] Evans, J. J., Pasnik, D. J., Brill, G. C., & Klesius, P. H. (2006).
    Un‐ionized Ammonia Exposure in Nile Tilapia: Toxicity, Stress Response, and
    Susceptibility to Streptococcus agalactiae. North American Journal of
    Aquaculture, 68(1), 23-33.
"""


###############################################################################
#   Model Components (Public Functions)
###############################################################################

def register_lettuce(m, **kwargs):
    return m


def register_fish_autotrophic(m, **kwargs):
    """The autotrophic nutrient sector component of the fish model from [2].

    Parameters
    ----------
    m : GEKKO
        The gekko model.

    Available kwargs
    ----------------
    kswitch : int > 0
        The quality of the switch/max/min functions. Larger is better.
    AFe_0 : float, default = 0 (kcal / pond)
        Initial AFe.
    AFp_0 : float, default = 0 (g protein / pond)
        Initial AFp.
    AFC_0 : float, default = 0 (kcal / day / pond)
        Initial AFC.
    HFe_0 : float, default = 0 (kcal / pond)
        Initial HFe.
    INC_0 : float, default = 1.63 (mg N / l)
        Initial INC.
    IPC_0 : float, default = 16.86 (mg P / l)
        Initial IPC.
    T_0 : float, default = 20
        The initial temperature.
    I0_0 : float, default = 2.5 / 4.184 (10^6 cal / m^2 / day)
        The initial I0.

    Fish-specific Available kwargs
    ------------------------------
    (coming soon, see paper)

    Connecting Gekko Variables Registered
    -------------------------------------
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

    Returns
    -------
    m : GEKKO
        The gekko model with the fish component added.
    """

    # ------------------
    # Initial Conditions
    # ------------------

    T_0 = kwargs.get('T_0', 20)
    I0_0 = kwargs.get('I0_0', 2.5 / 4.184)

    AFe_0 = kwargs.get('AFe_0', 0)
    AFp_0 = kwargs.get('AFp_0', 0)
    AFC_0 = kwargs.get('AFC_0', 0)
    HFe_0 = kwargs.get('HFe_0', 0)

    INC_0 = kwargs.get('INC_0', 1.63)
    IPC_0 = kwargs.get('IPC_0', 16.86)

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
    a = m.Param(value=kwargs.get('a', 0.000017))
    b = m.Param(value=kwargs.get('b', 0.000015))
    kT1 = m.Param(value=kwargs.get('kT1', 0.004))
    kT2 = m.Param(value=kwargs.get('kT2', 0.008))

    Ir = m.Param(value=kwargs.get('Ir', 6.547))
    Topta = m.Param(value=kwargs.get('Topta', 30))

    # --------------------
    # Connecting Variables
    # --------------------

    T = _register_connecting(m, 'T', m.MV, T_0)
    I0 = _register_connecting(m, 'I0', m.MV, I0_0)
    AFC = _register_connecting(m, 'AFC', m.MV, AFC_0)
    HFe = _register_connecting(m, 'HFe', m.MV, HFe_0, lb=0)
    INC = _register_connecting(m, 'INC', m.MV, INC_0)
    IPC = _register_connecting(m, 'IPC', m.MV, IPC_0)

    AFe = _register_connecting(m, 'AFe', m.SV, AFe_0, lb=0)
    AFp = _register_connecting(m, 'AFp', m.SV, AFp_0, lb=0)
    AFM = _register_connecting(m, 'AFM', m.SV, 0)

    # ---------------------------
    # Equations and Intermediates
    # ---------------------------

    # Elementary nutrient limitation
    fNP = m.min(INC / (INC + hN), IPC / (IPC + hP), k=kswitch)
    # Light limitation
    fI = I0 * m.exp(-(a * AFe + b * HFe)) / Ir
    # Phytoplankton growth vs. temperature relationship
    f1T = m.switch(
        m.exp(-kT1 * (T - Topta) ** 2),
        m.exp(-kT2 * (Topta - T) ** 2),
        T, Topta, k=kswitch
    )

    # Autotrophic food growth due to phytoplankton growth (kcal / day / pond)
    AFG = mumax * AFe * fNP * fI * f1T
    # Autotrophic food loss rate due to phytoplankton respiration
    #   (kcal / day / pond)
    AFR = AFe * kr
    # Autotrophic food entering heterotrophic food pool (kcal / day / pond)
    m.Equation(AFM == AFe * km1)

    # Autotrophic food quantity dynamics
    m.Equation(AFe.dt() == AFG - AFC - AFR - AFM)
    m.Equation(AFp.dt() == AFe.dt() * kp1)

    return m


def register_fish_heterotrophic(m, **kwargs):
    """The heterotrophic nutrient sector component of the fish model from [2].

    Parameters
    ----------
    m : GEKKO
        The gekko model.

    Available kwargs
    ----------------
    kswitch : int > 0
        The quality of the switch/max/min functions. Larger is better.
    DO_0 : float, default = 5
    DOC_0 : float, default = 5 (mg DO / l)
    HFe_0 : float, default = 0 (kcal / pond)
    HFp_0 : float, default = 0 (kcal / pond)
    AFM_0 : float, default = 0 (kcal / day / pond)
    FW_0 : float, default = 0 (kcal / day / pond)
    HFC_0 : float, default = 0 (kcal / day / pond)

    Fish-specific Available kwargs
    ------------------------------
    (coming soon, see paper)

    Connecting Gekko Variables Registered
    -------------------------------------
    DO : m.MV
        Input. Receives from environment.
        Dissolved oxygen.
    DOC : m.MV (mg DO / l)
        Input. Receives from environment.
        Dissolved oxygen concentration.
    AFM : m.MV (kcal / day / pond)
        Input. Receives from autotrophic food pool.
        Autotrophic food entering heterotrophic food pool.
    FW : m.MV (kcal / day / pond)
        Input. Receives from fish anabolism.
        Fish fecal wastes.
    HFC : m.MV (kcal / day / pond)
        Input. Receives from fish anabolism.
        Heterotrophic food loss from fish grazing.
    HFe : m.SV (kcal / pond)
        Ouput.
        Quantity of food nutrients in terms of energy.
    HFp : m.SV (kcal / pond)
        Output.
        Quantity of food nutrients in terms of protein.

    Returns
    -------
    m : GEKKO
        The gekko model with the fish component added.
    """

    # ------------------
    # Initial Conditions
    # ------------------

    DO_0 = kwargs.get('DO', 5)
    DOC_0 = kwargs.get('DOC', 5)
    HFe_0 = kwargs.get('HFe_0', 0)
    HFp_0 = kwargs.get('HFp_0', 0)
    AFM_0 = kwargs.get('AFM_0', 0)
    FW_0 = kwargs.get('FW_0', 0)
    HFC_0 = kwargs.get('HFC_0', 0)

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

    DO = _register_connecting(m, 'DO', m.MV, DO_0)
    DOC = _register_connecting(m, 'DO', m.MV, DOC_0)
    AFM = _register_connecting(m, 'AFM', m.MV, AFM_0)
    FW = _register_connecting(m, 'FW', m.MV, FW_0)
    HFC = _register_connecting(m, 'HFC', m.MV, HFC_0)
    HFe = _register_connecting(m, 'HFe', m.SV, HFe_0, lb=0)
    HFp = _register_connecting(m, 'HFp', m.SV, HFp_0, lb=0)

    # ---------------------------
    # Equations and Intermediates
    # ---------------------------

    # Decomposition of heterotrophic particles
    f1DO = m.switch(m.exp(-kDO * (4 - DOC) ** 2), 1, DO, 4, k=kswitch)
    # Heterotrophic food loss rate due to sedimentation
    HFS = HFe * ks
    # Heterotrophic food loss rate due to decomposition
    HFD = HFe * kd * f1DO
    # Heterotrophic food quantity dynamics
    m.Equation(HFe.dt() == AFM + FW - HFC - HFS - HFD)
    m.Equation(HFp.dt() == HFe.dt() * kp2)

    return m


def register_fish_growth(m, **kwargs):
    """The heterotrophic nutrient sector component of the fish model from [2].

    The fish mortality coefficient is computed using [4]. Data shows 0%
    mortality per day with UIA < 0.5, 100% mortality per day with UIA > 2, and
    approximately 10-20% mortality per day with UIA = 1. We fit a logistics
    curve of the form 1 / (1 + exp(-k(x - x0))), with x = UIA, x0 = 1.25, and
    k = 6.

    Parameters
    ----------
    m : GEKKO
        The gekko model.

    Available kwargs
    ----------------
    kfloor : int > 0
        Quality of the approximation of the floor function, larger is better.
    kswitch : int > 0
        The quality of the switch/max/min functions. Larger is better.
    FBi : float, default = 15 * 15
        Average fish biomass at stocking.
    FP_0 : int, default = 10
    FB_0 : float, default = FBi * FP0
    FPs_0 : int, default=0
    FA_0 : float, default=0
    FC_0 : float, default=0
    INC_0 : float, default=1.63

    Fish-specific Available kwargs
    ------------------------------
    (coming soon, see paper)

    Connecting Gekko Variables Registered
    -------------------------------------
    FPs : m.MV (fish / day / pond)
        Input, user-supplied.
        Fish stocking number.
    FA : m.MV (kcal / day / pond)
        Input. Receives from fish anabolism.
        Fish anabolism.
    FC : m.MV (kcal / day / pond)
        Input. Receives from fish catabolism.
        Fish catabolism.
    INC : m.MV (mg N / l)
        Input. Receives from nitrogen connection.
        Total inorganic nitrogen concentration.
    FB : m.SV (kcal / pond)
        Output.
        Total fish biomass.
    FP : m.SV (fish / pond)
        Output.
        Total fish population.

    Returns
    -------
    m : GEKKO
        The gekko model with the fish component added.
    """

    # ------------------
    # Initial Conditions
    # ------------------

    FBi_val = kwargs.get('FBi', 15 * 15)
    FP_0 = kwargs.get('FP_0', 10)
    FB_0 = kwargs.get('FB_0', FBi_val * FP_0)
    FPs_0 = kwargs.get('FPs_0', 0)
    FA_0 = kwargs.get('FA_0', 0)
    FC_0 = kwargs.get('FC_0', 0)
    INC_0 = kwargs.get('INC_0', 1.63)

    # ----------
    # Parameters
    # ----------

    kfloor = kwargs.get('kfloor', 10)
    kswitch = kwargs.get('kswitch', 100)

    # Tilapia params from paper
    FBi = m.Param(value=FBi_val)
    kNH = m.Param(value=kwargs.get('kNH', 0.365))

    # --------------------
    # Connecting Variables
    # --------------------

    FB = _register_connecting(m, 'FB', m.SV, FB_0, lb=0)
    FP = _register_connecting(m, 'FP', m.SV, FP_0, lb=0)

    FPs = _register_connecting(m, 'FPs', m.MV, FPs_0)
    FA = _register_connecting(m, 'FA', m.MV, FA_0)
    FC = _register_connecting(m, 'FC', m.MV, FC_0)
    INC = _register_connecting(m, 'INC', m.MV, INC_0)

    FPdt = _register_connecting(m, 'FPdt', m.SV, value=0)

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

    return m


def register_anabolism(m, **kwargs):
    """The fish anabolism component of the fish model from [2].

    Parameters
    ----------
    m : GEKKO
        The gekko model.

    Available kwargs
    ----------------
    kswitch : int > 0
        The quality of the switch/max/min functions. Larger is better.
    FA_0 : float, default = 0
    AFC_0 : float, default = 0
    HFC_0 : float, default = 0


    Fish-specific Available kwargs
    ------------------------------
    (coming soon, see paper)

    Connecting Gekko Variables Registered
    -------------------------------------
    FB : m.MV (kcal / pond)
        Input. Receives from fish population.
        Fish biomass.
    FA : m.SV (kcal / day / pond)
        Output.
        Fish anabolism.
    AFC : m.SV (kcal / day / pond)
        Output.
        Autotrophic food consumption.
    HFC : m.SV (kcal / day / pond)
        Output.
        Heterotrophic food consumption.

    Returns
    -------
    m : GEKKO
        The gekko model with the fish component added.
    """

    # ------------------
    # Initial Conditions
    # ------------------

    FB_0 = kwargs.get()
    FA_0 = kwargs.get('FA_0', 0)
    AFC_0 = kwargs.get('AFC_0', 0)
    HFC_0 = kwargs.get('HFC_0', 0)

    # ----------
    # Parameters
    # ----------

    kswitch = kwargs.get('kswitch', 100)

    # Tilapia params from paper

    # --------------------
    # Connecting Variables
    # --------------------

    FB = _register_connecting(m, 'FB', m.MV, FB_0)

    FA = _register_connecting(m, 'FA', m.SV, FA_0)
    AFC = _register_connecting(m, 'AFC', m.SV, AFC_0)
    HFC = _register_connecting(m, 'HFC', m.SV, HFC_0)

    # ---------------------------
    # Equations and Intermediates
    # ---------------------------

    # Autotrophic food consumption
    AFC = FB
    # Total food consumption (kcal / day / pond)
    TFC = AFC + HFC + SFC
    # Anabolism
    m.Equation(FA == kmaxa * TFC * fFQ)

    return m


def old_register_fish(m, **kwargs):
    """OBSOLETE!!

    Registers the fish component.

    Model taken from [2].

    Not that in order to use this sub-module, the solver should be set to 1 in
    order to get integer population counts. However, solver=3 can be used for a
    relaxed (floating point) population count.

    Parameters
    ----------
    m : GEKKO
        The gekko model.

    Available kwargs
    ----------------
    T_0 : float or array, default = 20 (deg C)
        The temperature over time, if a float, the temperature is constant
        across time.
    FB_0 : float, default = (kcal/tank)
        The initial fish biomass.
    FP_0 : int, default = (fish/tank)
        The initial fish population.
    FBi : float, default = (kcal/fish)
        Average individual fish biomass at stocking.
    floor_n : int > 0
        The precision for computing the floor function (higher is better).

    Fish-specific Available kwargs
    ------------------------------
    km2 : float, default = (1/day)
        Fish mortality coefficient.
    kmaxa : float, default = (unitless)
        Maximum assimalation coaefficient.
    Tmaxf : float, default = (deg C)
        Fish maximum temperature.
    Toptf : float, default = (deg C)
        Fish optimum temperature.
    Q10 : float, default = (unitless)
        Relative increase in the rate of a biological activity with an increase
        in temperature of 10 degC.
    kPE : float, default = ((g protein / kcal)^-2)
        Coefficient of PE on food assimilation

    Connecting Gekko Variables Registered
    -------------------------------------
    T : MV (deg C)
        Receives from environment.
        The temperature over time.
    FB : m.SV (kcal/tank)
        Measurable output.
        Total fish biomass.
    FP : m.SV (fish/tank)
        Measurable output.
        Total fish population.
    FPs : m.MV (fish/tank/day)
        Input.
        Fish stocking number.

    Returns
    -------
    m : GEKKO
        The gekko model with the fish component added.
    """
    # ------------------
    # Initial Conditions
    # ------------------

    T_0 = kwargs.get('T0', 0)
    FB_0 = kwargs.get('FB_0', 0)  # TODO
    FP_0 = int(kwargs.get('FP_0', 0))  # TODO

    # ----------
    # Parameters
    # ----------

    # Tilapia params from paper
    kp1 = m.Param(value=kwargs.get('kp1', 0.14))
    mumax = m.Param(value=kwargs.get('mumax', 1.6))
    km1 = m.Param(value=kwargs.get('km1', 0.6))
    kr = m.Param(value=kwargs.get('kr', 0.1))
    hN = m.Param(value=kwargs.get('hN', 0.2))
    hP = m.Param(value=kwargs.get('hP', 0.02))
    ks = m.Param(value=kwargs.get('ks', 0.14))
    kd = m.Param(value=kwargs.get('kd', 0.12))
    ha = m.Param(value=kwargs.get('ha', 0.51))
    hh = m.Param(value=kwargs.get('hh', 0.05))
    a = m.Param(value=kwargs.get('a', 0.000017))
    b = m.Param(value=kwargs.get('b', 0.000015))
    kT1 = m.Param(value=kwargs.get('kT1', 0.004))
    kT2 = m.Param(value=kwargs.get('kT2', 0.008))
    kp2 = m.Param(value=kwargs.get('kp2', 0.12))
    kmaxa = m.Param(value=kwargs.get('kmaxa', 0.75))
    kfeed = m.Param(value=kwargs.get('kfeed', 0.31))
    kfast = m.Param(value=kwargs.get('kfast', 0.005))
    n = m.Param(value=kwargs.get('n', -0.12))
    m = m.Param(value=kwargs.get('m', -0.3))
    FBs = m.Param(value=kwargs.get('FBs', 20))
    Q10 = m.Param(value=kwargs.get('Q10', 2.37))
    c = m.Param(value=kwargs.get('c', 0.59))
    d = m.Param(value=kwargs.get('d', 0.027))
    DOcrit = m.Param(value=kwargs.get('DOcrit', 1.0))
    PEmin = m.Param(value=kwargs.get('PEmin', 0.025))
    PEopt = m.Param(value=kwargs.get('PEopt', 0.09))
    kPEval = kwargs.get('kPE', 0.45)
    kPE = m.Param(value=kPEval)
    kn = m.Param(value=kwargs.get('kn', 0.01))
    knf = m.Param(value=kwargs.get('knf', 0.47))
    FAPPmax = m.Param(value=kwargs.get('FAPPmax', 0.17))
    kn1 = m.Param(value=kwargs.get('kn1', 0.17))
    ksn = m.Param(value=kwargs.get('ksn', 0.003))
    kap = m.Param(value=kwargs.get('kap', 0.001))
    kpr = m.Param(value=kwargs.get('kpr', 0.06))
    kps = m.Param(value=kwargs.get('kps', 28))
    DICO = m.Param(value=kwargs.get('DICO', 2))
    AWFT = m.Param(value=kwargs.get('AWFT', 0.018))
    kao = m.Param(value=kwargs.get('kao', 0.28))
    kho = m.Param(value=kwargs.get('kho', 0.28))
    kfo = m.Param(value=kwargs.get('kfo', 0.28))
    ksno = m.Param(value=kwargs.get('ksno', 0.28))
    khp = m.Param(value=kwargs.get('khp', 0.001))
    kDOT = m.Param(value=kwargs.get('kDOT', 4.0))
    kNHT = m.Param(value=kwargs.get('kNHT', 4.0))
    kDOf = m.Param(value=kwargs.get('kDOf', 2.5))
    kfn = m.Param(value=kwargs.get('kfn', 0.017))
    kfp = m.Param(value=kwargs.get('kfp', 0.002))
    kNH = m.Param(value=kwargs.get('kNH', 0.365))
    s = m.Param(value=kwargs.get('s', 21.08))
    kDO = m.Param(value=kwargs.get('kDO', 0.14))
    kTI = m.Param(value=kwargs.get('kTI', 0.012))

    Tmaxf = m.Param(value=kwargs.get('Tmaxf', 41))
    Toptf = m.Param(value=kwargs.get('Toptf', 30))

    # FBi = m.Param(value=kwargs.get('FBi', 0))  # TODO
    # km2 = m.Param(value=kwargs.get('km2', 0.2))  # TODO
    # kmaxa = m.Param(value=kwargs.get('kmaxa', 0))  # TODO
    # Tmaxf = m.Param(value=kwargs.get('Tmaxf', 0))  # TODO
    # Toptf = m.Param(value=kwargs.get('Toptf', 0))  # TODO
    # FC = m.Param(value=1)  # TODO move to variable

    # --------------------
    # Connecting Variables
    # --------------------

    T = _register_connecting(m, 'T', m.MV, T_0)
    FB = _register_connecting(m, 'FB', m.SV, FB_0, lb=0)
    FP = _register_connecting(m, 'FP', m.SV, FP_0, lb=0)
    FPs = _register_connecting(m, 'FPs', m.MV, 0)

    # ------------------
    # Internal Variables
    # ------------------

    # Total dissolved inorganic nitrogen (gN / pond)
    TIN = m.SV(value=0, lb=0)
    # Total nitrogen in sediment (gN / pond)
    TNS = m.SV(value=0, lb=0)
    # Total dissolved phosphorus (gP / pond)
    TIP = m.SV(value=0, lb=0)
    # Total phosphorus in sediment (gP / pond)
    TPS = m.SV(value=0, lb=0)
    # Dissolved oxygen quantity in pond (g DO / pond)
    QDO = m.SV(value=0, lb=0)

    # ---------------------------
    # Equations and Intermediates
    # ---------------------------

    # Mean fish biomass (kcal/fish)
    FBm = m.Intermediate(FB / FP)

    # ..............
    # Fish Anabolism
    # ..............

    # Effect of temperature on fish food consumption
    V = m.Intermediate((Tmaxf - T) / (Tmaxf - Toptf))
    S1 = m.Intermediate(m.log(Q10 * (Tmaxf - Toptf)))
    S2 = m.Intermediate(m.log(Q10 * (Tmaxf - Toptf + 2)))
    x = m.Intermediate((S1 ** 2 * (1 + (1 + 40 / S2) ** 0.5) ** 2) / 400)
    f2T = m.Intermediate(V * m.exp(x * (1 - V)))
    # Toxicity index
    TI = m.Intermediate(kDOT * (1 - f1DO) + kNHT * INC * kNH)
    fWQ = m.Intermediate(m.exp(-kTI * TI ** 2))
    # Fish appetite satiation
    FAPP = m.Intermediate(FAPPmax * FB * fsm * fWQ)
    # Autotrophic food availability
    fa = m.Intermediate(1 - m.exp(-s * (AFe / FB) ** 2.2))
    # Heteotrophic food availability
    fh = m.Intermediate(1 - m.exp(-s * (HFe / FB) ** 2.2))
    # The effect of fish size on food consumption
    fsm = m.Intermediate((FBm / FBs) ** m)
    # Autotrophic Fish Consumption (kcal/day/pond)
    AFC = m.Intermediate(FB * ha * fa * fsm * f2T * fWQ)
    # Heteotrophic Fish Consumption (kcal/day/pond)
    HFC = m.Intermediate(FB * hh * fh * fsm * f2T * fWQ)
    # Total Fish Consumption (kcal/day/pond)
    TFC = m.Intermediate(AFC + HFC + SFC)
    # Requirement for supplementary feed
    RSF = m.switch(FAPP - (AFC + HFC), 0, FAPP, AFC + HFC)
    # Actual consumption of supplementary feed
    SFC = m.min(RSF, SFA)
    # Fecal Waste
    FW = m.Intermediate(SFA - SFC + TFC - FA)
    # P:E ratio of total consumed food by fish and requirement of protein
    #   content from supplementary feed (g protein / kcal)
    PE = m.Intermediate((AFC * kp1 + HFC * kp2 + SFC * kp3) / TFC)
    # Required protein content of supplementary feed (g protein / kcal)
    Rkp3 = m.Intermediate((PEopt * TFC - AFC * kp1 - HFC * kp2) / SFC)
    # Effect of P:E ratio on food assimilation
    fFQ_up = m.switch(
        m.exp(-kPE * ((PEopt - PE) / (PEopt - PEmin)) ** 0.85),
        1, PE, PEopt
    )
    fFQ = m.switch(m.exp(-kPE), fFQ_up, PEmin, PE)
    # Fish Anabolism (kcal/day/pond)
    FA = m.Intermediate(kmaxa * TFC * fFQ)

    # ...............
    # Fish Catabolism
    # ...............

    # Fish respiration os a function of DO and T
    f2DO = m.switch(
        m.exp(-kDOf * (DOcrit - DOC) ** 2),
        1, DO, DOcrit
    )
    f3T = m.Intermediate(c + d * T)
    # Fish catabolism
    fsn = m.Intermediate((FBm / FBs) ** n)
    FC = m.Intermediate(TFC * kfeed + FB * kfast * fsn * f2DO * f3T)

    # ..........................
    # Elementary nutrient sector
    # ..........................

    # N-fixation rate
    FIXN = kn * AFe * m.exp(-knf * INC ** 2)
    # Nitrogen and phosphorus dynamics
    m.Equation(TIN.dt() == FC * kfn + AFR * kan + HFD * khn + TNS * ksn +
               FTN - (AFG * kan - FIXN) - TIN * kn1)
    m.Equation(TNS.dt() == HFS * khn - TNS * ksn)
    m.Equation(TIP.dt() == FC * kfp + AFR * kap + HFD * khp + TPS * kpr / dw +
               FTP - AFG * kap - TIP * kps / dw)
    m.Equation(TPS.dt() == TIP * kps / dw + HFS * khp - TPS * kpr / dw)

    # .......................
    # Dissolved oxygen sector
    # .......................

    # Oxygen exchange rate between air and water
    DOdf = (DICO * (DOC - DOS) / (AWFT * dw) * Vpond) / 1000
    # Dissolved oxygen dynamics
    m.Equation(QDO.dt() == AFG * kao + DOdf - FC * kfo - AFR * kao -
               HFD * kho - TNS * ksn / khn * ksno)

    # ....................
    # Food nutrient sector
    # ....................

    # Function representing elementary nutrient limitation
    fNP = m.min(INC / (INC + hN), IPC / (IPC + hP))
    # Function representing elementary light limitation
    fI = m.Intermediate(I0 + m.exp(-(a * AFe + b * HFe)) / Ir)
    # Phytoplankton growth vs. temperature relationship
    f1T = m.switch(
        m.exp(-kT1 * (T - Topta) ** 2),
        m.exp(-kT2 * (Topta - T) ** 2),
        T, Topta
    )
    # Decomposition of heterotrophic particles
    f1DO = m.switch(m.exp(-kDO * (4 - DOC) ** 2), 1, DO, 4)
    # Autotrophic food nutrient energy and protein dynamics
    AFG = m.Intermediate(mumax * AFe * fNP * fI * f1T)
    AFR = m.Intermediate(AFe * kr)
    AFM = m.Intermediate(AFe * km1)
    m.Equation(AFe.dt() == AFG - AFC - AFR - AFM)
    m.Equation(AFp.dt() == AFe.dt() * kp1)
    # Heterotrophic food nutrient energy and protein dynamics
    HFS = m.Intermediate(HFe * ks)
    HFD = m.Intermediate(HFe * kd * f1DO)
    m.Equation(HFe.dt() == AFM + FW - HFC - HFS - HFD)
    m.Equation(HFp.dt() == HFe.dt() * kp2)

    # ...................
    # Fish biomass growth
    # ...................

    m.Equation(FB.dt() == FPs * FBi + FA - FC - FP * km2 * FBm)
    m.Equation(FP.dt() == FPs - m.floor(FP * km2, kwargs.get('floor_n', 10)))

    return m


def register_bacteria(m, **kwargs):
    """Registers the bacteria component.

    Model taken from [1].

    Parameters
    ----------
    m : GEKKO
        The gekko model.

    Available kwargs
    ----------------
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
    Em : float, default = 0.05 (unitless)
        Production rate of dry mass of nitrosomonas from ammonia.
    Eb : float, default = 0.02 (unitless)
        Production rate of dry mass of nitrobacter from nitrite.
    fm : float in [0,1], default = 0.99
        The ratio of the mass of nitrite formed to that of ammonia oxidized.
    fn : float in [0, 1], default = 0.99
        The ratio of the mass of nitrate formed to that of nitrite oxidized.

    Connecting Gekko Variables Registered
    -------------------------------------
    T : MV (deg C)
        Receives from environment.
        The temperature over time.
    NH3_exc : MV (mg/l)
        Receives from fish.
        The *total* (integral over time) concentration of ammonia produced by
        the fish.
    NH3 : MV (mg/l)
        Sends to fish.
        The total concentration of ammonia.
    NO2 : MV (mg/l)
        Sends to fish
        The total concentration of nitrite.
    NO3 : MV (mg/l)
        Sends to lettuce.
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

    Returns
    -------
    m : GEKKO
        The gekko model with the bacteria component added.
    """

    # ------------------
    # Initial Conditions
    # ------------------

    T_0 = kwargs.get('T0', 0)
    Cm_0 = kwargs.get('Cm_0', 0.0025)
    Cb_0 = kwargs.get('Cb_0', 0.0005)
    NH3_0 = kwargs.get('NH3_0', 0)
    NO2_0 = kwargs.get('NO2_0', 0)
    NO3_0 = kwargs.get('NO3_0', 0)

    # ----------
    # Parameters
    # ----------

    Em = m.Param(value=kwargs.get('Em', 0.05))
    Eb = m.Param(value=kwargs.get('Eb', 0.02))
    fm = m.Param(value=kwargs.get('fm', 0.99))
    fn = m.Param(value=kwargs.get('fn', 0.99))

    # --------------------
    # Connecting Variables
    # --------------------

    T = _register_connecting(m, 'T', m.MV, T_0)
    NH3_exc = _register_connecting(m, 'NH3_exc', m.MV, 0)
    NH3 = _register_connecting(m, 'NH3', m.SV, NH3_0, lb=0)
    NO2 = _register_connecting(m, 'NO2', m.SV, NO2_0, lb=0)
    NO3 = _register_connecting(m, 'NO3', m.SV, NO3_0, lb=0)
    NO2_up = _register_connecting(m, 'NO2_up', m.MV, 0)
    NO3_up = _register_connecting(m, 'NO3_up', m.MV, 0)
    Cm = _register_connecting(m, 'Cm', m.SV, Cm_0, lb=0)
    Cb = _register_connecting(m, 'Cb', m.SV, Cb_0, lb=0)

    # ---------------------------
    # Equations and Intermediates
    # ---------------------------

    # Growth Constants
    km = m.Intermediate(10 ** (0.0413 * T - 0.944))
    kb = m.Intermediate(10 ** (0.0255 * T - 0.492))
    # Saturation constants
    NH3_m = m.Intermediate(10 ** (0.051 * T - 1.158))
    NO2_m = m.Intermediate(10 ** (0.063 * T - 1.149))
    # Growth of Nitrosomonas
    m.Equation(Cm.dt() == km * Cm * NH3 / (NH3 + NH3_m))
    # Growth of Nitrobacter
    m.Equation(Cb.dt() == kb * Cb * NO2 / (NO2 + NO2_m))
    # Oxidation of ammonia
    m.Equation(Cm - Cm_0 == Em * (NH3_0 + NH3_exc - NH3))
    # m.Equation(NH3.dt() == -Cm.dt() / Em + NH3_exc.dt())
    # # Oxidation of nitrite
    m.Equation(
        Cb - Cb_0 == Eb * (NO2_0 - NO2_up + fm * (NH3_0 + NH3_exc - NH3) - NO2)
    )
    # m.Equation(NO2.dt() == -Cb.dt() / Eb - NO2_up.dt() - fm * NH3.dt())
    # # Balance of nitrate
    # m.Equation(NO3 == NO3_0 + fm * (NH3_0 + NH3_exc - NH3) +
    #            fn * (NO2_0 - NO2_up - NO2) - NO3_up)
    m.Equation(NO3 == NO3_0 - NO3_up + fn * (
        NO2_0 - NO2_up + fm * (NH3_0 + NH3_exc - NH3) - NO2)
    )

    return m

###############################################################################
#   Private Helpers
###############################################################################


def _register_connecting(m, name, type, value, **kwargs):
    """Registers a connecting variable with the model, if not already
    registered. If registered, but type is `m.SV`, this will overwrite the old
    registration.

    This variable will be available as an MV at `m.aqua[name]`, and will only
    be created if `m.aqua[name]` does not already exist.

    Parameters
    ----------
    m : GEKKO
        The gekko model to which to register the connecting variable.
    name : str
        The name of the connecting variable.
    type : m.Var, m.MV, etc
        The type of variable to register.
    value : obj
        The value to assign to the connecting variable.
    kwargs : key word arguments
        Any additional arguments to pass in to the variable constructor.

    Returns
    -------
    var : MV
        A handle to the connecting variable, whether it was created or already
        exists.
    """
    if not hasattr(m, 'aqua'):
        m.aqua = {}

    if name not in m.aqua or type == m.SV:
        m.aqua[name] = type(value=value, **kwargs)

    return m.aqua[name]
