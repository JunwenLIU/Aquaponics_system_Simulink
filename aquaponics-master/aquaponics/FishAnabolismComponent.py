from .CoreComponent import CoreComponent


class FishAnabolismComponent(CoreComponent):
    """Registers the heterotrophic nutrient sector component of the fish model.

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
    FB : MV (kcal / pond)
        Input. Receives from fish growth.
        Total fish biomass.
    FP : MV (fish / pond)
        Input. Receives from fish growth.
        Total fish population.
    AFe : MV (kcal / pond)
        Input. Receives from autotrophic.
        The autotrophic food quantity in terms of energy.
    HFe : MV (kcal / pond)
        Input. Receives from heterotrophic.
        The heterotrophic food quantity in terms of energy.
    f1DO : MV
        Input. Receives from heterotrophic.
        Decomposition of heterotrophic particles.
    INC : MV
        Input. Receives from nitrogen connection.
        Incorporated nitrogen.
    FA : m.SV (kcal / day / pond)
        Output.
        Fish anabolism.
    AFC : m.SV (kcal / day / pond)
        Output.
        Autotrophic food consumption.
    HFC : m.SV (kcal / day / pond)
        Output.
        Heterotrophic food consumption.
    FW : m.SV (kcal / day / pond)
        Output.
        Fish fecal waste.
    TFC : m.SV (kcal / day / pond)
        Output.
        Total food consumption.
    """

    def __init__(self, m, aqua):
        CoreComponent.__init__(self, m, aqua)

    def register_connecting(self, **kwargs):
        """Registers all connecting variables associated with this component to
        the gekko model.

        kwargs
        ------
        FBi : float, default = 15 * 15
        FP_0 : float, default = 10
        FB_0 : float, default = FBi * FP_0
        AFe_0 : float, default = 0
        HFe_0 : float, default = 0
        FA_0 : float, default = 0
        AFC_0 : float, default = 0
        HFC_0 : float, default = 0
        f1DO_0 : float, default = 0
        INC_0 : float, default = 1.63 (mg N / l)
        FW_0 : float, default = 0
        """
        # Aliases
        a = self.aqua
        MV = self.m.MV
        SV = self.m.SV

        # Initial Conditions
        T_0 = kwargs.get('T_0', 20)
        FBi_val = kwargs.get('FBi', 15 * 15)
        FP_0 = kwargs.get('FP_0', 10)
        FB_0 = kwargs.get('FB_0', FBi_val * FP_0)
        FA_0 = kwargs.get('FA_0', 0)
        AFC_0 = kwargs.get('AFC_0', 0)
        HFC_0 = kwargs.get('HFC_0', 0)
        AFe_0 = kwargs.get('AFe_0', 0)
        HFe_0 = kwargs.get('HFe_0', 0)
        f1DO_0 = kwargs.get('f1DO_0', 0)
        INC_0 = kwargs.get('INC_0', 1.63)
        FW_0 = kwargs.get('FW_0', 0)

        # Connecting Variables
        a.register_connecting('T', MV, T_0)
        a.register_connecting('FB', MV, FB_0)
        a.register_connecting('FP', MV, FP_0)
        a.register_connecting('AFe', MV, AFe_0)
        a.register_connecting('HFe', MV, HFe_0)
        a.register_connecting('f1DO', MV, f1DO_0)
        a.register_connecting('INC', MV, INC_0)

        a.register_connecting('FA', SV, FA_0)
        a.register_connecting('AFC', SV, AFC_0)
        a.register_connecting('HFC', SV, HFC_0)
        a.register_connecting('FW', SV, FW_0)
        a.register_connecting('TFC', SV, 0)

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
        s = m.Param(value=kwargs.get('s', 21.08))
        _m = m.Param(value=kwargs.get('m', -0.3))
        Q10 = m.Param(value=kwargs.get('Q10', 2.37))
        ha = m.Param(value=kwargs.get('ha', 0.51))
        hh = m.Param(value=kwargs.get('hh', 0.05))
        FBs = m.Param(value=kwargs.get('FBs', 20))
        kmaxa = m.Param(value=kwargs.get('kmaxa', 0.75))
        kDOT = m.Param(value=kwargs.get('kDOT', 4.0))
        kNHT = m.Param(value=kwargs.get('kNHT', 4.0))
        kNH = m.Param(value=kwargs.get('kNH', 0.365))
        kTI = m.Param(value=kwargs.get('kTI', 0.012))
        FAPPmax = m.Param(value=kwargs.get('FAPPmax', 0.17))
        PEmin = m.Param(value=kwargs.get('PEmin', 0.025))
        PEopt = m.Param(value=kwargs.get('PEopt', 0.09))
        kPEval = kwargs.get('kPE', 0.45)
        kPE = m.Param(value=kPEval)
        kp1 = m.Param(value=kwargs.get('kp1', 0.14))
        kp2 = m.Param(value=kwargs.get('kp2', 0.12))

        Tmaxf = m.Param(value=kwargs.get('Tmaxf', 41))
        Toptf = m.Param(value=kwargs.get('Toptf', 30))

        # --------------------
        # Connecting Variables
        # --------------------

        T = a.T
        FB = a.FB
        FP = a.FP
        FA = a.FA
        AFC = a.AFC
        HFC = a.HFC
        AFe = a.AFe
        HFe = a.HFe
        f1DO = a.f1DO
        INC = a.INC
        FW = a.FW
        TFC = a.TFC

        # ---------------------------
        # Equations and Intermediates
        # ---------------------------

        # Mean fish biomass
        FBm = m.Intermediate(FB / FP)
        # Natural food availability
        fa = m.Intermediate(1 - m.exp(-s * (AFe / FB) ** 2.2))
        fh = m.Intermediate(1 - m.exp(-s * (HFe / FB) ** 2.2))
        # Effect of fish size on food consumption
        fsm = m.Intermediate((FBm / FBs) ** _m)
        # Effect of temperature on fish food consumption
        S1 = m.Intermediate(m.log(Q10 * (Tmaxf - Toptf)))
        S2 = m.Intermediate(m.log(Q10 * (Tmaxf - Toptf + 2)))
        x = m.Intermediate(((S1 ** 2) * (1 + (1 + 40 / S2) ** 0.5) ** 2))
        V = m.Intermediate((Tmaxf - T) / (Tmaxf - Toptf))
        f2T = m.Intermediate(V * m.exp(x * (1 - V)))
        # Toxicity index
        TI = m.Intermediate(kDOT * (1 - f1DO) + kNHT * INC * kNH)
        fWQ = m.Intermediate(m.exp(-kTI * TI ** 2))
        # Fish appetite satiation
        FAPP = FAPPmax * FB * fsm * fWQ
        # Grazing rate on autotrophic and heterotrophic food sources
        m.Equation(AFC == FB * ha * fa * fsm * f2T * fWQ)
        m.Equation(HFC == FB * hh * fh * fsm * f2T * fWQ)
        # Natural food consumption
        RSF = m.switch(0, FAPP - (AFC + HFC), FAPP, AFC + HFC, k=kswitch)
        # Supplementary feed availability (set to RSF so that supply exactly
        #   meets demand)
        SFA = RSF
        # Supplementary feed consumption
        SFC = m.min(RSF, SFA)
        # Total food consumption (kcal / day / pond)
        m.Equation(TFC == AFC + HFC + SFC)
        # Fish fecal waste
        m.Equation(FW == SFA - SFC + TFC - FA)
        # P:E ratio of total consumed food by fish and the requirement of
        #   protein content from supplementary feed
        kp3 = m.Intermediate((PEopt * TFC - AFC * kp1 - HFC * kp2) / SFC)
        PE = m.Intermediate((AFC * kp1 + HFC * kp2 + SFC * kp3) / TFC)
        # Effect of P:E ratio on food assimilation
        fFQ_up = m.switch(
            m.exp(-kPE * ((PEopt - PE) / (PEopt - PEmin)) ** 0.85),
            1, PE, PEopt, k=kswitch
        )
        fFQ = m.switch(m.exp(-kPE), fFQ_up, PEmin, PE, k=kswitch)
        # Anabolism
        m.Equation(FA == kmaxa * TFC * fFQ)
