from .CoreComponent import CoreComponent


class FishNutrientComponent(CoreComponent):
    """Registers the elementary nutrient sector component of the fish model.

    Model taken from [2] (see Aquaponics.py for sources).

    Parameters
    ----------
    m : GEKKO
        The gekko model.

    Connecting Variables
    --------------------
    dw : FV (cm)
        User input.
        The fish tank depth.
    FTN : MV (gN / day / pond)
        User input. Should probably stay at zero.
        Inorganic nitrogen added from fertilization.
    FTP : MV (gP / day / pond)
        User input. Should probably stay at zero.
        Inorganic phosphorus added from fertilization.
    AFe : MV (kcal / pond)
        Input. Receives from autotrophic.
        The autotrophic food quantity in terms of energy.
    INC : MV (mg N / l)
        Input. Receives from nitrogen connection.
        Total inorganic nitrogen concentration.
    FC : float (kcal / day / pond)
        Input. Receives from catabolism.
        Fish catabolism.
    AFR : m.SV (kcal / day / pond)
        Input. Receives from autotrophic.
        Autotrophic food loass due to phytoplankton respiration.
    AFG : m.SV (kcal / day / pond)
        Input. Receives from autotrophic.
        Autotrophic food loass due to phytoplankton growth.
    HFS : SV (kcal / day / pond)
        Input. Receives from heterotrophic.
        Heterotrophic food loss rate due to sedimentation.
    HFD : SV (kcal / day / pond)
        Input. Receives from heterotrophic.
        Heterotrophic food loss rate due to decomposition.
    TIN : float (gN / pond)
        Output.
        Total dissolved inorganic nitrogen.
    TIP : float (gN / pond)
        Output.
        Total nitrogen.
    TNS : float (gP / pond)
        Output.
        Total dissolved inorganic phosphorus.
    TPS : float (gP / pond)
        Output.
        Total phosphorus.
    """

    def __init__(self, m, aqua):
        CoreComponent.__init__(self, m, aqua)

    def register_connecting(self, **kwargs):
        """Registers all connecting variables associated with this component to
        the gekko model.

        kwargs
        ------
        dw : float, default = 0
        TIN_0 : float, default = 0
        TNS_0 : float, default = 0
        TIP_0 : float, default = 0
        TPS_0 : float, default = 0
        AFe_0 : float, default = 0
        INC_0 : float, default = 0
        FC_0 : float, default = 0
        AFR_0 : float, default = 0
        AFG_0 : float, default = 0
        HFS_0 : float, default = 0
        HFD_0 : float, default = 0
        """
        # Aliases
        a = self.aqua
        MV = self.m.MV
        SV = self.m.SV
        FV = self.m.FV

        # Initial Conditions
        dw = kwargs.get('dw', 40)
        TIN_0 = kwargs.get('TIN_0', 0)
        TNS_0 = kwargs.get('TNS_0', 0)
        TIP_0 = kwargs.get('TIP_0', 0)
        TPS_0 = kwargs.get('TPS_0', 0)

        AFe_0 = kwargs.get('AFe_0', 0)
        INC_0 = kwargs.get('INC_0', 0)
        FTN_0 = kwargs.get('FTN_0', 0)
        FTP_0 = kwargs.get('FTP_0', 0)
        FC_0 = kwargs.get('FC_0', 0)
        AFR_0 = kwargs.get('AFR_0', 0)
        AFG_0 = kwargs.get('AFG_0', 0)
        HFS_0 = kwargs.get('HFS_0', 0)
        HFD_0 = kwargs.get('HFD_0', 0)

        # Connecting Variables
        a.register_connecting('dw', FV, dw)
        a.register_connecting('TIN', SV, TIN_0)
        a.register_connecting('TNS', SV, TNS_0)
        a.register_connecting('TIP', SV, TIP_0)
        a.register_connecting('TPS', SV, TPS_0)

        a.register_connecting('FTN', MV, FTN_0)
        a.register_connecting('FTP', MV, FTP_0)
        a.register_connecting('AFe', MV, AFe_0)
        a.register_connecting('INC', MV, INC_0)
        a.register_connecting('FC', MV, FC_0)
        a.register_connecting('AFR', MV, AFR_0)
        a.register_connecting('AFG', MV, AFG_0)
        a.register_connecting('HFS', MV, HFS_0)
        a.register_connecting('HFD', MV, HFD_0)

    def register_equations(self, **kwargs):
        """Registers all equations and intermediates associated with this
        component to the gekko model.

        kwargs
        ------
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
        kn = m.Param(value=kwargs.get('kn', 0.01))
        knf = m.Param(value=kwargs.get('knf', 0.47))
        kfn = m.Param(value=kwargs.get('kfn', 0.017))
        kfp = m.Param(value=kwargs.get('kfp', 0.002))
        kp1 = m.Param(value=kwargs.get('kp1', 0.14))
        kp2 = m.Param(value=kwargs.get('kp2', 0.12))
        ksn = m.Param(value=kwargs.get('ksn', 0.003))
        kn1 = m.Param(value=kwargs.get('kn1', 0.17))
        khp = m.Param(value=kwargs.get('khp', 0.001))
        kpr = m.Param(value=kwargs.get('kpr', 0.06))
        kps = m.Param(value=kwargs.get('kps', 28))
        kap = m.Param(value=kwargs.get('kap', 0.001))

        # --------------------
        # Connecting Variables
        # --------------------

        dw = a.dw

        TIN = a.TIN
        TNS = a.TNS
        TIP = a.TIP
        TPS = a.TPS

        AFe = a.AFe
        INC = a.INC
        FTN = a.FTN
        FTP = a.FTP
        FC = a.FC
        AFR = a.AFR
        AFG = a.AFG
        HFS = a.HFS
        HFD = a.HFD

        # ---------------------------
        # Equations and Intermediates
        # ---------------------------

        kan = m.Intermediate(kp1 / 6.25)
        khn = m.Intermediate(kp2 / 6.25)
        # N-fixation rate
        FIXN = m.Intermediate(kn * AFe * m.exp(-knf * INC) ** 2)
        # Nitrogen and phosphorus dynamics
        m.Equation(TIN.dt() == FC * kfn + AFR * kan + HFD * khn + TNS * ksn +
                   FTN - (AFG * kan - FIXN) - TIN * kn1)
        m.Equation(TNS.dt() == HFS * khn - TNS * ksn)
        m.Equation(TIP.dt() == FC * kfp + AFR * kap + HFD * khp +
                   TPS * kpr / dw + FTP - AFG * kap - TIP * kps / dw)
        m.Equation(TPS.dt() == TIP * kps / dw + HFS * khp - TPS * kpr / dw)
