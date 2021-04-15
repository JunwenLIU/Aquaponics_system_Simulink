

class CoreComponent:
    """Base class for aquaponics components (e.g. bacteria, etc.).

    Parameters
    ----------
    m : GEKKO
        The gekko model with which this component is built.
    aqua : Aquaponics
        The aquaponics system.
    """

    def __init__(self, m, aqua):
        self.m = m
        self.aqua = aqua

    def register_connecting(self, **kwargs):
        """Registers all connecting variables associated with this component to
        the gekko model.

        Should be implemented by the child class.

        kwargs
        ------
        Initial conditions.
        """
        raise NotImplementedError(
            'register_connecting must be implemented by the child class'
        )

    def register_equations(self, **kwargs):
        """Registers all equations and intermediates associated with this
        component to the gekko model.

        Should be implemented by the child class.

        kwargs
        ------
        Model parameters.
        """
        raise NotImplementedError(
            'register_equations must be implemented by the child class'
        )
