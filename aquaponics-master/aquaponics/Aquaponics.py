import numpy as np
from gekko import GEKKO

from .gekko_extensions import register_extensions
from .BacteriaComponent import BacteriaComponent
from .FishAutotrophicComponent import FishAutotrophicComponent
from .FishHeterotrophicComponent import FishHeterotrophicComponent
from .FishGrowthComponent import FishGrowthComponent
from .FishAnabolismComponent import FishAnabolismComponent
from .FishCatabolismComponent import FishCatabolismComponent
from .FishNutrientComponent import FishNutrientComponent
from .NitrogenConnectionComponent import NitrogenConnectionComponent
from .PlantComponent import PlantComponent
from .HydroPlantComponent import HydroPlantComponent
from .HydroNitrogenComponent import HydroNitrogenComponent


class Aquaponics:
        """
        Management system for loading and running various aquaponics models.

        Note: the model created by this system will be wrapped with
        `gekko_extensions`.

        Usage
        -----
        To load only the bacteria component:

            a = Aquaponics('bacteria')
            m = a.get_model()

        To load heterotrophic and autotrophic food source components:

            a = Aquaponics('heterotrophic', 'autotrophic')
            m = a.get_model()

        To load all models:

            a = Aquaponics('all')
            m = a.get_model()

        To load all models with component arguments specified:

            a = Aquaponics('all', kswitch=50, kp1=2, NH3_0=2)
            m = a.get_model()

        args
        ----
        The names of individual components to load, or 'all' if all components
        should be loaded. If 'all' is passed as even one argument, even with
        other arguments, then all components will be loaded, regardless of what
        is given.

        Avilable arguments include:
            - 'all'
            - 'bacteria'
            - 'autotrophic'
            - 'heterotrophic'
            - 'fishgrowth'
            - 'anabolism'
            - 'nutrient'
            - 'nitrogen'
            - 'plant'

        kwargs
        ------
        Key word arguments to pass into the components as they are registered.

        Sources
        -------
        [1] Knowles, G., Downing, A. L., & Barrett, M. J. (1965). Determination
            of kinetic constants for nitrifying bacteria in mixed culture, with
            the aid of an electronic computer. Microbiology, 38(2), 263-278.
        [2] Li, L., & Yakupitiyage, A. (2003). A model for food nutrient
            dynamics of semi-intensive pond fish culture. Aquacultural
            Engineering, 27(1), 9-38.
        [3] Zhang, K., Burns, I. G., & Turner, M. K. (2008). Derivation of a
            dynamic model of the kinetics of nitrogen uptake throughout the
            growth of lettuce: calibration and validation. Journal of plant
            nutrition, 31(8), 1440-1460.
        [4] Evans, J. J., Pasnik, D. J., Brill, G. C., & Klesius, P. H. (2006).
            Un‐ionized Ammonia Exposure in Nile Tilapia: Toxicity, Stress
            Response, and Susceptibility to Streptococcus agalactiae. North
            American Journal of Aquaculture, 68(1), 23-33.
        """

        def __init__(self,  *args, **kwargs):
            if 'server' in kwargs:
                self._m = register_extensions(
                    GEKKO(server=kwargs.get('server'))
                )
                kwargs.pop('server', None)
            else:
                self._m = register_extensions(GEKKO(
                    name='{}'.format(np.random.randint(1000000000))
                ))

            available_components = {
                'bacteria': BacteriaComponent,
                'autotrophic': FishAutotrophicComponent,
                'heterotrophic': FishHeterotrophicComponent,
                'fishgrowth': FishGrowthComponent,
                'anabolism': FishAnabolismComponent,
                'catabolism': FishCatabolismComponent,
                'nutrient': FishNutrientComponent,
                'nitrogen': NitrogenConnectionComponent,
                'plant': PlantComponent,
                'hydroplant': HydroPlantComponent,
                'hydronitrogen': HydroNitrogenComponent
            }
            if 'fish' in args:
                args = list(args)
                args.remove('fish')
                for comp in ['autotrophic', 'heterotrophic', 'fishgrowth',
                             'anabolism', 'catabolism', 'nutrient',
                             'nitrogen']:
                    if comp not in args:
                        args.append(comp)
            if 'all' in args:
                to_register = [
                    'bacteria', 'autotrophic', 'heterotrophic', 'fishgrowth',
                    'anabolism', 'catabolism', 'nutrient', 'nitrogen', 'plant'
                ]
                registered_components = [
                    available_components[component]
                    for component in to_register
                ]
            else:
                registered_components = [
                    available_components[component] for component in args
                ]
            self.registered_components = registered_components
            # print(registered_components)
            # print(len(registered_components))

            # Initialize components
            registered_components = [
                component(self._m, self) for component in registered_components
            ]

            # Register all connecting variables before equations
            for component in registered_components:
                component.register_connecting(**kwargs)
            for component in registered_components:
                component.register_equations(**kwargs)

        def get_model(self):
            """Returns the gekko model.

            Returns
            -------
            m : GEKKO
            """
            return self._m

        def register_connecting(self, name, type, value, **kwargs):
            """Registers a connecting variable with the model, if not already
            registered. If registered, but type is `m.SV`, this will overwrite
            the old registration.

            This variable will be available at `self.name`.

            Parameters
            ----------
            name : str
                The name of the connecting variable. Do not prefix this name
                with an '_'
            type : m.Var, m.MV, etc
                The type of variable to register.
            value : obj
                The value to assign to the connecting variable.
            kwargs : key word arguments
                Any additional arguments to pass in to the variable
                constructor.

            Returns
            -------
            var : type
                A handle to the connecting variable, whether it was created or
                already exists.
            """
            if not hasattr(self, name) or type == self._m.SV:
                setattr(self, name, type(value=value, **kwargs))

            return getattr(self, name)

        def solve(self, glamdring=False, imode=None, solver=3, **kwargs):
            """Wrapper around `m.solve()`.

            If glamdring=True, will connect to the idealabs glambdring machine
            and solve there. Otherwise, it will solve given the options.

            Parameters
            ----------
            glamdring : bool, default=False
                Whether or not to run on glamdring.
            imode : int or None
                The simulation mode, or None if pre-specified.
            solver : int
                The solver to use.
            **kwargs : key word arguments
                Arguments to pass into `m.solve()` (e.g. disp=False).
            """
            m = self._m
            m.options.solver = solver

            # Set glamdring options if specified
            if glamdring:
                # m.server = 'http://192.168.17.19:8341'
                m.server = 'https://idealabs.byu.edu/apm'
                m.solver_options = ['linear_solver ma97']
                kwargs['remote'] = True

            # Set imode options if specified
            if imode is not None:
                m.options.IMODE = imode

                if imode == 7:
                    for var in m.variables:
                        if hasattr(var, 'FSTATUS'):
                            var.FSTATUS = 1
                        var.value.change = False

                    m.options.CSV_READ = 1

            m.solve(**kwargs)
