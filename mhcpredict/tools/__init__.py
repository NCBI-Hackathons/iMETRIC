import ConfigParser
import importlib
import importlib.machinery
import glob

def get_predictors(names=None, config=None):
    """Create a predictor for each name given, or create instances
    of all available predictors if 'names' is None.
    
    Args:
        names: list of predictor names.
        config: ConfigParser instance, or path to config file.
    
    Returns:
        A list of MHCPeptidePredictors
    """
    if isinstance(config, str):
        config = ConfigParser.SafeConfigParser()
        config.read(config)
        
    modules = get_MHCPeptide_modules()
    
    if names is None:
        names = modules.keys()
    
    elif isinstance(names, str):
        names = [names]
    
    def create_predictor(name):    
        if name in modules:
            mod = modules[name]
            cls = mod.predictor_class
            return cls(config)
        
        else:
            raise Exception("No such predictor {}".format(name))
    
    return map(create_predictor, names)

## Internal ##

all_MHCPeptide_modules = None

def get_MHCPeptide_modules():
    if all_MHCPeptide_modules is None:
        spec = importlib.machinery.PathFinder.find_spec("mhcpredict.tools")
        mod_dir = spec.submodule_search_locations[0]
        mod_names = set(map(
            lambda path: os.path.splitext(os.path.basename(path))[0], 
            glob.glob("mp_*.py*")
        ))
        all_MHCPeptide_modules = dict(
            (mod_name, importlib.import_module(mod_name, package="mhcpredict.tools"))
            for mod_name in mod_names)

    return all_MHCPeptide_modules
