import configparser
import glob
import importlib
import os
import sys

import mhcpredict
import mhcpredict.util

def get_predictors(names=None, config=None):
    return ToolLoader().get_predictors(names, config)

class ToolLoader(object):
    def __init__(self):
        self.all_MHCPeptide_modules = None
    
    def get_predictors(self, names=None, config=None):
        """Create a predictor for each name given, or create instances
        of all available predictors if 'names' is None.
    
        Args:
            names: list of predictor names.
            config: ConfigParser instance, or path to config file.
    
        Returns:
            A list of MHCPeptidePredictors
        """
        if config is None:
            config = mhcpredict.util.DictConfig()
        
        modules = self.get_MHCPeptide_modules()
    
        if names is None:
            names = modules.keys()
    
        elif isinstance(names, str):
            names = [names]
    
        def create_predictor(name):    
            if name in modules:
                mod = modules[name]
                return mod.get_instance(config)
        
            else:
                raise Exception("No such predictor {}".format(name))
    
        return dict((name, create_predictor(name)) for name in names)

    def get_MHCPeptide_modules(self):
        if self.all_MHCPeptide_modules is None:
            
            mod_dir = os.path.join(
                os.path.dirname(os.path.abspath(mhcpredict.__file__)), 
                "tools"
            )
            
            #if sys.version_info >= (3, 4):
            #    import importlib.util
            #    spec = importlib.util.find_spec("mhcpredict.tools")
            #    mod_dir = spec.submodule_search_locations[0]
            
            #elif sys.version_info < (3, 0):
            #    import imp
            #    file, mod_dir, description = imp.find_module("tools", )
                
            mod_names = set(map(
                lambda path: os.path.splitext(os.path.basename(path))[0], 
                glob.glob(os.path.join(mod_dir, "mp_*.py*"))
            ))
            self.all_MHCPeptide_modules = dict(
                (mod_name, importlib.import_module("."+mod_name, package="mhcpredict.tools"))
                for mod_name in mod_names)

        return self.all_MHCPeptide_modules