from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
from collections import defaultdict
import tempfile

def create_temp_fasta(sequences, tempdir=None):
    if isinstance(sequences, str):
        sequences=[sequences]
    temp = tempfile.mkstemp(dir=tempdir)[1]
    with open(temp, "w") as out:
        for i, seq in enumerate(sequences):
            SeqIO.write(
                SeqRecord(Seq(seq), id="temp{}".format(i), description=""), 
                out, "fasta")
    return temp

def sort_by_length(seq):
    d = defaultdict(lambda: [])
    for s in seq:
        d[len(s)].append(s)
    return d

def rbind(dfs):
    df = None
    for d in dfs:
        if d is None:
            continue
        if not isinstance(d, pd.DataFrame):
            d = pd.DataFrame(d)
        if df is None:
            df = d
        else:
            df.append(d, ignore_index=True)
    return df

def load_config(config):
    if isinstance(config, str):
        config = configparser.SafeConfigParser()
        config.read(config)
    if isinstance(config, dict):
        config = DictConfig(config)
    else:
        config = CPConfig(config)
    return config

class CPConfig(object):
    """Wrapper for a ConfigParser"""
    def __init__(self, config):
        self.config = config
    
    def get_section(self, section):
        if self.config.has_section(section):
            return dict(self.config.items(section))
        else:
            return {}
    
    def get(self, section, option, default=None):
        if self.config.has_option(section, option):
            return self.config.get(section, option)
        else:
            return default

class DictConfig(object):
    def __init__(self, config={}):
        self.config = config
    
    def get_section(self, section):
        return self.config.get(section, {})
        
    def get(self, section, option, default=None):
        if section in self.config:
            return self.config[section].get(option, default)
