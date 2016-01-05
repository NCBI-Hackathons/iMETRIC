# TODO: refactor to move code that is common between
# this and mp_NetMHCIIPan.py to a shared base class

import csv
import os
import re
import subprocess
import pandas as pd
from mhcpredict.predict import MHCPeptidePredictor
from mhcpredict.util import create_temp_fasta, sort_by_length, rbind

def get_instance(config):
    #if IEDB locally installed
    return LocalNetMHCPanPredictor(**config.get_section("NetMHCpan"))
    #else
    #   return WebNetMHCIIPanPredictor(config)

class LocalNetMHCPanPredictor(MHCPeptidePredictor):
    """MHCPeptidePredictor that calls a locally installed 
    netMHCpan instance."""
    
    def init(self, **kwargs):
        self.executable = kwargs.get("executable", "netMHCpan")
        self.tempdir = kwargs.get("tempdir", None)

    def getPeptidePredictions(self, sequences, alleles, species):
        if len(sequences) == 0 or len(alleles) == 0:
            # TODO: warn
            return None
        seq_lengths = sort_by_length(sequences)
        results = []
        for seq_len, seqs in seq_lengths.items():
            results.extend(self._predict(seqs, [seq_len], alleles, species))
        return self._prepare_DataFrame(list(result[0] for result in results))
    
    def getProteinPredictions(self, sequences, lengths, alleles, species):
        if len(sequences) == 0 or len(alleles) == 0:
            # TODO: warn
            return None
        results = self._predict(sequences, lengths, alleles, species)
        return self._prepare_DataFrame(list(result[0] for result in results))
    
    def _predict(self, sequences, lengths, alleles, species):
        alleles = list(allele.replace("*", "_") for allele in alleles)
        lengths = ",".join(map(str, lengths))
        seq_file = create_temp_fasta(sequences, self.tempdir)
        
        try:
            return list(self._execute(seq_file, lengths, allele)
                for allele in alleles)
        finally:
            #os.remove(seq_file)
            pass
        
    def _execute(self, seq_file, lengths_str, allele):
        cmd = [
            self.executable, 
            "-l", lengths_str, 
            "-a", allele, 
            "-f", seq_file,
            "-tdir", self.tempdir
        ]
        output = subprocess.check_output(cmd)

        accuracies = []
        summaries = []
        rows = []
        div_count = 0
        
        for row in output.split("\n"):
            if row.startswith("#"):
                continue
            elif row.startswith("-"):
                # we found a table divider row
                div_count += 1
            elif div_count == 0:
                # we're before the table
                accuracies.append(row)
            elif div_count == 1:
                # this is the header - ignore
                pass
            elif div_count == 2:
                # we're in the main body of the table
                rows.append(row)
            else:
                # we're after the table
                summaries.append(row)
        
        return (
            list(row[0:7] for row in
                csv.reader(rows, delimiter=" ", skipinitialspace=True)),
            accuracies, 
            summaries
        )
    
    def _prepare_DataFrame(self, rows_list):
        df = rbind(rows_list)
        df.columns = [
            "pos", "allele", "peptide", "identity",
            "1-log50k(aff)", "affinity", "rank"
        ]
        df = df.drop(["identity", "rank"], 1)
        df["pos"] = pd.to_numeric(df["pos"], errors='coerce')
        df["1-log50k(aff)"] = pd.to_numeric(df["1-log50k(aff)"], errors='coerce')
        df["affinity"] = pd.to_numeric(df["affinity"], errors='coerce')
        df = df.dropna()
        df["rank"] = df["affinity"].rank(method="min", ascending=1)
        return df
    
    def listMHCAlleles(self):
        """Get available alleles"""
        cmd = [self.executable, "-listMHC"]
        temp = subprocess.check_output(cmd)
        alleles = temp.split("\n")[34:]
        return alleles
