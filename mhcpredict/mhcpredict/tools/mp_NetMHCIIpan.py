from mhcpredict.predict import MHCPeptidePredictor
from mhcpredict.util import create_temp_fasta, sort_by_length, rbind
import os
import re
import subprocess

def get_instance(config):
    #if IEDB locally installed
    return LocalNetMHCIIPanPredictor(**config.get_section("NetMHCIIpan"))
    #else
    #   return WebNetMHCIIPanPredictor(config)

class LocalNetMHCIIPanPredictor(MHCPeptidePredictor):
    """MHCPeptidePredictor that calls a locally installed 
    netMHCIIpan instance."""
    
    def init(self, **kwargs):
        self.executable = kwargs.get("executable", "netMHCIIpan")
        self.tempdir = kwargs.get("tempdir", None)

    def getPeptidePredictions(self, sequences, alleles, species):
        if len(sequences) == 0 or len(alleles) == 0:
            # TODO: warn
            return None
        seq_lengths = sort_by_length(sequences)
        rows_list = []
        for seq_len, seqs in seq_lengths.items():
            rows_list.extend(self._predict(seqs, [seq_len], alleles, species))
        return self._prepare_DataFrame(rows_list)
    
    def getProteinPredictions(self, sequences, lengths, alleles, species):
        if len(sequences) == 0 or len(alleles) == 0:
            # TODO: warn
            return None
        rows_list = self._predict(sequences, lengths, alleles, species)
        return self._prepare_DataFrame(rows_list)
    
    def _predict(self, sequences, lengths, alleles, species):
        alleles = list(allele.split("-")[1].replace("*", "_") for allele in alleles)
        lengths = ",".join(map(str, lengths))
        seq_file = create_temp_fasta(sequences, self.tempdir)

        try:
            return list(self._execute(seq_file, lengths, allele)
                for allele in alleles)

        finally:
            os.remove(seq_file)
        
    def _execute(self, seq_file, lengths_str, allele):
        cmd = [
            self.executable, 
            "-length", lengths_str, 
            "-a", allele, 
            "-f", seq_file,
            "-tdir", self.tempdir
        ]
        output = subprocess.check_output(cmd)
        output = output.split("\n")[19:]
        ignore = set(("Protein","pos",""))

        def parse_row(row):
            if not row.startswith("-"):
                row = re.split("\s*", row.strip())[:9]
                if len(row) == 9 and row[0] not in ignore:
                    return row
        
        return filter(None, map(parse_row, output))
    
    def _prepare_DataFrame(self, rows_list):
        df = rbind(rows_list)
        df.colnames = [
            "pos", "allele", "peptide", "identity", "pos", 
            "core", "1-log50k(aff)", "affinity", "rank"
        ]
        df = df.drop(["pos", "identity", "rank"], 1)
        pd.to_numeric(df[:,"pos"], errors='coerce')
        pd.to_numeric(df[:,"1-log50k(aff)"], errors='coerce')
        pd.to_numeric(df[:,"affinity"], errors='coerce')
        df = df.dropna()
        df["rank"] = df["affinity"].rank(method="min", ascending=1)
        return df
    
    def listMHCAlleles(self):
        """Get available alleles"""
        cmd = [self.executable, "-list"]
        temp = subprocess.check_output(cmd)
        alleles = temp.split("\n")[34:]
        return alleles
