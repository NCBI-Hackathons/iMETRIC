# TODO: refactor to move code that is common between
# this and mp_NetMHCIIPan.py to a shared base class

from mhcpredict.predict import MHCPeptidePredictor
from mhcpredict.util import create_temp_fasta, sort_by_length
import os
import subprocess

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
        print(seq_file)
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
        output = output.split("\n")[19:]
        ignore = set(("Protein","pos",""))

        def parse_row(row):
            if not row.startswith("-"):
                row = re.split("\s*", row.strip())[:9]
                if len(row) == 9 and row[0] not in ignore:
                    return row
        
        return filter(map(parse_row, output))
    
    def _prepare_DataFrame(self, rows_list):
        df = rbind(rows_list)
        df.colnames = [
            "pos", "allele", "peptide", "identity", "pos", 
            "core", "1-log50k(aff)", "affinity", "rank"
        ]
        df = df.convert_objects(convert_numeric=True)
        df = df.drop(["pos", "identity", "rank"], 1)
        df = df.dropna()
        df["rank"] = df["affinity"].rank(method="min", ascending=1)
        return df
    
    def listMHCAlleles(self):
        """Get available alleles"""
        cmd = [self.executable, "-list"]
        temp = subprocess.check_output(cmd)
        alleles = temp.split("\n")[34:]
        return alleles
