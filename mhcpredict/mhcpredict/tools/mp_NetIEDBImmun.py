# TODO: refactor to move code that is common between
# this and mp_NetMHCIIPan.py to a shared base class

from mhcpredict.predict import MHCPeptidePredictor
from mhcpredict.util import create_temp_fasta, sort_by_length
import os
import subprocess

def get_instance(config):
    #if IEDB locally installed
    return LocalNetMHCPanPredictor(**config.get_section("IEDBImmun"))
    #else
    #   return WebNetMHCIIPanPredictor(config)

class LocalIEDBImmunPredictor(MHCImmunoPredictor):
    """MHCImmunoPredictor that calls a locally installed 
    IEDB Immunonecity tool local instance."""
    
    def init(self, **kwargs):
        self.executable = kwargs.get("executable", "predict_immunogenicity.py")
        self.tempdir = kwargs.get("tempdir", None)


    def getEpitopePredictions(self, sequences):
        seq_lengths = sort_by_length(sequences)
        rows_list = []
#        for seq_len, seqs in seq_lengths.items():
            rows_list.extend(self._predict(seqs))
        return self._prepare_DataFrame(rows_list)



#    def getProteinPredictions(self, sequences, lengths, alleles, species):
#        rows_list = self._predict(sequences, lengths, alleles, species)
#        return self._prepare_DataFrame(rows_list)
    
 
    def _predict(self, sequences):
#        alleles = list(allele.split("-")[1].replace("*", "_") for allele in alleles)
#        lengths = ",".join(map(str, lengths))
        seq_file = create_temp_fasta(sequences, self.tempdir)

        try:
            return list(self._execute(seq_file))


        finally:
            os.remove(seq_file)


    def _predict(self, sequences):
#        alleles = list(allele.split("-")[1].replace("*", "_") for allele in alleles)
#        lengths = ",".join(map(str, lengths))
        seq_file = create_temp_fasta(sequences, self.tempdir)

        try:
            return list(self._execute(seq_file))


        finally:
            os.remove(seq_file)




    def _execute(self, seq_file, lengths_str, allele):
        cmd = [
            self.executable, 
            seq_file
        ]
        print(cmd)
        output = subprocess.check_output(cmd)
        output = output.split("\n")[1:]
        return output
    
    def _prepare_DataFrame(self, rows_list):
        df = rbind(rows_list)
        df.colnames = [
            "peptide", "length", "score"
        ]            
        df = df.convert_objects(convert_numeric=True)
        df = df.drop(["length"], 1)
        df = df.dropna()
    
    def listMHCAlleles(self):
        """Get available alleles"""
        cmd = [self.executable, "-list"]
        temp = subprocess.check_output(cmd)
        alleles = temp.split("\n")[34:]
        return alleles
