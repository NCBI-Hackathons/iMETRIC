import mhcpredict.tools

def predictPeptides(sequences, alleles=None, species=None, 
                    methods=None, config=None, **kwargs):
    
    predictors = mhcpredict.tools.get_MHCPeptidePredictors(methods, config)
    results = dict((name, pred.predictPeptides(sequences, alleles, species, **kwargs))
        for name, pred in predictors.items())
    return results
    # TODO: generate consensus table from results

def predictProteins(sequences, lengths=None, alleles=None, species=None, 
                    methods=None, config=None, **kwarg):
    
    predictors = mhcpredict.tools.get_MHCPeptidePredictors(methods, config)
    results = dict((name, pred.predictProteins(sequences, lengths, alleles, species, **kwargs))
        for name, pred in predictors.items())
    return results
    # TODO: generate consensus table from results

class MHCPeptidePredictor(object):
    """Base class for tools that predict MHC:peptide binding strength.
    """
    def __init__(self, **kwargs):
        # set defaults that can be over-ridden by subclass init()
        self.all_alleles = None
        self.all_species = None
        self.min_peptide_length = 8
        self.max_peptide_length = 15
        self.init(**kwargs)
    
    def predictPeptides(self, sequences, alleles=None, species=None, **kwargs):
        """Predict binding between one or more peptide and one or more MHC alleles.
        
        Args:
            peptide: The peptide amino acid sequences. May be of different lengths.
            alleles: List of MHC alleles, or None if all alleles should be queried.
            species: Species name, or None if all species should be queried.
                     TODO: what to use for species names?
            kwargs: Additional arguments that are passed through to subclasses.
        
        Returns:
            A pandas DataFrame with the following columns:
            TODO
        """
        sequences, alleles, species = self._validate_args(sequences, alleles, species,
            self.min_peptide_length, self.max_peptide_length)
        
        return self.getPeptidePredictions(sequences, alleles, species, **kwargs)
    
    def predictProteins(self, sequences, lengths=None, alleles=None, species=None, **kwarg):
        """Predict binding between peptides within a protein sequence and one 
        or more MHC alleles. Each tool provides it's own method for deriving
        peptides from a protein; tools that do not provide such ability will
        raise an Exception.
        
        Args:
            sequence: The protein amino acid sequences.
            lengths: List of peptide lengths for which to make prediction.
            alleles: List of MHC alleles, or None if all alleles should be queried.
            species: Species name, or None if all species should be queried.
                     TODO: what to use for species names?
            kwargs: Additional arguments that are passed through to subclasses.
        
        Returns:
            A pandas DataFrame with the following columns:
            TODO
        """
        sequences, alleles, species = self._validate_args(sequences, alleles, species,
            self.min_peptide_length)
        
        if lengths is None:
            lengths = range(self.min_peptide_length, self.max_peptide_length)
        
        else:
            if isinstance(lengths, int):
                lengths = [lengths]
            
            for l in lengths:
                if l < self.min_peptide_length or l > self.max_peptide_length:
                    raise Exception("All lengths must be between {0} and {1}".format(
                        self.min_peptide_length, self.max_peptide_length))
        
        return self.getProteinPredictions(sequences, lengths, alleles, species, **kwargs)
    
    def _validate_args(self, sequences, alleles, species, min_seq_len=None, max_seq_len=None):
        if sequences is None or len(sequences) == 0:
            raise Exception("No sequences given")
        if isinstance(sequences, str):
            sequences = [sequences]
        
        valid_sequences = []
        for seq in sequences:
            if ((min_seq_len is None or len(seq) >= min_seq_len) and
                    (max_seq_len is None or len(seq) <= max_seq_len)):
                valid_sequences.append(seq)
            else:
                # TODO: warn about this
                pass
        
        all_alleles = self.getAllMHCAlleles()
        if alleles is None:
            valid_alleles = all_alleles
        else:
            if isinstance(alleles, str):
                alleles = [alleles]
            
            valid_alleles = set(alleles) & set(all_alleles)
            invalid_alleles = set(alleles) - valid_alleles
            if len(invalid_alleles) > 0:
                # TODO: warn about this
                pass
            valid_alleles = list(valid_alleles)
        
        all_species = self.getAllSpecies()
        if species is None:
            valid_species = all_species
        else:
            if isinstance(species, str):
                species = [species]
            
            valid_species = set(species) & set(all_species)
            invalid_species = set(species) - valid_species
            if len(invalid_species) > 0:
                # TODO: warn about this
                pass
            valid_species = list(valid_species)

        return (valid_sequences, valid_alleles, valid_species)
    
    def getAllMHCAlleles(self):
        """Enumerate all alleles supported by the predictor.
        
        Returns:
            A list of MHC allele names.
        """
        if self.all_alleles is None:
            self.all_alleles = self.listMHCAlleles()
        
        return self.all_alleles
    
    def getAllSpecies(self):
        """Enumerate all species supported by the predictor.
        
        Returns:
            A list of species identifiers.
        """
        if self.all_species is None:
            self.all_species = self.listSpecies()
        
        return self.all_species
    
    ## Internal methods to be implemented by subclasses ##
    
    def init(self, **kwargs):
        pass
    
    def getPeptidePredictions(self, sequences, alleles, species, **kwargs):
        raise NotImplemented()
    
    def getProteinPredictions(self, sequences, lengths, alleles, species, **kwargs):
        raise NotImplemented()
    
    def listMHCAlleles(self):
        raise NotImplemented()
    
    def listSpecies(self):
        """List all species supported by the predictor. Returns an
        empty list by default, since some tools do not require species
        information."""
        return []

def predictEpitopes(sequences, methods=None, config=None, **kwargs):
    predictors = mhcpredict.tools.get_MHCImmunoPredictors(methods, config)
    results = dict((name, pred.predictEpitopes(sequences, **kwargs))
        for name, pred in predictors.items())
    return results
    # TODO: generate consensus table from results

class MHCImmunoPredictor(object):
    def predictEpitopes(self, sequences, **kwargs):
        """Predict immunogenicity of one or more epitopes.
        
        Args:
            sequences: The peptide amino acid sequences. May be of different lengths.
            kwargs: Additional arguments that are passed through to subclasses.
        
        Returns:
            A pandas DataFrame with the following columns:
            TODO
        """
        sequences = self._validate_args(sequences)
        self.getEpitopePredictions(sequences)
    
    def _validate_args(self, sequences, alleles, species, min_seq_len=None, max_seq_len=None):
        if sequences is None or len(sequences) == 0:
            raise Exception("No sequences given")
        if isinstance(sequences, str):
            sequences = [sequences]
        return sequences
                    
    ## Internal methods to be implemented by subclasses ##
    def getEpitopePredictions(self, sequences, **kwargs):
        raise NotImplemented()