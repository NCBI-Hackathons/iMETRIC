import mhcpredict.tools

def predictPeptides(sequences, alleles=None, species=None, 
                    methods=None, config=None, **kwargs):
    
    predictors = mhcpredict.tools.get_predictors(methods, config)
    results = dict((name, pred.predictPeptides(sequences, alleles, species, **kwargs))
        for name, pred in predictors.items())
    return results
    # TODO: generate consensus table from results

def predictProteins(sequences, lengths=None, alleles=None, species=None, 
                    methods=None, config=None, **kwarg):
    
    predictors = mhcpredict.tools.get_predictors(methods, config)
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
        
        self.getPeptidePredictions(sequences, alleles, species, **kwargs)
    
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
        
        self.getProteinPredictions(sequences, lengths, alleles, species, **kwargs)
    
    def _validate_args(self, sequences, alleles, species, min_seq_len=None, max_seq_len=None):
        if sequences is None or len(sequences) == 0:
            raise Exception("No sequences given")
        if isinstance(sequences, str):
            sequences = [sequences]
        for seq in sequences:
            if ((min_seq_len is None or len(seq) < min_seq_len) and
                    (max_seq_len is None or len(seq) > max_seq_len)):
                raise Exception("Sequence length must be between {0} and {1}".format(
                    min_seq_len or 0, max_seq_len or "Inf"))
        
        if alleles is None:
            alleles = self.getAllMHCAlleles()
        elif isinstance(alleles, str):
            alleles = [alleles]
        
        if species is None:
            species = self.getAllSpecies()
        elif isinstance(species, str):
            species = [species]

        return (sequences, alleles, species)
    
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
