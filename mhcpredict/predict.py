import mhcpredict.tools

def predictPeptides(self, sequences, alleles=None, species=None, 
                    methods=None, config=None, **kwarg):
    
    predictors = mhcpredict.tools.get_predictors(methods, config)
    results = (pred.predictPeptides(sequences, alleles, species, **kwargs)
        for pred in predictors)
    # TODO: generate consensus table from results

def predictProteins(self, sequence, lengths=None, alleles=None, species=None, 
                    methods=None, config=None, **kwarg):
    
    predictors = mhcpredict.tools.get_predictors(methods, config)
    results = (pred.predictProteins(sequences, lengths, alleles, species, **kwargs)
        for pred in predictors)
    # TODO: generate consensus table from results

MIN_PEPTIDE_LENGTH = 8
MAX_PEPTIDE_LENGTH = 15

class MHCPeptidePredictor(object):
    """Base class for tools that predict MHC:peptide binding strength.
    """
    def __init__(self, **kwargs):
        self.all_alleles = None
        self.all_species = None
        self.init(**kwargs)
    
    def predictPeptides(self, sequences, alleles=None, species=None, **kwarg):
        """Predict multiple peptides. By default, this just calls predictPeptide
        for each sequence, but subclasses may override to implement a more 
        efficient batch method.
        """
        return list(self.predictPeptide(seq, alleles, species, **kwarg)
            for seq in sequences)
    
    def predictPeptide(self, sequence, alleles=None, species=None, **kwarg):
        """Predict binding between a single peptide and one or more MHC alleles.
        
        Args:
            peptide: The peptide amino acid sequence.
            alleles: List of MHC alleles, or None if all alleles should be queried.
            species: Species name, or None if all species should be queried.
                     TODO: what to use for species names?
            kwargs: Additional arguments that are passed through to subclasses.
        
        Returns:
            A pandas DataFrame with the following columns:
            TODO
        """
        if (sequence is None or 
                len(sequence) < MIN_PEPTIDE_LENGTH or 
                len(sequence) > MAX_PEPTIDE_LENGTH):
            raise Exception("Peptide length must be between {0} and {1}".format(
                MIN_PEPTIDE_LENGTH, MAX_PEPTIDE_LENGTH))
        
        if alleles is None:
            alleles = self.getAllMHCAlleles()
        if isinstance(alleles, str):
            alleles = [alleles]
        
        if species is None:
            species = self.getAllSpecies()
        if isinstance(species, str):
            species = [species]
        
        self.getPeptidePredictions(sequence, alleles, species, **kwargs)
    
    def predictProteins(self, sequences, lengths=None, alleles=None, species=None, **kwarg):
        """Predict multiple proteins. By default, this just calls predictProtein
        for each sequence, but subclasses may override to implement a more 
        efficient batch method.
        """
        return list(self.predictProtein(seq, lengths, alleles, species, **kwarg)
            for seq in sequences)
    
    def predictProtein(self, sequence, lengths=None, alleles=None, species=None, **kwarg):
        """Predict binding between peptides within a protein sequence and one 
        or more MHC alleles. Each tool provides it's own method for deriving
        peptides from a protein; tools that do not provide such ability will
        raise an Exception.
        
        Args:
            sequence: The protein amino acid sequence.
            lengths: List of peptide lengths for which to make prediction.
            alleles: List of MHC alleles, or None if all alleles should be queried.
            species: Species name, or None if all species should be queried.
                     TODO: what to use for species names?
            kwargs: Additional arguments that are passed through to subclasses.
        
        Returns:
            A pandas DataFrame with the following columns:
            TODO
        """
        if sequence is None or len(sequence) < MIN_PEPTIDE_LENGTH:
            raise Exception("Peptide length must be between {0} and {1}".format(
                MIN_PEPTIDE_LENGTH, MAX_PEPTIDE_LENGTH))
            
        if alleles is None:
            alleles = self.getAllMHCAlleles()
        if isinstance(alleles, str):
            alleles = [alleles]
        
        if species is None:
            species = self.getAllSpecies()
        if isinstance(species, str):
            species = [species]
        
        if lengths is None:
            lengths = range(8, 15)
        elif isinstance(lengths, int):
            lengths = [lengths]
        
        self.getProteinPredictions(sequence, lengths, alleles, species, **kwargs)
    
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
    
    def getPeptidePredictions(sequence, alleles, species, **kwargs):
        raise NotImplemented()
    
    def getProteinPredictions(sequence, lengths, alleles, species, **kwargs):
        raise NotImplemented()
    
    def listMHCAlleles(self):
        raise NotImplemented()
    
    def listSpecies(self):
        raise NotImplemented()
