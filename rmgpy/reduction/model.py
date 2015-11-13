from collections import Counter

class ReductionReaction(object):
    """
    A class that enhances RMG-Py's  Reaction object
    by providing storage for the forward (kf) and backward
    (kb) rate coefficient.

    Once k is computed, it is stored and fetched
    when requested.

    """
    def __init__(self, rmg_reaction):
        super(ReductionReaction, self).__init__()
        self.rmg_reaction = rmg_reaction
        self.reactants = rmg_reaction.reactants
        self.products = rmg_reaction.products
        self.kf = None
        self.kb = None
        self.stoichio = {}
        self.create_stoichio()
    
    def __str__(self):
        return str(self.rmg_reaction)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (self.__class__, (self.rmg_reaction, ))


    def getRateCoefficient(self, T,P):
        if self.kf is None:
            self.kf = self.rmg_reaction.getRateCoefficient(T,P)
            return self.kf
        else: return self.kf
    
    def getReverseRateCoefficient(self, T, P):
        if self.kb is None:
            kf = self.getRateCoefficient(T,P) 
            self.kb = kf / self.rmg_reaction.getEquilibriumConstant(T)
            return self.kb
        else: return self.kb
    
    def create_stoichio(self):
        c_reactants = Counter([mol.label for mol in self.reactants])
        self.stoichio['reactant'] = c_reactants

        c_products = Counter([mol.label for mol in self.products])
        self.stoichio['product'] = c_products

    def get_stoichiometric_coefficient(self, spc_i, reactant_or_product):       
        return self.stoichio[reactant_or_product][spc_i.label]