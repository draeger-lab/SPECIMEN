"""Classes to import, handle, manipulate and save media.
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

import cobra
import copy
from importlib.resources import files
import pandas as pd
import re
import warnings

################################################################################
# variables (part 1)
################################################################################

PATH_TO_MEDIA_DB = files('medium').joinpath('media.csv')

################################################################################
# classes
################################################################################

class Compound:
    """The Compound objects describe the compounds of a Medium object.

    The name, BiGG ID, molecule formula and exchange_flux in the medium are saved as attributes
    of the Compound objects. While not directly needed for the object, when using a
    Compound object as a part of a medium object, at least the formula attribute
    should be set.

    :param name: The name of the compound.
    :type name: string
    :param bigg: The BiGG ID (without compartment e.g. '_e') of the compound.
    :type bigg: string
    :param formula: The molecule formula (without charges, e.g. NH4 for ammonia).
    :type formula: string
    :param exchange_flux: The exchange_flux of the compound in the medium.
    :type exchange_flux: double

    :var name: The name of the compound.
    :vartype name: string
    :var bigg: The BiGG ID (without compartment e.g. '_e') of the compound.
    :vartype bigg: string
    :var formula: The molecule formula (without charges, e.g. NH4 for ammonia).
    :vartype formula: string
    :var exchange_flux: The exchange_flux of the compound in the medium.
    :vartype exchange_flux: double
    """

    def __init__(self, name=None, bigg=None, formula=None, exchange_flux=None):
        self.name = name
        self.bigg = bigg
        self.formula = formula
        self.exchange_flux = exchange_flux


    def __eq__(self, other):

        if isinstance(other, self.__class__):
            if not (pd.isnull(other.name) and pd.isnull(other.name)) and not other.name == self.name:
                return False
            if not (pd.isnull(other.bigg) and pd.isnull(other.bigg)) and not other.bigg == self.bigg:
                return False
            if not (pd.isnull(other.formula) and pd.isnull(other.formula)) and not other.formula == self.formula:
                return False
            if not (pd.isnull(other.exchange_flux) and pd.isnull(other.exchange_flux)) and not other.exchange_flux == self.exchange_flux:
                return False
        return True


    def __ne__(self,other):
        return not self.__eq__(other)


    def get_elements(self):
        """Count the chemical elements of the Compound.

        :returns: A dictionary of the chemical elements and their counts.
        :rtype: dict
        """
        elements = {}

        if self.formula and not self.formula == '-' and not self.formula == '':

            symbol = None
            amount = ''

            for char in self.formula:
                if char.isupper():
                    if symbol:
                        if amount == '':
                            elements[symbol] = 1
                        else:
                            elements[symbol] = int(amount)
                    symbol = char
                    amount = ''
                elif char.isnumeric():
                    amount += char
                else:
                    symbol += char

            if symbol:
                if amount == '':
                    elements[symbol] = 1
                else:
                    elements[symbol] = int(amount)

        else:
            warnings.warn('No formula found, returning empty list.')

        return elements


    def bigg_to_ex(self):
        """Tranform the BiGG ID into the BiGG ID of corresponding the exchange reactions.
        If not BiGG ID is found, tries to use the formula instead.

        :raises: :class:`ValueError`: 'Both bigg and formula attributes are not set, conversion not possible.'

        :returns: The BiGG ID of the corresponding exchange reactions.
        :rtype: string
        """

        if not pd.isnull(self.bigg):
            return 'EX_' + self.bigg + '_e'
        elif not pd.isnull(self.formula):
            warnings.warn(F'No bigg id, using formula for conversion instead: {self.formula}')
            return 'EX_' + self.formula + '_e'
        else:
            raise ValueError('Both bigg and formula attributes are not set, conversion not possible.')


class Medium:
    """The Medium class objects contain the name and composition of a medium
    for genome-scale metabolic modelling.

    :param name: Name of the medium.
    :type  name: string
    :param compounds: Name (ideally BiGG identifier) and exchange_fluxs of the compounds of the medium.
    :type compounds: list of Compounds

    :var name: Name of the medium.
    :vartype name: string
    :var compounds: Name (ideally BiGG identifier) and exchange_fluxs of the compounds of the medium.
    :vartype compounds: dict, string (identifer) as key and float (exchange_flux in mM) as value.
    """

    def __init__(self, name=None, description=None, compounds=[]):
        self.name = name
        self.description = description
        self.compounds = {}
        for c in compounds:
            self.add_compound(c)


    def add_compound(self, c):
        """Add a compound to the medium.

        :param c: The compound to be added.
        :type c: Compound

        :raises: :class:`ValueError`: 'Neither bigg id or formula found. Cannot add to medium.'
        """
        if not pd.isnull(c.bigg) and not c.bigg == '' and not c.bigg == '-':
            self.compounds[c.bigg] = c
        else:
            if c.formula:
                self.compounds[c.formula] = c
            else:
                raise ValueError('Neither bigg id nor formula found. Cannot add to medium.')


    def has_compound(self, c):
        """Checks, if a compound is part of the medium.

        The check can either be performed with a Compound object or a name of a compound.

        :param c: Either name or Compound to search the Medium for.
        :type c: Compound or string

        :raises: :class:`TypeError`: ''Unknown type for c (compound)'
        """

        if isinstance(c, Compound):
            return c in self.compounds.values()
        elif isinstance(c, str):
            return c in self.compounds.keys()
        else:
            raise TypeError('Unknown type for c (compound): ',type(c))


    def get_compound_with_attribute(self, name=None, bigg=None, formula=None, exchange_flux=None):
        """Get all compounds that match certain attribute(s).

        Search for compounds that have the given name/bigg/formula/exchange_flux.
        Default for all attributes is None, which means a skip for this attribute

        :param name: The name of the compound.
        :type name: string
        :param bigg: The BiGG ID (without compartment e.g. '_e') of the compound.
        :type bigg: string
        :param formula: The molecule formula (without charges, e.g. NH4 for ammonia).
        :type formula: string
        :param exchange_flux: The exchange_flux of the compound in the medium.
        :type exchange_flux: double

        :returns: The compounds that match the given attributes (all!).
            Return an empty list if no match is found.
        :rtype: list, of Compound objects
        """

        same_name = None
        same_bigg = None
        same_formula = None
        same_exchange_flux = None

        if name:
            same_name = [_ for _ in self.compounds.values() if _.name == name]
        if bigg:
            same_bigg = [_ for _ in self.compounds.values() if _.bigg == bigg]
        if formula:
            same_formula = [_ for _ in self.compounds.values() if _.formula == formula]
        if exchange_flux:
            same_exchange_flux = [_ for _ in self.compounds.values() if _.exchange_flux == exchange_flux]

        intersec = ()
        test = [_ for _ in [same_bigg, same_name, same_formula, same_exchange_flux] if _]
        if len(test) > 1:
            intersec = test[0]
            for result in test[1:]:
                intersec = [_ for _ in intersec if _ in result]
            return intersec
        elif len(test) == 1:
            return test[0]
        else:
            return []


    def remove_compound(self, c):
        """Remove a compound from the medium.

        The removal can either be performed with a Compound object or the name of a compound.

        :param c: Either name or Compound to search the Medium for.
        :type c: Compound or string
        """

        if isinstance(c, Compound):
            if c in self.compounds.values():
                k = [_ for _ in self.compounds.keys() if self.compounds[_] == c][0]
                del self.compounds[k]
        elif isinstance(c, str):
            if c in self.compounds.keys():
                del self.compounds[c]


    def get_source_of(self, element):
        """Get the names of the compounds that are a source for a given element,
        meaning that the given element is found in the formula.

        :param element: The element to get the sources of.
        :type element: string

        :returns: The names of the compounds that are a source of the given element.
        :rtype: list of strings
        """

        # for all compounds
        #     get chemical formula
        #     get elements
        #     if element in elements -> source
        found = []
        for k,v in self.compounds.items():
            if element in v.get_elements().keys():
                found.append(k)

        return found


    def set_source_of(self, element, new_source, remove_all=False):
        """Set the source of an element to a given new source (Compound).

        This function replaces either all old sources with the new one (remove_all=True)
        or just the compounds that are not the exclusive source for any ofther element
        (remove_all=False).

        :param element: The chemical element.
        :type element: string
        :param new_source: The new source for the element.
        :type new_source: Compound
        :param remove_all: Option to eitherr remove all possible sources
            of the element or only those that serve NOT as the exclusive source
            for another element.
            Default is the latter option (set to False).
        :type remove_all: bool
        """

        # identify current source of element
        current = self.get_source_of(element)

        # if remove_all = False, check to see if a source from another
        # element is simultaniously removed
        #     if that is the case, do not remove that compound
        if not remove_all:
            tmp = current
            current = []
            # for every 'old' source of the given elements
            for comp in tmp:
                required = False
                comp_elements = self.compounds[comp].get_elements().keys()
                # check all elements
                for el in comp_elements:
                    # if its the elements to be replaced with a new source, skip
                    if el == element:
                        continue
                    # else get all possible sources for the other element
                    # and check if additional sources for that element exists
                    alternatives = self.get_source_of(el)
                    alternatives.remove(comp)
                    if not alternatives or all(_ in current for _ in alternatives):
                        required = True
                        break
                # only set for removal if another source for every other element
                # - save for the one to be replaced -
                # exists
                if not required:
                    current.append(comp)

        # delete old sources (based on remove_all rules)
        for comp in current:
            self.remove_compound(comp)

        # add new source
        self.add_compound(new_source)


    def is_aerobic(self):
        """Check if the medium is aerobic = contains O2.

        Checks the medium for the name 'o2' and the formula 'O2'

        :returns: The result of the check (True if medium is aerobic).
        :rytpe: bool
        """
        # if 'pure' O2 is in medium, it can be considered aerobic
        if self.has_compound('o2') or len(list(self.get_compound_with_attribute(formula='O2'))) > 0:
            return True
        else:
            return False


    def make_aerobic(self, c=10.0):
        """Make the medium aerobic (if it isn't already) by adding O2.

        :param c: The exchange_flux of O2.
            Default is 10.0
        :type c: double
        """

        # check if condition is already fullfilled
        if self.is_aerobic():
            return
        # else add elemental oxygen to medium
        else:
            self.add_compound(Compound(name='o2', bigg='o2', formula='O2', exchange_flux=c))
            return


    def make_anaerobic(self):
        """Make the medium anaerobic by removing elemental oxygen from it.

        Removes oxygen by searching for compounds with the formula 'O2'.

        @TODO
            add option:
            remove all additional oxygen sources, that are not needed for other elements?
        """
        if self.is_aerobic():
            oxygen = self.get_compound_with_attribute(formula='O2')
            while(len(oxygen) != 0):
                self.remove_compound(oxygen.pop())


    def combine(self,other):
        """Combine the medium with another medium.

        Add all elements of the second medium, which keys are not in the first
        medium to the first one.

        :returns: The combined medium (new object).
        :rtype: Medium
        """
        combined_medium = copy.deepcopy(self)
        combined_medium.name = self.name + '_' + other.name
        combined_medium.description = 'combined ' + self.name + ' + ' + other.name
        if isinstance(other, self.__class__):
            for c in other.compounds.keys():
                if c not in combined_medium.compounds.keys():
                    combined_medium.add_compound(copy.deepcopy(other.compounds[c]))
        else:
            raise TypeError("Parameter of type 'medium' expected, found", type(other))
        return combined_medium


    def __add__(self,other):
        return self.combine(other)


    def export_to_cobra(self, default_c=10.0):
        """Export the medium into a cobra model medium format.

        :returns: The medium in a format that can be directly set as a medium of a model.
        :rytpe: dict
        """
        cobra_medium = {}

        for compound in self.compounds.values():
            if not pd.isnull(compound.exchange_flux):
                cobra_medium[compound.bigg_to_ex()] = compound.exchange_flux
            else:
                cobra_medium[compound.bigg_to_ex()] = default_c

        return cobra_medium

    # @TODO
    def validate_for_model(self, model, add=False):
        """Validate a medium to be used for a cobra model.

        Checks, if all exchange reactions created from the medium are present in
        the model. If not,

        - add=False: removes corresponding compounds from the medium
        - @TODO add=True: add the missing reaction (and metabolite to the model)

        """
        if not add:
            ex_reac = self.export_to_cobra().keys()
            model_reac = [_.id for _ in model.reactions]
            for reac in ex_reac:
                if not reac in model_reac:
                    self.remove_compound(re.search('EX_(.*?)_e',reac).group(1))
        else:
            #@TODO
            #    add the exchange reactions and possible metabolites to model
            pass



################################################################################
# functions
################################################################################


# importing media
# ---------------

def from_table(table):
    """Construct a medium from a table.

    The table needs to have the following columns:
    medium,description,compound,bigg_id,formula, exchange_flux

    :param table: The table.
    :type table: pd.DataFrame

    :returns: The constructed medium.
    :rtype: Medium
    """
    # ...................................
    # @TODO:
    #    check if table format is correct
    # ...................................
    name = list(table['medium'])[0]
    description = list(table['description'])[0]
    compounds = table.apply(lambda row: Compound(row['compound'], row['bigg_id'], row['formula'], row['exchange_flux']), axis=1)
    return Medium(name, description, compounds)


def load_media_db(path=None):
    """Load a media database into a dictionary of media (name:medium).

    The database should be a csv file with ';' as seperators.
    The first line should be the following header:
    medium;description;compound;bigg_id;formula;exchange_flux

    :param path: The path to the database file.
        If None is given, reads the in-build media database.
    :type path: string

    :returns: The media found in the database as a dictionary of names and Media.
    :rtype: dict, {string:Medium}
    """
    if pd.isnull(path):
        path = PATH_TO_MEDIA_DB
    db = pd.read_csv(path, sep=';', header=0)
    media = {}
    for abbrev, medium in db.groupby('medium'):
        media[abbrev] = from_table(medium)
    return media


def import_medium_from_cobra(model):
    """Import a medium from a cobra model into a Medium ob object.

    :param model: The cobra model.
    :type model: cobra.Model

    :returns: The imported medium.
    :rtype: Medium
    """

    medium = Medium()
    medium.name = F'{model.id}_medium'
    medium.description = F'medium imported from cobra model {model.id}'
    for ex_name,c in model.medium.items():
        model_comp = list(model.reactions.get_by_id(ex_name).metabolites.keys())[0]
        if 'bigg.metabolite' in model_comp.annotation.keys():
            bigg = model_comp.annotation['bigg.metabolite']
        else:
            bigg = None
        if pd.isnull(model_comp.formula) or model_comp.formula == '*':
            warnings.warn(F'No formula or only wildcard found. Metabolite {model_comp.id} will not be added to medium.')
            continue
        comp = Compound(model_comp.id, bigg, model_comp.formula, c)
        medium.add_compound(comp)

    return medium


# ......
#@ TODO:
#    see function body
# ......
def model_minimal_medium(model, objective='flux', growth_rate=0.5, open_exchanges=False):
    """Get the minimal medium based on different objectives:
    - 'flux':      find the minimal fluxes based in current medium.
    - 'medium':    find the minimal number of compounds for the current medium.
    - 'exchanges': find the minimal number of compounds in a medium based on all avaiblae exchange reactions in the model.
    Note: there may be multiple solution for the minimisation, but only 1 will be returned

    :param model: The model to minimize the medium for.
    :type  model: cobra.Model
    :param objective: Objective for the minimisation task. Options listed above.
    :type  objective: string
    :param growth_rate: Minimum growth rate the model has to archieve.
    :type  growth_rate: double, only for objectives medium and exchanges
    :param open_exchanges: If set to True assigns large upper bound to all import reactions.
        Author's note: setting is to true leads to an endlessly running program for me.
    :type  open_exchanges: bool

    :raises: :class:`ValueError`: 'Unknown objective for minimisation.'

    :returns: The medium that is a solution for the minimisation task.
    :rtype: Medium
    """

    # ...............................
    # @TODO:
    #    make running multiple iterations possible
    # set iterations to True if no concrete number is given
    # PROBLEM: will leads to a different output
    #if not iterations:
    #    iterations = True
    # minimize_components = iteration
    # ...............................

    # minimise the fluxes od the current medium
    if objective == 'flux':
        max_growth = model.slim_optimize()
        min_medium = dict(cobra.medium.minimal_medium(model, max_growth))

    # minimise components of current medium
    elif objective == 'medium':
        warnings.warn('Warning: cobrapy.minimal_medium uses MIP formulation. This may take some time.')
        min_medium = dict(cobra.medium.minimal_medium(model, growth_rate, minimize_components=True, open_exchanges=open_exchanges))

    # get minimal number of medium components
    # based on exchange reaction possible in the model
    # note 1: can be time consuming
    # note 2: can lead to different results if run only once each time
    elif objective == 'exchanges':
        # create cobra medium from all available exchange reactions
        ex_medium = {_.id:1000.0 for _ in model.exchanges}
        # perform minimisation
        model.medium = ex_medium
        warnings.warn('Warning: cobrapy.minimal_medium uses MIP formulation. This may take some time.')
        min_medium = dict(cobra.medium.minimal_medium(model, growth_rate, minimize_components=True, open_exchanges=open_exchanges))

    else:
        raise ValueError('Unknown objective for minimisation.')

    # create a Medium object from the minimal medium
    with model as tmp_model:
        tmp_model.medium = min_medium
        medium = import_medium_from_cobra(tmp_model)

    return medium

# .........................................................................,
# @TODO
#    make an option that takes the current Medium
#    and expands using exchange reactions until the model growths (or all have been added)
#    problem: how to set the combination ?
# .........................................................................,
#def expand_medium_to_enable_growth():
#    pass

# save media
# ----------

def db_to_table(db):
    """Transform a medium database (dictionary of name:Medium) into a DataFrame.

    :param db: The medium dictionary.
    :type db: dict, {name:Medium}

    :returns: The table.
    :rtype: pd.DataFrame
    """
    for medium in db.values():
        compounds = []
        for compound in medium.compounds.values():
            compounds.append((medium.name,medium.description,compound.name,compound.bigg,compound.formula,compound.exchange_flux))
    table = pd.DataFrame(compounds,columns=['medium','description','compound','bigg_id','formula','exchange_flux'])

    return table


def save_db(db, out_path):
    """Save a medium database (dictionary of name:Medium) into a csv file
    (seperator is ';').

    :param db: The medium dictionary.
    :type db: dict, {name:Medium}
    :param out_path: Path to the file for the output.
    :type out_path: string
    """
    table = db_to_table(db)
    table.to_csv(out_path, sep=';',header=True, index=False)


def add_medium_to_db(medium, db_path):
    """Add a medium to a media database file.

    :param medium: The medium to be added.
    :type medium: Medium
    :param db_path: Path to the database.
    :type db_path: string
    """
    compounds = []
    for compound in medium.compounds.values():
        compounds.append((medium.name,medium.description,compound.name,compound.bigg,compound.formula,compound.exchange_flux))
    table = pd.DataFrame(compounds,columns=['medium','description','compound','bigg_id','formula','exchange_flux'])
    table.to_csv(db_path, mode='a', sep=';', header=False, index=False)


################################################################################
# variables (part 2)
################################################################################

CASAMINO_ACIDS = Medium(name='CAS',
                        description='Casamino acids, based in USBiological Life Sciences',
                        compounds=[Compound(name='L-Alanine', bigg='ala__L', formula='C3H7NO2'),
                                   Compound(name='L-Argenine', bigg='arg__L', formula='C6H15N4O2'),
                                   Compound(name='L-Aspartate', bigg='asp__L', formula='C4H6NO4'),
                                   Compound(name='L-Glutamate', bigg='glu__L', formula='C5H8NO4'),
                                   Compound(name='Glycine', bigg='gly', formula='C2H5NO2'),
                                   Compound(name='L-Histidine', bigg='his__L', formula='C6H9N3O2'),
                                   Compound(name='L-Isoleucine', bigg='ile__L', formula='C6H13NO2'),
                                   Compound(name='L-Leucine', bigg='leu__L', formula='C6H13NO2'),
                                   Compound(name='L-Lysine', bigg='lys__L', formula='C6H15N2O2'),
                                   Compound(name='L-Methionine', bigg='met__L', formula='C5H11NO2S'),
                                   Compound(name='L-Phenylalanine', bigg='phe__L', formula='C9H11NO2'),
                                   Compound(name='L-Proline', bigg='pro__L', formula='C5H9NO2'),
                                   Compound(name='L-Serine', bigg='ser__L', formula='C3H7NO3'),
                                   Compound(name='L-Threonine', bigg='thr__L', formula='C4H9NO3'),
                                   Compound(name='L-Tryptophan', bigg='trp__L', formula='C11H12N2O2'),
                                   Compound(name='L-Tyrosine', bigg='tyr__L', formula='C9H11NO3'),
                                   Compound(name='L-Valine', bigg='val__L', formula='C5H11NO2'),
                                   Compound(name='L-Cystine', bigg='cystine__L', formula='C6H12N2O4S2')])

AMINO_ACIDS = Medium(name = 'AA',
               description = 'the 20 proteinogenic amino acids',
               compounds=[Compound(name='L-Alanine', bigg='ala__L', formula='C3H7NO2'),
               Compound(name='L-Argenine', bigg='arg__L', formula='C6H15N4O2'),
               Compound(name='L-Asparagine', bigg='asn__L', formula='C4H8N2O3'),
               Compound(name='L-Aspartate', bigg='asp__L', formula='C4H6NO4'),
               Compound(name='L-Cysteine', bigg='cys__L', formula='C3H7NO2S'),
               Compound(name='L-Glutamine', bigg='gln__L', formula='C5H10N2O3'),
               Compound(name='L-Glutamate', bigg='glu__L', formula='C5H8NO4'),
               Compound(name='Glycine', bigg='gly', formula='C2H5NO2'),
               Compound(name='L-Histidine', bigg='his__L', formula='C6H9N3O2'),
               Compound(name='L-Isoleucine', bigg='ile__L', formula='C6H13NO2'),
               Compound(name='L-Leucine', bigg='leu__L', formula='C6H13NO2'),
               Compound(name='L-Lysine', bigg='lys__L', formula='C6H15N2O2'),
               Compound(name='L-Methionine', bigg='met__L', formula='C5H11NO2S'),
               Compound(name='L-Phenylalanine', bigg='phe__L', formula='C9H11NO2'),
               Compound(name='L-Proline', bigg='pro__L', formula='C5H9NO2'),
               Compound(name='L-Serine', bigg='ser__L', formula='C3H7NO3'),
               Compound(name='L-Threonine', bigg='thr__L', formula='C4H9NO3'),
               Compound(name='L-Tryptophan', bigg='trp__L', formula='C11H12N2O2'),
               Compound(name='L-Tyrosine', bigg='tyr__L', formula='C9H11NO3'),
               Compound(name='L-Valine', bigg='val__L', formula='C5H11NO2')])
