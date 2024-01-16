"""Classes to import, handle, manipulate and save different objects etc.
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

from Bio.KEGG import REST, Enzyme
import cobra
import os.path
import urllib.error
import warnings

################################################################################
# functions
################################################################################

# @DELETE
# added to refinegems
def read_model_cobra(path: str):
    """Reads in a model, depending on the given file type.

    :param path: Path to the model.
    :type  path: string

    :raises:  :class:`ValueError`: Unknown file extension for model.

    :returns: The model read in with cobra.
    :rtype: cobra.Model
    """

    extension = os.path.splitext(path)[1].replace('.','')

    match extension:
        case 'xml':
            data = cobra.io.read_sbml_model(path)
        case 'json':
            data = cobra.io.load_json_model(path)
        case 'yml':
            data = cobra.io.load_yaml_model(path)
        case 'mat':
            data = cobra.io.load_matlab_model(path)
        case _:
            raise ValueError('Unknown file extension for model: ', extension)
            sys.exit(1)

    return data


#@TODO
def kegg_reaction_parser(rn_id):
    """Get the KEGG reaction entry for a KEGG reaction ID
    and parse it to retrieve information in form of a dictionary.

    :param rn_id: The KEGG reaction ID that shall be retrived and parsed.
    :type  rn_id: string
    :returns:     Information about name, equationand database links of the entry (reaction).
    :rtype:       dict, keys include 'name','db','rc' and 'equation'
    """

    # get KEGG reaction entry
    # .................
    #@TODO
    #     implement try-catch for errors (just in case)
    #     although they (theoreticlly) should not occur
    # current solution: do it outside this function
    # .................
    kegg_reac = REST.kegg_get(F'rn:{rn_id}')
    kegg_reac = kegg_reac.read()

    # parse the entry for necessary information
    features = {}
    db_entries = []
    pathways = []
    rc = []
    references = {}
    collect = False
    for line in kegg_reac.split('\n'):
        if line:
            if line.startswith('NAME'):
                features['name'] = line.replace('NAME','',1).strip()
            elif line.startswith('EQUATION'):
                features['equation'] = line.replace('EQUATION','',1).strip()
            elif line.startswith('ENZYME'):
                references['ec-code'] = line.replace('ENZYME','',1).strip()
            elif line.startswith('RCLASS'):
                rc.append(line.replace('RCLASS','',1).strip().split(' ')[0])
                collect = True
            elif line.startswith('PATHWAY'):
                pathways.append(line.replace('PATHWAY','',1).strip().split(' ')[0])
                collect = True
            elif line.startswith('DBLINKS'):
                db_entries.append(line.replace('DBLINKS','',1).strip())
                collect = True
            elif collect == True and line[0] != '/':
                if len(db_entries) == 0:
                    if line[0].isupper():
                        collect = False
                    else:
                        line = line.strip()
                        if line.startswith('RC'):
                            rc.append(line.split(' ')[0])
                        else:
                            pathways.append(line.split(' ')[0])
                else:
                    db_entries.append(line.strip())
            else:
                continue

    # parse references
    for entry in db_entries:
        db, identifier = entry.split(':')
        db = db.strip().lower()
        if db in references:
            references[db] = references[db].append(identifier)
        else:
            references[db] = [identifier.strip()]
    if len(pathways) > 0:
        references['kegg.pathway'] = pathways

    if len(references) > 0:
        features['db'] = references
    if len(rc) > 0:
        features['rc'] = rc

    return features
