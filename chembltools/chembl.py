"""
module docstring
"""

from chembl_webresource_client.new_client import new_client


def inspect_synonyms(entry, term):
    """
    docstring
    """
    return term.lower() in {
        x["molecule_synonym"].lower() for x in entry["molecule_synonyms"]
    }


def get_chembl_id(compounds):
    """
    docstring
    """
    molecule = new_client.molecule
    ids = dict()
    if isinstance(compounds, str):
        compounds = [compounds]
    for compound in compounds:
        res = molecule.filter(pref_name__iexact=compound)
        if len(res) == 0:
            print("No exact match found for {}, trying a search for synonyms".format(compound))
            res = molecule.search(compound)
            res = [i for i in res if inspect_synonyms(i, compound)]
            if len(res) == 0:
                print("No match found for {}, skipping ...".format(compound))
                continue
        if len(res) > 1:
            print("Found {} exact matches for {}, picking first one".format(
                len(res), compound))
        ids[compound] = res[0]["molecule_chembl_id"]
    return ids


def get_similar_molecules(chembl_ids, similarity=90, show_similarity=False):
    """
    Get similar molecules based on a structural similarity search.
    Given a list of chembl ID's

    Note:
    -----
    Similarity matches may appear to be ordered by descening similarity score,
    though this is not guaranteed and should not be relied upon.

    Parameters:
    -----------
    chembl_ids: list-like
        CHEMBL ID codes

    similarity: numeric, between 0 and 100 (default = 90)
        minimum structural similarity threshold to determine a match

    show_similarity: Boolean (default = False)
        whether or not to return the similarity score alongside
        the matching CHEMBL ID's


    Returns:
    --------
    A dictionary of lists, keys are query molecules.

    e.g if `show_similarity is False`:
        {search_chembl_id_0: [similar1, similar2, ...],
         search_chembl_id_1: [similar2, similar3, ...]}

    e.g if `show_similarity is True`:
        {search_chembl_id_0: [(similar1, 90.2), (similar2, 94.0), ...],
         search_chembl_id_1: [(similar2, 98.4), (similar3, 92.1), ...]}
    """
    similar_molecules = {}
    similarity_query = new_client.similarity
    for compound in chembl_ids:
        similars = []
        res = similarity_query.filter(chembl_id=compound, similarity=similarity)
        for entry in res:
            if entry["molecule_chembl_id"] == compound:
                # ignore query compound similarity matching itself
                continue
            if show_similarity:
                similars.append((entry["molecule_chembl_id"],
                                 float(entry["similarity"])))
            else:
                similars.append(entry["molecule_chembl_id"])
        similar_molecules[compound] = similars
    return similar_molecules


def get_similar_molecules_smile(smiles, similarity=90, show_similarity=False):
    """
    Get similar moleules based on a structural similarity.
    Given a list of structures in SMILE format.

    Note:
    -----
    Similarity matches may appear to be ordered by descening similarity score,
    though this is not guaranteed and should not be relied upon.


    Parameters:
    -----------
    smiles: list-like
        SMILE strings

    similarity: numeric, between 0 and 100 (default = 90)
        minimum structural similarity threshold to determine a match

    show_similarity: Boolean (default = False)
        whether or not to return the similarity score alongside
        the matching CHEMBL ID's


    Returns:
    --------
    A dictionary of lists, keys are query molecules.

    e.g if `show_similarity is False`:
        {smile1: [similar1, similar2, ...],
         smile2: [similar2, similar3, ...]}

    e.g if `show_similarity is True`:
        {smile1: [(similar1, 90.2), (similar2, 94.0), ...],
         smile2: [(similar2, 98.4), (similar3, 92.1), ...]}
    """
    similar_molecules = {}
    similarity_query = new_client.similarity
    for compound in smiles:
        similars = []
        res = similarity_query.filter(smiles=compound, similarity=similarity)
        for entry in res:
            # TODO: check if 'similar' is actually an exact match, if so should
            # this be returned? Can determine this by similarity == 100 or
            # identical canonical SMILE strings.
            if show_similarity:
                similars.append((entry["molecule_chembl_id"],
                                 float(entry["similarity"])))
            else:
                similars.append(entry["molecule_chembl_id"])
        similar_molecules[compound] = similars
    return similar_molecules


def get_target_ids(chembl_ids, organism="Homo sapiens", ignore_empty=False,
                   standard_value_threshold=None):
    """
    TODO: docstring

    Parameters:
    -----------
    chembl_ids: list-like

    organism: string (default = "Homo sapiens")

    ignore_empty: Boolean (default = False)

    standard_value_threshold: numeric (default = None/infinity)



    Returns:
    ---------
    Dictionary
    """
    compound_target_dict = _get_target_ids_as_chembl(
            chembl_ids, organism, standard_value_threshold
    )
    chunk_size = 50

    for key, val in compound_target_dict.items():
        lval = list(val)
        uniprots = set()
        for i in range(0, len(val), chunk_size):
            targets = new_client.target\
                .filter(target_chembl_id__in=lval[i:i + chunk_size])
            uniprots_tmp = set()
            for target in targets:
                if target["target_type"] == "SINGLE PROTEIN":
                    for component in target["target_components"]:
                        uniprots_tmp.add(component["accession"])
                uniprots_tmp = set(uniprots_tmp)
                uniprots = uniprots.union(uniprots_tmp)
            compound_target_dict[key] = uniprots
    # FIXME: pretty inefficient. shouldn't be creating the empty sets in the
    #        first place
    if ignore_empty:
        compound_target_dict = {key: value for key, value in \
            compound_target_dict.items() if len(value) > 0}
    return compound_target_dict


def _get_target_ids_as_chembl(chembl_ids, organism="Homo sapiens",
                              standard_value_threshold=None):
    """
    internal function: Returns target ID, as CHEMBL_ID's
    used within get_target_ids

    Returns:
    --------
    dictionary
    """
    if isinstance(chembl_ids, str):
        chembl_ids = [chembl_ids]
    # set default  threshold to infinity, so all values are below it,
    # therefore all targets are returned
    if standard_value_threshold is None:
        standard_value_threshold = float("Inf")

    chembl_ids = list(chembl_ids)
    chunk_size = 50
    compounds2targets = {chembl: set() for chembl in chembl_ids}

    ID_forms = dict()
    for x in chembl_ids:
        ID_forms[x] = set()

    for i in range(0, len(chembl_ids), chunk_size):
        for form in new_client.molecule_form\
            .filter(parent_chembl_id__in=chembl_ids[i:i + chunk_size]):
            ID_forms[form['parent_chembl_id']]\
                .add(form['molecule_chembl_id'])

    for i in range(0, len(chembl_ids), chunk_size):
        for form in new_client.molecule_form\
            .filter(molecule_chembl_id__in=chembl_ids[i:i + chunk_size]):
            ID_forms[form['molecule_chembl_id']]\
                .add(form['parent_chembl_id'])

    values = []
    for x in ID_forms.values():
        values.extend(x)
    forms_to_ID = dict()
    for x in values:
        forms_to_ID[x] = set()

    for k in forms_to_ID:
        for parent, molecule in ID_forms.items():
            if k in molecule:
                forms_to_ID[k] = parent

    for i in range(0, len(values), chunk_size):
        activities = new_client.activity\
            .filter(molecule_chembl_id__in=values[i:i + chunk_size])\
            .filter(target_organism__istartswith=organism)
            # FIXME: could potentially filter the standard value here
        for act in activities:
            # if assay concentration readout is less than the threshold
            # then add as a target
            try:
                # sometimes returns a NoneType if there is no standard value
                # so place in the try:except block
                standard_val = float(act["standard_value"])
                if standard_val < standard_value_threshold and act["standard_units"] == "nM":
                    compounds2targets[forms_to_ID[act['molecule_chembl_id']]]\
                        .add(act['target_chembl_id'])
            except TypeError:
                # if we get a type error then there isn't a standard value
                # associated with this target, so skip it
                continue


    return compounds2targets
