"""
docstring
"""

import urllib.request
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


def get_similar_molecules(chembl_ids, similarity=90):
    """
    Get similar moleules based on a structural similarity search.
    Given a list of chembl ID's

    Returns:
    --------
    A dictionary of lists, keys are query molecules.
    e.g:
        {search_chembl_id_0: [similar1, similar2, ...],
         search_chembl_id_1: [similar2, similar3, ...]}
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
            similars.append(entry["molecule_chembl_id"])
        similar_molecules[compound] = similars
    return similar_molecules


def get_similar_molecules_smile(smiles, similarity=90):
    """
    Get similar molecules based on structural similarity search.
    Given a list of structures in SMILE format.

    Returns:
    --------
    A dictionary of lists, keys are query molecules.
    e.g:
        {search_smile_0: [similar_chembl_ids ...],
         search_smile_1: [similar_chembl_ids ...]}
    """
    similar_molecules = []
    similarity_query = new_client.similarity
    for compound in smiles:
        similars = []
        res = similarity_query.filter(smiles=compound, similarity=similarity)
        for entry in res:
            similars.append(entry["molecule_chembl_id"])
        similar_molecules[compound] = similars
    return similar_molecules


def get_target_ids(chembl_ids, organism="Homo sapiens"):
    """
    docstring
    """
    chembl_ids = list(chembl_ids)
    chunk_size = 50
    compounds2targets = {chembl: set() for chembl in chembl_ids}

    ID_forms = dict()
    for x in chembl_ids:
        ID_forms[x] = set()

    for i in range(0, len(chembl_ids), chunk_size):
        for form in new_client.molecule_form.filter(parent_chembl_id__in=chembl_ids[i:i + chunk_size]):
            ID_forms[form['parent_chembl_id']].add(form['molecule_chembl_id'])

    for i in range(0, len(chembl_ids), chunk_size):
        for form in new_client.molecule_form.filter(molecule_chembl_id__in=chembl_ids[i:i + chunk_size]):
            ID_forms[form['molecule_chembl_id']].add(form['parent_chembl_id'])

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
        for act in activities:
            compounds2targets[forms_to_ID[act['molecule_chembl_id']]]\
                .add(act['target_chembl_id'])

    for key, val in compounds2targets.items():
        lval = list(val)
        uniprots = set()
        for i in range(0, len(val), chunk_size):
            targets = new_client.target.filter(target_chembl_id__in=lval[i:i + chunk_size])
            uniprots = uniprots.union(
                set(sum([[comp["accession"] for comp in t["target_components"]] for t in targets],[]))
            )
        compounds2targets[key] = uniprots
    return compounds2targets


def get_uniprot_name(uniprot_ids):
    """
    docstring
    """
    uniprot_name_dict = {}
    for identifier in uniprot_ids:
        data = _get_go_data(identifier)
        for line in data:
            line = line.decode("utf-8")
            if line.startswith("DE   RecName:"):
                name = line.split("Full=")[-1]
                break
        uniprot_name_dict[identifier] = name
    return uniprot_name_dict


def get_uniprot_info(uniprot_ids):
    """
    get the entire uniprot website text for a list of uniprot ID's
    """
    uniprot_info_dict = {}
    for identifier in uniprot_ids:
        info = []
        data = _get_go_data(identifier)
        for line in data:
            info.append(line.decode("utf-8"))
        uniprot_info_dict[identifier] = info
    return uniprot_info_dict


def get_go_name_from_uniprot_id(uniprot_ids):
    """
    get descriptive go terms
    """
    go_term_dict = {}
    for identifier in uniprot_ids:
        go_terms = []
        data = _get_go_data(identifier)
        for line in data:
            line = line.decode("utf-8")
            if line.startswith("DR   GO"):
                go_term = line.split(":")[-2].split(";")[-2]
                go_terms.append(go_term)
        go_term_dict[identifier] = go_terms
    return go_term_dict


def get_go_code_from_uniprot_id(uniprot_ids):
    """
    get go codes e.g "GO:003674"
    """
    go_code_dict = {}
    for identifier in uniprot_ids:
        go_codes = []
        data = _get_go_data(identifier)
        for line in data:
            line = line.decode("utf-8")
            if line.startswith("DR   GO"):
                go_code = line.split(";")[1].strip()
                go_codes.append(go_code)
        go_code_dict[identifier] = go_codes
    return go_code_dict


def get_go_from_uniprot_id(uniprot_ids):
    """
    get both go code and descriptive name for a list of uniprot ID's
    e.g.
        {uniprot_id: [go_code, go_name]}
    """
    go_dict = {}
    for identifier in uniprot_ids:
        go_codes, go_names = [], []
        data = _get_go_data(identifier)
        for line in data:
            line = line.decode("utf-8")
            if line.startswith("DR   GO"):
                go_codes.append(line.split(";")[1].strip())
                go_names.append(line.split(":")[-2].split(";")[-2])
        go_dict[identifier] = list(zip(go_codes, go_names))
    return go_dict


def _get_go_data(uniprot_code):
    """
    internal function to fetch all uniprot text data for a single uniprot ID
    """
    url = "http://www.uniprot.org/uniprot/{}.txt".format(uniprot_code)
    data = urllib.request.urlopen(url)
    return data


