"""
module docstring
"""


import urllib.request


def get_uniprot_data(uniprot_code):
    """
    internal function to fetch all uniprot text data for a single uniprot ID
    """
    url = "http://www.uniprot.org/uniprot/{}.txt".format(uniprot_code)
    return urllib.request.urlopen(url)


def warn_missing_uniprot(identifier):
    """
    internal function to warn if identifier is not found in uniprot
    """
    msg = "Warning: Could not find {} entry on uniprot".format(identifier)
    print(msg)


def get_uniprot_name(uniprot_ids):
    """
    Given a list of uniprot IDs / accession codes, return a dictionary
    mapping {uniprot_id: protein_name}

    Parameters:
    -----------
    uniprot_ids: list-like
        uniprot_id / accession codes


    Returns:
    --------
    Dictionary {uniprot_id: protein_name}
    """
    uniprot_name_dict = {}
    for identifier in uniprot_ids:
        try:
            data = get_uniprot_data(identifier)
        except urllib.error.HTTPError as err:
            if err.code == 404 or err.code == 300:
                warn_missing_uniprot(identifier)
                continue
            else:
                raise err
        for line in data:
            line = line.decode("utf-8")
            if line.startswith("DE   RecName:"):
                name = line.split("Full=")[-1].strip(";\n")
                break # already have the name, no need to keep iterating
        uniprot_name_dict[identifier] = name
    return uniprot_name_dict


def get_uniprot_info(uniprot_ids):
    """
    get the entire uniprot website text for a list of uniprot ID's
    """
    uniprot_info_dict = {}
    for identifier in uniprot_ids:
        info = []
        try:
            data = get_uniprot_data(identifier)
        except urllib.error.HTTPError as err:
            if err.code == 404 or err.code == 300:
                warn_missing_uniprot(identifier)
                continue
            else:
                raise err
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
        try:
            data = get_uniprot_data(identifier)
        except urllib.error.HTTPError as err:
            if err.code == 404 or err.code == 300:
                warn_missing_uniprot(identifier)
                continue
            else:
                raise err
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
        try:
            data = get_uniprot_data(identifier)
        except urllib.error.HTTPError as err:
            if err.code == 404 or err.code == 300:
                warn_missing_uniprot(identifier)
                continue
            else:
                raise err
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
        try:
            data = get_uniprot_data(identifier)
        except urllib.error.HTTPError as err:
            # handle 404 errors by skipping that identifier
            if err.code == 404 or err.code == 300:
                warn_missing_uniprot(identifier)
                continue
            else:
                raise err
        for line in data:
            line = line.decode("utf-8")
            if line.startswith("DR   GO"):
                go_codes.append(line.split(";")[1].strip())
                go_names.append(line.split(":")[-2].split(";")[-2])
        go_dict[identifier] = list(zip(go_codes, go_names))
    return go_dict

