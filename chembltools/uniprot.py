"""
Functions for returning UNIPROT data, typically GO-terms.
"""


import urllib.request
import re


def get_uniprot_data(uniprot_code, fasta=False, decode=False):
    """
    internal function to fetch all uniprot text data for a single uniprot ID
    """
    url = "http://www.uniprot.org/uniprot/{}.txt".format(uniprot_code)
    if fasta:
        url = url.replace("txt", "fasta")
    if decode is False:
        return urllib.request.urlopen(url)
    else:
        # decode into string
        data = urllib.request.urlopen(url)
        return [line.decode("utf-8") for line in data]


def warn_missing_uniprot(identifier, err_code):
    """
    internal function to warn if identifier is not found in uniprot
    """
    msg  = "Warning: Could not find '{}' entry on uniprot".format(identifier)
    msg2 = "-- error code = {}.".format(err_code)
    print(msg, msg2)


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
                warn_missing_uniprot(identifier, err.code)
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
                warn_missing_uniprot(identifier, err.code)
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
                warn_missing_uniprot(identifier, err.code)
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
                warn_missing_uniprot(identifier, err.code)
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
                warn_missing_uniprot(identifier, err.code)
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


def accession_to_gene_name(codes, cache=True):
    """
    Convert uniprot accession codes to gene-names
    """
    if isinstance(codes, str):
        codes = [codes]
    if cache:
        result = []
        id_cache = {}
        for uniprot_id in codes:
            try:
                gene_name = id_cache[uniprot_id]
            except KeyError:
                try:
                    uniprot_data = get_uniprot_data(
                        uniprot_id, fasta=True, decode=True
                    )
                    gene_name = _get_gene_name(uniprot_data[0])
                    id_cache[uniprot_id] = gene_name
                except urllib.error.HTTPError as err:
                    if err.code == 404 or err.code == 300:
                        warn_missing_uniprot(uniprot_id, err.code)
                        continue
                    else:
                        raise err
            result.append(gene_name)
        return result
    elif cache is False:
        fasta_data = [get_uniprot_data(i, fasta=True, decode=True)[0] for i in codes]
        return [_get_gene_name(i) for i in fasta_data]
    else:
        raise TypeError("cache requires a Boolean")


def _get_gene_name(messy_string):
    """
    Get string between "GN=" and first ";" or " ".
    Used for getting the gene-name out of a fasta file.
    """
    result = re.search("(?<=GN=)(.*?)(?=;| )", messy_string)
    if out is None:
        raise RuntimeError("Could not parse gene name from FASTA string")
    return result.group(1)

