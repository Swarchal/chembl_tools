"""
module docstring
"""


def common(list_of_lists):
    """
    Find intersection of elements between list of lists.
    Useful for finding common GO-terms.

    Returns:
    ---------
    set of GO terms
    """
    return set(list_of_lists[0]).intersection(*list_of_lists[1:])

