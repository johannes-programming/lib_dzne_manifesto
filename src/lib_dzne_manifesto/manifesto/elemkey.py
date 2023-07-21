"""This module contains some values and function """\
"""that are used by multiple modules in bindBASE.manifesto. """


import sys

import pandas as pd


# ToDo: replace with the gene order from igblast
def _part2float(s):
    """Converting the part of of a gene/family-name into float. """
    assert len(s) < 5
    if s == 'N/A':
        return -1
    if s == 'NL1':
        return 1.9
    if s.endswith('D'):
        s = s[:-1] + ".001"
    return float(s)
def main(s):
    """Converting a gene/family-name into float. """\
    """This is done so that the ticklabels can be in a sensible order. """
    if s == "":# should this even be allowed to happen?
        return -1.0e24
    ans = 0.0
    for i, x in enumerate(s.split('-')):
        ans += _part2float(x) * (1.0e12 ** (-i))
    return ans









