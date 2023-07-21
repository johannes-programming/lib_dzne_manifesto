"""This module concerns itself with the creating of the masterfile """\
"""- an Excel-file that is populated with the data from the mTable. """




import sys

from ... import Config


def conv(x):
    """Convert value to best possible datatype. """
    for t in (int, float):
        try:
            return t(x)
        except ValueError as exc:
            continue
    return x

def set(*, cell, value):
    """Setting value of cell. """
    cell.value = value
    #cell.alignment = Alignment()#horizontal='general')

def main(*, table):
    """Writing data from table into masterfile-template. """
    masterrow = Config.get_config()['m']['masterrow']
    if masterrow == "":
        masterrow = 1
    else:
        masterrow = int(masterrow)
    wb = Config.get_workbook()
    ws = wb.active
    for colnum in range(1, ws.max_column + 1):
        c = ws.cell(column=colnum, row=masterrow)
        v = c.value
        if type(v) is not str:
            continue
        if v.startswith('='):
            continue
        if v == "":
            continue
        if v not in table.columns:
            set(cell=c, value="")
            continue
        for i, x in enumerate(table[v].tolist()):
            y = conv(x)
            _c = ws.cell(column=colnum, row=masterrow+i)
            set(cell=_c, value=y)
    return wb






