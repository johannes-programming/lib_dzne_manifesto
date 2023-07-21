

import BASE_dzne.libBASE.table as tbl

from ... import Config


def main(*, table, cfg):
    """Making a simplified mTable. """
    drop = cfg.get('drop')
    if drop is not None:
        drop = tbl.identify_columns(table=table, patterns=drop)
        table = table.drop(columns=drop)
    keep = cfg.get('keep')
    if keep is not None:
        keep = tbl.identify_columns(table=table, patterns=keep)
        table = table[keep]
    return table





