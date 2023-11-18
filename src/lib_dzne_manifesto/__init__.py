#import sys

#import lib_dzne_sqlite.exec
#import BASE_dzne.heart.comparing
#import BASE_dzne.heart.database.lookup.primary_key
#import BASE_dzne.heart.database.lookup.best_id
#import BASE_dzne.heart.database.lookup.parent_id
#from BASE_dzne.heart.database.Connection import Connection
#from BASE_dzne.heart.analyzing.Seqread import Seqread

import na_quantors as _na_quantors


def get_antibody_data(antibody_id):
    with Connection() as connection:
        cursor = connection.cursor()
        df = lib_dzne_sqlite.exec.select(
            cursor=cursor,
            table='antibodies',
            columns=['*'],
            where={'id':antibody_id},
        )
    antibodyrow, = [row.to_dict() for i, row in df.iterrows()]
    ans = dict()
    ans['antibody'] = dict()
    for k, v in antibodyrow.items():
        newkey = k.replace('_', '-')
        if newkey in ans['antibody'].keys():
            raise KeyError()
        ans['antibody'][newkey] = v
    for cw in ('heavy', 'light'):
        ans[cw] = get_chain_data(
            collection_id=antibodyrow['collection_id'],
            nr=antibodyrow['nr'],
            chain_type=antibodyrow[f'{cw}_chain_type'],
        )
    return ans

def get_chain_data(*, collection_id, nr, chain_type):
    chain_key = dict(
        collection_id=collection_id,
        nr=nr,
        chain_type=chain_type,
    )
    ans = dict()
    best_id = BASE.heart.database.lookup.best_id.by_chain_key(**chain_key)
    if _na_quantors.isna(best_id):
        return ans
    best = Seqread(best_id)
    best.calc()
    ans['seqread'] = best.data
    parent_id = BASE.heart.database.lookup.parent_id.by_child_id(best_id)
    if _na_quantors.isna(parent_id):
        return ans
    parent = Seqread(parent_id)
    parent.calc()
    ans['parent'] = parent.data
    comparison = BASE.heart.comparing.main(
        parent=parent,
        child=best,
    )
    #comparison.pop('primers')
    #comparison.pop('regions')
    ans['comparison'] = comparison
    print(ans, file=sys.stderr)
    return ans






 
