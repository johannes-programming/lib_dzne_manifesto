"""This module contains the function that produces a chain-based df. """


import pandas as pd

gene_col_axis_labels = {
    'TOP_V_GENE': "V{ct} gene",
    'TOP_V_FAM': "V{ct} gene family",
    'TOP_J_GENE': "J{ct} gene"
}




def main(df):
    """Deconstructing the table so that one row corresponds to one chain. """
    non_chain_cols = ['group', 'antibody id', 'original cell id']
    chain_cols = ['REGION_CDR3_TR_LEN'] + list(gene_col_axis_labels.keys())
    rz = list()
    for i, row in df.iterrows():
        for cw in ('HEAVY', 'LIGHT'):
            rl = dict()
            if cw == 'HEAVY':
                rl['chain type'] = 'HC'
            else:
                rl['chain type'] = row['LIGHT_CHAIN_TYPE']
            for col in non_chain_cols:
                rl[col] = row[col]
            for col in chain_cols:
                rl[col] = row.get(f"{cw}_BEST_{col}", "")
            rz.append(rl)
    return pd.DataFrame(rz)


