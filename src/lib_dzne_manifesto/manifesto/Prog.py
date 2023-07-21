"""This module contains the Prog-class. """



import os
import string
import sys
from collections import defaultdict

import lib_dzne_basetables.table as tbl
import lib_dzne_math.statistics
import lib_dzne_math.statistics.plots
import pandas as pd
from BASE_dzne import config
from BASE_dzne.binding import levenshtein

from . import chain, elemkey, masterfile, slim


class Prog:
    """This class creates the manifesto. """
    def __init__(self, *, table):
        """Creating all the variables. """

        # a manifesto can only be created from an mFile.
        # In truth this is not the only condition, but it is easily checked for. 
        # Therefore it is tested here.
        tbl.check(table=table, filetype='m')
        self.table = table
        config = config.get_config()
        self.cfg = config['m']
        self.config = config
        self.stats_rz = list()

        # basic preparing of mTable
        ### we are retaining an unaltered copy of the original table
        ### our working duplicate is named as a DataFrame 
        ### as we are about to diviate from the standards we have set for tables
        self.dfB = tbl.make(table)
        ### we are adding the IDs
        rz = list()
        for i, row in self.dfB.iterrows():
            row['antibody id'] = self.cfg['IDsep'].join(row.get(col, "") for col in ['COLLECTION', 'NR', 'LIGHT_CHAIN_TYPE'])
            row['original cell id'] = self.cfg['IDsep'].join(row.get(col, "") for col in ['COLLECTION', 'NR'])
            rz.append(row)
        self.dfB = pd.DataFrame(rz)
        ### We demand that a group-column exists:
        ### the GROUP column is renamed
        ### so that we will not need to alter the legends
        ### of the plots after their creation
        if 'GROUP' in self.dfB.columns:
            self.dfB = self.dfB.rename(columns={'GROUP': 'group'})
            self.dfB = self.dfB.sort_values(by=['group'])
        else:
            self.dfB['group'] = [""] * table.shape[0]


        # make sure that each hue always corresponds to the same group
        self.hue_order = list(set(self.dfB['group']))
        self.hue_order.sort()
        assert 'all' not in self.hue_order

        # dfA
        self.dfA = self.dfB.copy()
        self.dfA['group'] = ['all'] * self.dfA.shape[0]
        self.dfA = pd.concat(
            [self.dfB, self.dfA], 
            axis=0, 
            ignore_index=True,
        )

        # create chain_df
        # (the conversion differs from the m2d-conversion
        # because one heavy chain ends up as two rows
        # if it has both a kappa and a lambda chain)
        self.chain_dfA = chain.main(self.dfA)
        self.chain_dfA = self.chain_dfA.rename(columns={'REGION_CDR3_TR_LEN': 'cdr3 length'})
        self.chain_dfB = self.chain_dfA['all' != self.chain_dfA['group']]

        #
        self.figsize = self.cfg.get('figsize', "")











    
    def run(self, /):
        """This generator ultimately creates the manifesto. """
        """It is the endgame of BASE. """
        """Every yield is a Page-object which will ultimately be converted into one output-file. """
        for s in ['x', 'wb', 'slim', 'AB', 'chainAB', 'cdr3l', 'cdr3l_scatter', 'count', 'vj', 'levenshtein', 'stats']:
            f = getattr(self, f"run_{s}")
            for x in f():
                yield x




    def run_x(self, /):
        """Yielding the original mFile for the sake of completeness. """
        yield DataEntry.make(
            filename="x.mbase",
            data=self.table,
        )




    def run_wb(self, /):
        """Yielding the masterfile. """
        wb = masterfile.main(table=self.table)
        yield DataEntry.make(
            filename='masterfile.xlsx',
            data=wb,
        )


    def run_slim(self, /):
        """With the help of some config-values """\
        """a human-readable shortversion of the input-file is created. """
        try:
            c = self.cfg['slim']
        except KeyError as exc:
            return
        slim_table = slim.main(
            table=self.table,
            cfg=c,
        )
        yield DataEntry.make(
            filename='slime.mbase', 
            data=slim_table,
        )



    def run_AB(self, /):
        """Yielding the chain-tables that are used """\
        """by most of the other run-generators. """
        yield DataEntry.make(
            data=self.dfA,
            filename='dfA.csv',
        )
        yield DataEntry.make(
            data=self.dfB,
            filename='dfB.csv',
        )


    def run_chainAB(self, /):
        yield DataEntry.make(
            data=self.chain_dfA,
            filename="chain_dfA.csv",
        )
        yield DataEntry.make(
            data=self.chain_dfB,
            filename="chain_dfB.csv",
        )


    def run_cdr3l(self, /):
        # cdr3l
        cdr3l_dfA = self.chain_dfA.copy()
        cdr3l_dfA = cdr3l_dfA[~cdr3l_dfA['cdr3 length'].isin(["", 'N/A'])]
        cdr3l_dfA['cdr3 length'] = cdr3l_dfA['cdr3 length'].map(int)
        cdr3l_dfB = cdr3l_dfA['all' != cdr3l_dfA['group']]
        cdr3l_catplot_infos = (
            ('h', 'cdr3 length', 'chain type'),
            ('v', 'chain type', 'cdr3 length')
        )
        addhuess = {
            'A': ['all'],
            'B': []
        }
        for switch, x, y in cdr3l_catplot_infos:
            for kind in ('bar', 'box'):
                for letter, addhues in addhuess.items():
                    yield DataEntry(
                        data=dict(
                            plot=sns.catplot,
                            figsize=self.figsize,
                            title="",
                            data=cdr3l_dfA,
                            x=x,
                            y=y,
                            hue='group',
                            hue_order=self.hue_order+addhues,
                            order=['HC', 'KC', 'LC'],
                            kind=kind,
                        ),
                        filename=f"CDR3_TR_LEN_{kind}_{switch}_{letter}.png",
                    )


    def run_cdr3l_scatter(self, /):
        """Creating a scatterplot of cdr3 lengths and cdr3 gravy scores """\
        """(heavy chain versus light chain). """
        
        AB_scatterplot_infos = (
            ('REGION_CDR3_TR_LEN', 'cdr3 length', int),
            ('REGION_CDR3_GRAVY', 'cdr3 gravy score', float)
        )
        for col, humancol, datatype in AB_scatterplot_infos:
            scatterplot_df = self.dfB.copy()
            for cw in ('HEAVY', 'LIGHT'):
                title = f'{cw.lower()} chain'
                scatterplot_df = scatterplot_df.rename(columns={
                    f'{cw}_BEST_{col}': title
                })
                scatterplot_df = scatterplot_df[~scatterplot_df[title].isin(["", 'N/A'])]
                scatterplot_df[title] = scatterplot_df[title].map(datatype)
            scatterplot_df = scatterplot_df[['heavy chain', 'light chain', 'group', 'antibody id']]
            yield DataEntry.make(
                data=scatterplot_df, 
                filename=f"{col}_heavy_light_scatterplot.csv",
            )
            yield DataEntry.make(
                data=dict(
                    figsize=self.figsize,
                    plot=sns.scatterplot,
                    data=scatterplot_df,
                    x="heavy chain",
                    y="light chain",
                    hue='group',
                ),
                filename=f"{col}_heavy_light_scatterplot.png",
            )




    def run_count(self, /):
        #return
        """Count the items in certain columns and display their counts in a barplot. """
        addhuess = {'A': ['all'], 'B': []}
        countchars = {'n': 'n', '%': 'p'}
        for CT in ('HC', 'KC', 'LC'):
            df_CT = self.chain_dfA.loc[self.chain_dfA['chain type'] == CT]
            for col, axislabelf in chain.gene_col_axis_labels.items():
                axislabel = axislabelf.format(ct=CT[0])
                l = ['group', col]
                rz = list()
                totals = df_CT['group'].value_counts()
                for k, n in df_CT[l].value_counts().to_dict().items():
                    rl = dict(zip(l, k))
                    rl['n'] = n
                    total = totals[rl['group']]
                    rl['total'] = total
                    rl['%'] = 100 * n / total
                    rz.append(rl)
                df_q = pd.DataFrame(rz)
                df_q = df_q.rename(columns={col: axislabel})
                order = df_q[axislabel].tolist()
                order = list(set(order))
                order.sort(key=elemkey.main)
                yield DataEntry.make(
                    data=df_q, 
                    filename=f"counter_{CT}_{col}.csv",
                )
                for counttype in "n%":
                    for huechar in "AB":
                        for switch in "vh":# vertical or horizontal bars
                            yield DataEntry.make(
                                filename="counter_{CT}_{col}_{h}_{s}_{c}.png".format(
                                    CT=CT,
                                    col=col,
                                    h=huechar,
                                    s=switch,
                                    c=countchars[counttype]
                                ),
                                data=dict(
                                    plot=sns.catplot,
                                    figsize=self.figsize,
                                    title="",
                                    kind='bar',
                                    x=axislabel if switch == 'v' else counttype,
                                    y=counttype if switch == 'v' else axislabel,
                                    order=order,
                                    hue='group',
                                    hue_order=self.hue_order+addhuess[huechar],
                                    data=df_q,
                                    ci=None,
                                ),
                            )



    def run_vj(self, /):
        """Plotting the vj-pairings into heatmaps. """
        for ct in "HKL":
            CT = ct + 'C'
            ct_df = self.chain_dfA[CT == self.chain_dfA['chain type']]
            v_col = 'TOP_V_FAM'
            j_col = 'TOP_J_GENE'
            try:
                v_labels = ct_df[v_col]
                j_labels = ct_df[j_col]
            except:
                continue
            v_labels = list(set(v_labels.tolist()))
            j_labels = list(set(j_labels.tolist()))
            v_labels.sort(key=elemkey.main)
            j_labels.sort(key=elemkey.main)
            for group_n, group in enumerate(self.hue_order):
                groupcode = ""
                for ch in group:
                    if ch in (string.ascii_letters + string.digits):
                        groupcode += ch
                    else:
                        groupcode += '_'
                groupcode += '_' + str(group_n)
                group_ct_df = ct_df[group == ct_df['group']]
                rz = {v: {j: 0 for j in j_labels} for v in v_labels}
                for i, row in group_ct_df.iterrows():
                    rz[row[v_col]][row[j_col]] += 1
                vj_df = pd.DataFrame(rz)
                filename = f"vj-heatmap_{groupcode}_{CT}.csv"
                yield DataEntry.make(
                    data=vj_df, 
                    filename=filename,
                )
                heat_dfs = {"vj": vj_df.copy(), "jv": vj_df.transpose()}
                for k, heat_df in heat_dfs.items():
                    yield DataEntry.make(
                        filename=f"{k}-heatmap_{groupcode}_{CT}.png",
                        data=dict(
                            figsize=self.figsize,
                            plot=sns.heatmap,
                            data=heat_df,
                            square=True,
                            annot=True,
                            fmt="d",
                            cbar=False,
                            xlabel=f"{k[0].upper()}{ct}",
                            ylabel=f"{k[1].upper()}{ct}",
                            title=f"{k}-pairing for {group}",
                        ),
                    )


    def run_levenshtein(self, /):
        """Utilizing the levenshtein program. """
        for subseq_name in ('VAR', 'REGION_CDR3'):
            for chain_type in ('Kappa', 'Lambda'):
                CT = chain_type.upper()[0] + 'C'
                table = self.dfB[CT == self.dfB['LIGHT_CHAIN_TYPE']]
                IDs = table['original cell id'].tolist()
                difference_matrices = dict()
                for cw in 'HEAVY', 'LIGHT':
                    col = f'{cw}_BEST_{subseq_name}_TR'  # Jakob wanted TR instead of SEQ
                    if col not in table.columns:
                        continue
                    seqs = table[col].tolist()
                    seqs = [(x if (x != 'N/A') else "") for x in seqs]
                    cw_difference_matrix = levenshtein.get_matrix(args=seqs)
                    df = pd.DataFrame(cw_difference_matrix.tolist(), columns=IDs, index=IDs)
                    yield DataEntry.make(
                        data=df,
                        filename=f"{chain_type.lower()}-ABs_{cw.lower()}_{subseq_name}.csv",
                    )
                    difference_matrices[cw] = cw_difference_matrix
                difference_matrix = difference_matrices['HEAVY'] + difference_matrices['LIGHT']
                df = pd.DataFrame(difference_matrix.tolist(), columns=IDs, index=IDs)
                yield DataEntry(
                    data=df,
                    filename=f"{chain_type.lower()}-ABs_{subseq_name}.csv",
                )
                pca_difference_matrix = lib_dzne_math.statistics.pca(difference_matrix, 2)
                df = pd.DataFrame(data=pca_difference_matrix.tolist(), columns=["x", "y"], index=IDs)
                df['antibody id'] = IDs
                df['group'] = table['group'].tolist()
                yield DataEntry(
                    data=df,
                    filename=f"{chain_type.lower()}-ABs_{subseq_name}_pca.csv",
                )
                yield DataEntry(
                    filename=f"{chain_type.lower()}-ABs_{subseq_name}_pca.png",
                    data=dict(
                        figsize=self.figsize,
                        plot=sns.scatterplot,
                        data=df,
                        x="x",
                        y="y",
                        hue='group',
                        hue_order=self.hue_order,
                        title=f"Antibodies with {chain_type.lower()}-chains",
                        xlabel="",
                        ylabel="",
                    ),
                )


    def run_stats(self, /):
        stats_table = pd.DataFrame(self.stats_rz)
        check(stats_table)
        yield DataEntry.make(
            data=stats_table,
            filename='stats_table.csv',
        )








