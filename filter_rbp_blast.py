import argparse
import pandas as pd


def filtre_rbp(rep, out):
    pd_bl = pd.read_csv(rep, sep='\t', index_col="assembly")
    pd_bl_grouped = pd_bl.groupby("assembly")

    filtered_rbps = pd.DataFrame()
    for ass, rows in pd_bl_grouped:
        print(f"{ass}")
        # df = rows.sort_values(['sstart', 'send'], ascending=[True, False])
        # df['group'] = ((df['sstart'].cummax().shift() <= df['sstart'] and df['send'].cummax().shift() <= df['send'])).cumsum()
        # df.sort_index(inplace=True)

        ass_sseqid_preds = rows.groupby(['sseqid'])
        for sseqid, group in ass_sseqid_preds:
            if len(group.index) == 1:
                filtered_rbps = pd.concat([filtered_rbps, group])
            else:
                # same seqid with mutiple coordinates
                # add key column
                group['combined'] = group['sstart'].astype(str) + '_' + group['send'].astype(str)

                # need to merge overlapping coords
                def merge_overlapping(row):
                    overlapped_by = find_overlapping(row, group)
                    return overlapped_by

                group = group.sort_values(['pident', 'qcovs'], ascending=[False, False])
                group['pred_overlap'] = group.apply(merge_overlapping, axis=1)
                group = group.loc[group['pred_overlap'].isnull()]
                # if multiple none found then take the first pred of group sorted by qcov and ident
                group = group.sort_values(['pident', 'qcovs'], ascending=[False, False])
                coord_group = group.groupby(['combined'])
                for sseqid, group in coord_group:
                    if len(group.index) == 1:
                        group = group.drop(['pred_overlap', 'combined'], axis=1)

                    else:
                        group = pd.DataFrame(group.head(1))
                        group = group.drop(['pred_overlap', 'combined'], axis=1)

                    filtered_rbps = pd.concat([filtered_rbps, group])
    filtered_rbps.to_csv(
        out,
        sep="\t")


def find_overlapping(row, group):
    for index, no in group.iterrows():
        if row['qacc'] != no['qacc']:
            x = check_overlap(row, no)
            if x is not None:
                return x

    return None


def check_overlap(pred, no):
    # if pred['sacc'] != no['sacc']:
    #     return False
    # included in pred coord pred
    if pred['sstart'] >= no['sstart'] and pred['send'] <= no['send']:
        if check_qcov_ident(pred, no):
            return no['combined']
    # case prior start to no start and pred end after no start
    elif pred['sstart'] <= no['sstart'] <= pred['send'] and pred['send'] <= no['send']:
        if check_qcov_ident(pred, no):
            return no['combined']
    # case pred start after no start and pred end after no end
    elif no['sstart'] <= pred['sstart'] <= no['send'] and pred['send'] >= no['send']:
        if check_qcov_ident(pred, no):
            return no['combined']
    # case complete overlap
    elif pred['sstart'] <= no['sstart'] and pred['send'] >= no['send']:
        if check_qcov_ident(pred, no):
            return no['combined']
    return None


# check qcov and ident within same range or not
def check_qcov_ident(pred, no):
    if pred['pident'] <= (no['pident'] - 5):
        return True

    # if pred['qcovs'] <= (no['qcovs'] - 5):
    #     return True

    return False


if __name__ == '__main__':
    argParser = argparse.ArgumentParser()

    # mandatory
    argParser.add_argument(
        "-r",
        "--blast_report",
        help="blast output report with assembly label",
        required=True
    )

    argParser.add_argument(
        "-o",
        "--out",
        help="output filepath",
        required=True
    )

    args = argParser.parse_args()

    filtre_rbp(
        args.blast_report,
        args.out
    )
