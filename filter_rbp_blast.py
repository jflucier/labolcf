import argparse
import pandas as pd


def filtre_rbp(rep):
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

                # aggregation_functions = {'price': 'sum', 'amount': 'sum', 'name': 'first'}
                # df_new = df.groupby(df['id']).aggregate(aggregation_functions)
                # merged = pd.merge(group,
                #                   group,
                #                   left_on='pred_overlap',
                #                   right_on='name',
                #                   suffixes=['_df1', '_df2'],
                #                   how='outer')  # assuming you want to keep unmatched records
        # no_overlap_preds = process_assembly_predictions(rows)
        # filtered_rbps = pd.concat([filtered_rbps, no_overlap_preds])
        # filtered_rbps = filtered_rbps.append([no_overlap_preds],ignore_index=True)
        # for p in no_overlap_preds:
        #     n = pd.DataFrame([p])
        #     n.insert(0, 'assembly', ass)
        #     filtered_rbps = pd.concat([filtered_rbps, n], ignore_index=True)
        #     # print('zzzz')
    filtered_rbps.to_csv(
        "/storage/Documents/service/externe/lcfortier/20240130_cdiff_assembly/RBP/blastout_tblastn_qcov80_besthit.fixcoord.filtered.tsv",
        sep="\t")
    # for pred in no_overlap_preds:
    #     pred.to_csv("/storage/Documents/service/externe/lcfortier/20240130_cdiff_assembly/RBP/test.tsv", sep="\t",mode='a')

    # print(no_overlap_preds)

    # tmp_rbps = rows.groupby(
    #         (
    #             (rows.sstart >= rows.sstart.shift())
    #             &
    #             (rows.send <= rows.send.shift())
    #         ).cumsum().apply(lambda x: x)
    #     )
    # # x = tmp_rbps.apply(lambda x: x.sort_values(by="sstart",ascending=False))
    #
    # print(tmp_rbps)
    # #.apply(lambda clustering_df: clustering_df.sort_values(by=['sstart'])))
    # # tmp_rbps = (rows.agg({'sstart': min, 'send': max}))
    # #     # .merge(df1.item, left_on='i_max', right_index=True, how='left')
    # #     # .drop(columns='i_max')
    # #     # .reset_index(drop=True)
    # # filtered_rbps.append(tmp_rbps)
    # # print(tmp_rbps.to_csv("/storage/Documents/service/externe/lcfortier/20240130_cdiff_assembly/RBP/test.tsv", sep="\t"))
    #
    # for xx in x:
    #     print(f"{xx}")
    #     # print(f"{y}")


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


def process_assembly_predictions(rows):
    non_overlapping = []
    for ind, pred in rows.iterrows():
        overlapping = False
        for index, no in enumerate(non_overlapping):
            if check_overlap(pred, no):
                overlapping = True
                break

        if not overlapping:
            non_overlapping.append(pred)

    best_nonoverlapping = pd.DataFrame()
    for index, no1 in enumerate(non_overlapping):
        pred_same_acc_no_overlap = [no2 for no2 in non_overlapping if
                                    no1['sacc'] == no2['sacc'] and not check_overlap(no1, no2)]
        for index, same_acc_no_overlap in enumerate(pred_same_acc_no_overlap):
            best_nonoverlapping = pd.concat([best_nonoverlapping, same_acc_no_overlap.to_frame().T])

        pred_same_acc_with_overlap = [no2 for no2 in non_overlapping if
                                      no1['sacc'] == no2['sacc'] and check_overlap(no1, no2)]
        best_one = None
        for index, same_acc in enumerate(pred_same_acc_with_overlap):
            if no1['pident'] <= same_acc['pident']:
                best_one = same_acc

        # best_nonoverlapping.append(best_one,ignore_index=True)
        if best_one is not None:
            best_nonoverlapping = pd.concat([best_nonoverlapping, best_one.to_frame().T])
        # best_nonoverlapping.loc[len(best_nonoverlapping)] = best_one
        #remove_best_one_from non_overlapping
    x = best_nonoverlapping.drop_duplicates()
    return x


if __name__ == '__main__':
    argParser = argparse.ArgumentParser()

    # mandatory
    argParser.add_argument(
        "-r",
        "--blast_report",
        help="blast output report with assembly label",
        required=True
    )

    args = argParser.parse_args()

    filtre_rbp(
        args.blast_report
    )
