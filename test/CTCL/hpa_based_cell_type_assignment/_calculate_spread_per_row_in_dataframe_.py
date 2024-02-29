import pandas as pd
def _calculate_spread_per_row_in_dataframe_(dataframe):
    list_of_spread_scores_allColumns_per_row=[]
    column_values_desc_for_all_rows=[]
    for index in dataframe.index:
        row=dataframe.loc[[index]]
        update_df = dataframe.drop(index)
        row_complement=pd.DataFrame(update_df.max()).T.rename(index={0: index})
        spread_scores_allColumns_for_index=row.sub(row_complement) # Spread score (x-max(-x)) for index=index, across all columns present in HPA data.
        index_with_highest_spread=spread_scores_allColumns_for_index.index[0]
        spread_scores_allColumns_for_index_desc_dict=spread_scores_allColumns_for_index[spread_scores_allColumns_for_index.iloc[-1,:].sort_values(ascending=False).index].T[index_with_highest_spread].to_dict()
        spread_scores_allColumns_for_index_desc_list=list(spread_scores_allColumns_for_index_desc_dict.keys())
        len_spread_scores_allColumns_for_index_desc_list=len(spread_scores_allColumns_for_index_desc_list)
        list_of_spread_scores_allColumns_per_row.append(spread_scores_allColumns_for_index_desc_list)
        column_values_desc_for_all_rows.append(list(spread_scores_allColumns_for_index_desc_dict.values()))
    return len_spread_scores_allColumns_for_index_desc_list, list_of_spread_scores_allColumns_per_row, column_values_desc_for_all_rows