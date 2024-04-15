def _transfer_values_from_local_to_global_dict_(local_dict, global_dict):
    for key, value in local_dict.items():
        if key in global_dict.keys():
            global_dict[key].append(value)
        else:
            global_dict[key]=[value]
    return global_dict

def _save_dict_of_dataframes_to_csv_files_(dictionary, filename_prefix=None):
    for key, value in dictionary.items():
        if filename_prefix:
            exec_str=f'value.to_csv("{filename_prefix}_{key}.csv")'
        else:
            exec_str=f'value.to_csv("{key}.csv")'
        exec(exec_str)

def _substitute_value_in_list_from_dict_(list_, lookup_dict):
    list_=[lookup_dict.get(item,item)  for item in list_]
    return list_

def _substitute_value_in_dataframe_from_dict_(dataframe, column_name, lookup_dict):
    column=dataframe[column_name]
    column=_substitute_value_in_list_from_dict_(column, lookup_dict)
    dataframe[column_name]=column
    return dataframe