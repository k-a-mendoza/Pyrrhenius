import numpy as np
import pandas as pd
from . import model, mechanisms
def _id_row_rename(col):
    col = col.replace(' ', '_').replace('(', '').replace(')', '')
    return col.lower()

class Database:
    """
    Master database for electric conductivity data

    >> location_of_database = 'database.csv'
    >> ECData = Database(location_of_database)
    >> ECData.load_models()


    """

    def __init__(self, csv):
        self.database = pd.read_csv(csv, encoding_errors='replace').dropna(how='all')
        self.database.rename(_id_row_rename,inplace=True,axis=1)
        self.models = {}

    def load_models(self):
        self.database['ec_model']= self.database.apply(lambda row : create_model_from_row(row),axis=1)
        self.create_anisotropic_models()


    def create_anisotropic_models(self):
        subframe = self.database[self.database['crystal_direction']!='isotropic']
        subframe['grouping_id'] = subframe['entry_id'].str.slice(stop=-5)
        groups = list(subframe.groupby(['grouping_id'])['ec_model'])
        for g in groups:

            base_row= dict(vars(g[1].values[0].metadata))
            base_row['entry_id'] = g[0]+'_anisomodel'
            base_row['ec_model'] = model.AnisotropicModel(g[1].values)
            appending_series= pd.Series(base_row)
            self.database = pd.concat([self.database,appending_series.to_frame().T],ignore_index=True,axis=0)

    def get_model_names(self):
        return self.database['entry_id'].unique()

    def get_model_properties(self,entry_id):
        return self.database[self.database['entry_id']==entry_id]

    def get_phases(self):
        assert 'ec_model' in self.database.columns,"please load database with .load_data()"
        return list(self.database['phase_type'].unique())

    def get_model_list_for_phase(self,phase):
        assert 'ec_model' in self.database.columns, "please load database with .load_data()"
        return list(self.database[self.database['phase_type']==phase]['entry_id'])
def _create_multiple_mechanisms(row):
    potential_mechanisms = [x.strip() for x in row['eq_id'].split('+')]
    constants = create_constants_from_row(row)[::-1]
    mech_list = []
    for pot_mech in potential_mechanisms:
        n_constants = mechanisms.model_dict[pot_mech].n_constants
        specific_constants = [constants.pop() for n in range(n_constants)]
        mech_list.append(mechanisms.model_dict[pot_mech](*specific_constants))
    return mech_list

def _create_single_mechanism(row):
    mechanism = row['eq_id'].strip()
    target_mechanism = mechanisms.model_dict[mechanism]
    constants        = create_constants_from_row(row)

    assert len(constants) == target_mechanism.n_constants, "Incorrect constant number defined for mechanism. "+\
        "Should have: " +str(target_mechanism.n_constants) + " but found " + str(len(constants)) +" in file\n"+\
                                                            "Check database entry "+row['entry_id']+" against "+\
                                                            str(target_mechanism)
    return [target_mechanism(*constants)]
def create_mechanism_from_row(row) -> list:
    mechanism        = row['eq_id']
    if '+' in mechanism:
        return _create_multiple_mechanisms(row)
    else:
        return _create_single_mechanism(row)

def get_const_rows(row):
    letter_rows = filter(lambda x: len(x) == 1,row.keys())
    return sorted(filter(lambda x: len(x) == 1 and ~np.isnan(row[x]), row.keys()),key=lambda x : x.lower())
def create_constants_from_row(row: pd.Series):
    letter_columns = get_const_rows(row)
    try:
        variables = [mechanisms.StochasticConstant(row[f'{x}'],
                                               row[f'{x}_uncertainty'],
                                               row[f'{x}_description'].split('_')[0],
                                               'log' in row[f'{x}_description'])
                 for x in letter_columns]
    except Exception as e:
        print(f'problems initializing ' + row['entry_id'])
        print(e)
        variables = None
    return variables

def create_model_from_row(row : pd.Series):
    mechs = create_mechanism_from_row(row)
    metadata_row_keys = list(set(row.keys()) -  set(get_const_rows(row)))
    return model.Model(mechs, row[metadata_row_keys])