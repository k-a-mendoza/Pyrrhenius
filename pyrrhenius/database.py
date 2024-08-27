import os 
import numpy as np
import pandas as pd
from . import model, mechanisms

DEFAULT_DATABASE_PATH = os.path.join('..', 'database', 'publication_database.csv')



def _id_row_rename(col):
    col = col.replace(' ', '_').replace('(', '').replace(')', '')
    return col.lower()


def calc_average_pressure(row):
    if pd.isna(row['pressure_average_gpa']):
        return (row['pressure_min_gpa'] + row['pressure_max_gpa']) / 2
    return row['pressure_average_gpa']


def create_group(x):
    for x1 in [0,1]:
        for y1 in [0,1]:
            for z1 in [0,1]:
                x=x.replace(f'[{x1}{y1}{z1}]','[]')
    return x
class Database:
    """
    Master database for electric conductivity data

    >> location_of_database = 'database.csv'
    >> ECData = Database(location_of_database)
    >> ECData.load_models()


    """
    isotropic_name = '_isotropic'

    def __init__(self, csv=DEFAULT_DATABASE_PATH):
        database = pd.read_csv(csv, encoding_errors='replace').dropna(how='all')

        data_rows = [model.create_model_from_row(x).get_row() for i, x in database.iterrows() ]
        self.database = pd.DataFrame(data_rows)

    def create_isotropic_models(self):
        self.database['grouping_id'] = self.database['entry_id'].apply(lambda row: create_group(row))
        new_models = []
        for _, g in self.database.groupby(['grouping_id']):
            if len(g)>1 and 'isotropic' not in _:
                new_model = model.IsotropicMixture([x['ec_model'] for index, x in g.iterrows()])
                new_models.append(new_model)
        self.register_new_model(new_models)

    def register_new_model(self, ecmodel: model.ModelInterface or list):
        if isinstance(ecmodel, list):
            self.database = pd.concat([self.database] + [x.get_row().to_frame().T for x in ecmodel],
                                      ignore_index=True)
        else:
            db_row = ecmodel.get_row().to_frame().T
            self.database = pd.concat([self.database, db_row], ignore_index=True)

    def get_model_names(self):
        return self.database['entry_id'].unique()

    def get_model_properties(self, entry_id):
        return self.database[self.database['entry_id'] == entry_id]

    def get_model(self, entry_id):
        return self.database[self.database['entry_id'] == entry_id]['ec_model'].values[0]


    def get_phases(self):
        assert 'ec_model' in self.database.columns, "please load database with .load_data()"
        return list(self.database['phase_type'].unique())

    def get_model_list_for_phase(self, phase):
        assert 'ec_model' in self.database.columns, "please load database with .load_data()"
        return list(self.database[self.database['phase_type'] == phase]['entry_id'])

