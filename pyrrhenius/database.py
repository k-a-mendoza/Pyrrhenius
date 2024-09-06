import os 
import numpy as np
import pandas as pd
from . import model, mechanisms

DEFAULT_DATABASE_PATH = os.path.join('..', 'database', 'publication_database.csv')

def calc_average_pressure(row):
    """Calculate the Average Pressure of a model entry

    Parameters
    ----------
    row : pd.Series
        pandas series describing a model

    Returns
    -------
    modified_row
        a modified row with a new pressure_average_gpa key-value pair. 
    """
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
    The Database class is designed to manage and manipulate a master database of electric conductivity data.


    
    Attributes
    ----------
    isotropic_name : str
        A class-level attribute that assigns the name suffix for isotropic models.
    database : pd.DataFrame
        A pandas DataFrame that holds the loaded and processed database entries.
    
    Methods
    -------
    __init__(csv=DEFAULT_DATABASE_PATH)
        Initializes the Database object by loading data from a CSV file.
    
    create_isotropic_models()
        Creates isotropic models from the existing database entries and registers them.
    
    register_new_model(ecmodel)
        Registers a new model or a list of models into the database.
    
    get_model_names()
        Returns a list of unique model entry IDs in the database.
    
    get_model_properties(entry_id)
        Returns the properties of a model given its entry ID.
    
    get_model(entry_id)
        Returns the model object for a given entry ID.
    
    get_phases()
        Returns a list of unique phase types in the database.
    
    get_model_list_for_phase(phase)
        Returns a list of model entry IDs for a given phase type.
    
    Examples
    --------
    >>> database = Database()
    >>> phase_list = database.get_phases()
    >>> # pick the phase  you want to use
    >>> model_id_list = database.get_model_list_for_phase('olivine')
    >>> # look at the model_id's you could use
    >>> ecmodel = database.get_model(model_id_list[0])
    """
    isotropic_name = '_isotropic'

    def __init__(self, csv_file=DEFAULT_DATABASE_PATH):
        """
        Initializes the Database object by loading data from a CSV file.

        Parameters
        ----------
        csv_file : str, optional
            The path to the CSV file containing the database data. Defaults to the DEFAULT_DATABASE_PATH, the database shipped with Pyrrhenius. 
        """
        database = pd.read_csv(csv_file, encoding_errors='replace').dropna(how='all')

        data_rows = [model.create_model_from_row(x).get_row() for i, x in database.iterrows() ]
        self.database = pd.DataFrame(data_rows)

    def create_isotropic_models(self):
        """
        Creates isotropic models from the existing database entries and registers them.

        This method identifies groups of models that are identical except for their orientation,
        and creates isotropic mixtures from these groups. It then registers these new isotropic models.

        """
        self.database['grouping_id'] = self.database['entry_id'].apply(lambda row: create_group(row))
        new_models = []
        for _, g in self.database.groupby(['grouping_id']):
            if len(g)>1 and 'isotropic' not in _:
                new_model = model.IsotropicMixture([x['ec_model'] for index, x in g.iterrows()])
                new_models.append(new_model)
        self.register_new_model(new_models)

    def register_new_model(self, ecmodel: model.ModelInterface or list):
        """
        Registers a new model or a list of models into the database.

        Parameters
        ----------
        ecmodel : model.ModelInterface or list
            The model or list of models to be registered.
        """
        if isinstance(ecmodel, list):
            self.database = pd.concat([self.database] + [x.get_row().to_frame().T for x in ecmodel],
                                      ignore_index=True)
        else:
            db_row = ecmodel.get_row().to_frame().T
            self.database = pd.concat([self.database, db_row], ignore_index=True)

    def get_model_names(self):
        """
        Returns a list of unique model entry IDs in the database.

        Returns
        -------
        list[str]
            A list of unique model entry IDs.
        """
        return self.database['entry_id'].unique()

    def get_model_properties(self, entry_id):
        """
        Returns the properties of a model given its entry ID.

        Parameters
        ----------
        entry_id : str
            The entry ID of the model.

        Returns
        -------
        properties : pd.Series
            A pandas Series containing the properties of the model.
        """
        return self.database[self.database['entry_id'] == entry_id]

    def get_model(self, entry_id):
        """
        Returns the electric conductivity model object for a given entry ID.

        Parameters
        ----------
        entry_id : str
            The entry ID of the model.

        Returns
        -------
        ec_model : model.ModelInterface
            The electric conductivity model object.
        """
        return self.database[self.database['entry_id'] == entry_id]['ec_model'].values[0]


    def get_phases(self):
        """
        Returns a list of unique phase types in the database.

        Returns
        -------
        list[str]
            A list of unique phase types.
            
        """
        assert 'ec_model' in self.database.columns, "please load database with .load_data()"
        return list(self.database['phase_type'].unique())

    def get_model_list_for_phase(self, phase):
        """
        Returns a list of model entry IDs for a given phase type.

        Parameters
        ----------
        phase : str
            The phase type.

        Returns
        ------- 
        list[str]
            A list of model entry IDs for the given phase type.
        """
        assert 'ec_model' in self.database.columns, "please load database with .load_data()"
        return list(self.database[self.database['phase_type'] == phase]['entry_id'])

