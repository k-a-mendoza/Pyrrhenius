import pandas as pd
from . import model
class Database:
    """
    Master database for electric conductivity data

    >> location_of_database = 'database.csv'
    >> ECData = Database(location_of_database)
    >> ECData.load_models()


    """


    def __init__(self,csv):
        self.database = pd.read_csv(csv,encoding_errors= 'replace')
        self.models = {}

    def load_models(self):
        for phase in self.database['Phase Type'].unique():
            self.models[phase]={}
            self.make_single_models(phase)

            self.make_composite_models(phase)

    def make_single_models(self, phase):
        for index, row in self.database[(self.database['Phase Type'] == phase)].iterrows():
            self.models[phase][row['Author ID']] = model.Model(row)

    def make_composite_models(self, phase):
        authors = self.database[(self.database['Phase Type'] == phase)]['Author ID'].unique()
        anisotropic_authors = filter(lambda x: '[' in x, authors)
        group_models = {}
        for author in anisotropic_authors:
            author_stripped = author[:-4]
            if author_stripped not in group_models.keys():
                group_models[author_stripped] = []
            group_models[author_stripped].append(author)
        for k, v in group_models:
            self.models[phase][k + '_iso'] = model.CompositeModel([self.models[phase][x] for x in v])
            
    def get_publication_list(self):
        return self.database['Author ID'].unique()
    
    
    def get_models_from_pubs(self,pub):
        return False

