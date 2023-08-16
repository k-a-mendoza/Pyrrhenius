import pandas as pd


class Model:
    def __init__(self,pandas_row):
        self.pandas_row = pandas_row


class CompositeModel:

    def __init__(self, models):
        self.models = models