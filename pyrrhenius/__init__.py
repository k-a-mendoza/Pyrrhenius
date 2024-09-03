import json
import os

# Determine the path to the config.json file
constants_path = os.path.join(os.path.dirname(__file__), 'constants.json')

# Load the configuration
with open(constants_path, 'r') as config_file:
    CONSTANTS = json.load(config_file)