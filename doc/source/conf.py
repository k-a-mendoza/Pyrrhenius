project = 'Pyrrhenius'
copyright = '2024, Kevin A. Mendoza'
author = 'Kevin A. Mendoza'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['numpydoc',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx_copybutton',
    "jupyter_sphinx",
    'sphinx.ext.todo',
]

templates_path = ['_templates']
exclude_patterns = []

numfig = True
math_numfig = True

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

