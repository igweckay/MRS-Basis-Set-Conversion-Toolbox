"""
core/formats.py
Basis Set Converter — format and tool compatibility definitions

Central registry mapping each detected format to:
    - The software tools that can use that format
    - A short display label
    - The file extension(s)

Used by:
    detect.py     — adds 'compatible_tools' to detection result
    screen1_upload.py — displays compatible tools in the result card
    screen2_output.py — groups tools by format compatibility
"""

############## Format → compatible tools mapping #############
#
# Key   : detected format string (must match exactly what detect.py returns)
# Value : list of tool names that can read/use that format

FORMAT_TOOLS = {

    ############## MARSS #############
    'MARSS .mat': [
        'MARSS',
    ],

    'MARSS .raw': [
        'MARSS',
        'LCModel',
        'SPANT',
        'TARQUIN',
    ],

    ############## LCModel #############
    'LCModel .BASIS': [
        'LCModel',
        'SPANT',
        'TARQUIN',
    ],

    'LCModel .raw': [
        'LCModel',
        'SPANT',
        'TARQUIN',
    ],

    ############## Osprey #############
    'Osprey .mat': [
        'Osprey',
    ],

    ############## FSL-MRS #############
    'FSL-MRS .json': [
        'FSL-MRS',
    ],

    ############## INSPECTOR #############
    'INSPECTOR .mat': [
        'INSPECTOR',
    ],

    ############## MRSCloud #############
    'MRSCloud .mat': [
        'MRSCloud',
    ],

    ############## FID-A #############
    'FID-A .mat': [
        'FID-A',
    ],

    ############## NIfTI-MRS #############
    'NIfTI-MRS': [
        'FSL-MRS',
        'Osprey',
        'SPANT',
    ],

    ############## jMRUI / NMRScopeB / QUEST #############
    'jMRUI .txt': [
        'jMRUI',
        'QUEST',
        'NMRScopeB',
    ],

    ############## ProFit #############
    'ProFit .mat': [
        'ProFit',
    ],

    ############## PyAMARES / OXSA #############
    'PyAMARES .csv': [
        'PyAMARES',
        'OXSA',
    ],

    ############## VeSPA #############
    'VeSPA .xml': [
        'VeSPA',
    ],
}


def get_compatible_tools(fmt):
    """
    Return the list of compatible tools for a detected format string.
    Returns empty list if format is not in the registry.

    Parameters
    ----------
    fmt : str
        Detected format string e.g. 'LCModel .BASIS'

    Returns
    -------
    list of str
    """
    return FORMAT_TOOLS.get(fmt, [])


def get_tools_display(fmt):
    """
    Return a human-readable string of compatible tools.
    e.g. 'LCModel · SPANT · TARQUIN'

    Parameters
    ----------
    fmt : str

    Returns
    -------
    str
    """
    tools = get_compatible_tools(fmt)
    if not tools:
        return ''
    return '  ·  '.join(tools)


############## Output tools — what screen2_output offers and core.write writes #############
#
# id      : the tool, as core.write.write_basis takes it
# formats : its output formats (ids), if it has more than one

FITTING_TOOLS = [
    {
        'id': 'lcmodel', 'name': 'LCModel', 'ext': '.BASIS / .raw',
        'dimmed': False,
        'formats': [
            {'id': 'both',  'label': '.BASIS + .raw files', 'sub': 'recommended'},
            {'id': 'basis', 'label': '.BASIS file only',    'sub': 'single combined file'},
            {'id': 'raw',   'label': '.raw files only',     'sub': 'one file per metabolite'},
        ]
    },
    {
        'id': 'osprey', 'name': 'Osprey', 'ext': '.mat',
        'dimmed': False,
        'formats': [
            {'id': 'mat',   'label': '.mat file',   'sub': 'Osprey native format'},
            {'id': 'basis', 'label': '.BASIS file', 'sub': 'LCModel-compatible'},
            {'id': 'nifti', 'label': 'NIfTI-MRS',   'sub': '.nii.gz'},
        ]
    },
    {'id': 'fsLmrs',    'name': 'FSL-MRS',   'ext': '.json folder', 'dimmed': False, 'formats': []},
    {'id': 'inspector', 'name': 'INSPECTOR', 'ext': '.mat',         'dimmed': False, 'formats': []},
    {
        'id': 'spant', 'name': 'SPANT', 'ext': '.BASIS / .raw / NIfTI',
        'dimmed': False,
        'formats': [
            {'id': 'basis', 'label': '.BASIS file',  'sub': 'recommended for spant'},
            {'id': 'raw',   'label': '.raw files',   'sub': 'LCModel-compatible folder'},
            {'id': 'nifti', 'label': 'NIfTI-MRS',   'sub': '.nii.gz'},
        ]
    },
    {
        'id': 'abfit', 'name': 'ABfit', 'ext': '.BASIS / .raw / NIfTI',
        'dimmed': False,
        'formats': [
            {'id': 'basis', 'label': '.BASIS file',  'sub': 'recommended for ABfit'},
            {'id': 'raw',   'label': '.raw files',   'sub': 'LCModel-compatible folder'},
            {'id': 'nifti', 'label': 'NIfTI-MRS',   'sub': '.nii.gz'},
        ]
    },
    {'id': 'jmrui',    'name': 'jMRUI',    'ext': '.txt folder', 'dimmed': False, 'formats': []},
    {'id': 'aqses',    'name': 'AQSES',    'ext': '.txt folder', 'dimmed': False, 'formats': []},
    {'id': 'quest',    'name': 'QUEST',    'ext': '.txt folder', 'dimmed': False, 'formats': []},
    {'id': 'tarquin',  'name': 'TARQUIN',  'ext': '.BASIS',      'dimmed': False, 'formats': []},
    {'id': 'profit',   'name': 'ProFit',   'ext': '.mat',        'dimmed': False, 'formats': []},
    {'id': 'niftimrs', 'name': 'NIfTI-MRS','ext': '.nii.gz',     'dimmed': False, 'formats': []},

    {'id': 'jet',      'name': 'JET',       'ext': 'SpinWizard',  'dimmed': False, 'formats': []},
    {'id': 'midas',    'name': 'MIDAS',     'ext': '.xml priors', 'dimmed': False, 'formats': []},
    {'id': 'gava',     'name': 'GAVA',      'ext': '.txt priors', 'dimmed': False, 'formats': []},
    {'id': 'fida_fit', 'name': 'FID-A',     'ext': '.mat',        'dimmed': False, 'formats': []},
    {'id': 'gannet',   'name': 'Gannet',   'ext': 'internal',    'dimmed': True,  'formats': []},
]

GENERATOR_TOOLS = [
    {
        'id': 'marss', 'name': 'MARSS', 'ext': '.mat / .raw',
        'dimmed': False,
        'formats': [
            {'id': 'mat_combined',  'label': '.mat combined',   'sub': 'all metabolites in one file'},
            {'id': 'mat_individual','label': '.mat individual', 'sub': 'one file per metabolite'},
            {'id': 'raw',           'label': '.raw files',      'sub': 'LCModel-compatible'},
        ]
    },
    {'id': 'mrscloud',  'name': 'MRSCloud',  'ext': '.mat',        'dimmed': False, 'formats': []},
    {'id': 'fida',      'name': 'FID-A',     'ext': '.mat',        'dimmed': False, 'formats': []},



    {'id': 'spinwizard','name': 'SpinWizard (Spinach)','ext': 'no extension', 'dimmed': False, 'formats': []},
    {'id': 'vespa_gen2','name': 'VeSPA',     'ext': '.xml priors', 'dimmed': False, 'formats': []},
    {'id': 'nmbscope',  'name': 'NMRScopeB', 'ext': '.txt folder', 'dimmed': False, 'formats': []},
]

PEAK_FITTING_TOOLS = [
    {'id': 'pyamares', 'name': 'PyAMARES', 'ext': '.csv priors',  'dimmed': False, 'formats': []},
    {'id': 'vespa_analysis', 'name': 'VeSPA Analysis', 'ext': '.xml priors', 'dimmed': False, 'formats': []},
    {'id': 'oxsa',     'name': 'OXSA',     'ext': '.m script',    'dimmed': False, 'formats': []},
]

ALL_TOOLS = FITTING_TOOLS + GENERATOR_TOOLS + PEAK_FITTING_TOOLS
