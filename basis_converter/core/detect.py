"""
core/detect.py
Basis Set Converter — format detection

Given a list of file paths, figures out what format they are,
reads one file to extract real metadata, and returns it for display.

Returns a dict with:
    'status'            : 'full' | 'partial' | 'ambiguous' | 'failed'
    'format'            : string name, e.g. 'LCModel .BASIS'
    'ambiguous_options' : list of possible format strings (only if ambiguous)
    'data'              : dict — core struct with real extracted values
    'missing'           : list of field names not found
    'message'           : human-readable explanation
"""

import os
import sys
from core.formats import get_compatible_tools

# Add parent dir to path so readers are importable
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def detect_format(paths):
    result = _detect_format_inner(paths)
    # Add compatible tools based on detected format
    fmt = result.get('format', '')
    result['compatible_tools'] = get_compatible_tools(fmt)
    return result


def _detect_format_inner(paths):
    if not paths:
        return _failed("No files provided.")

    # Filter hidden/system files and macOS resource files
    # Catches: .DS_Store, ._filename, Thumbs.db, desktop.ini
    # NOTE: extension-less files are allowed — SpinWizard/JET uses
    # files with no extension where filename = metabolite name
    IGNORED_NAMES  = {'Thumbs.db', 'desktop.ini'}
    IGNORED_NO_EXT = {'Icon'}   # macOS Icon file (has no extension)
    paths = [
        p for p in paths
        if not os.path.basename(p).startswith('.')
        and os.path.basename(p).strip() not in IGNORED_NAMES
        and not (os.path.splitext(p)[1] == ''
                 and os.path.basename(p).strip() in IGNORED_NO_EXT)
    ]

    if not paths:
        return _failed("No usable files found (only hidden or system files).")

    extensions = set(os.path.splitext(p)[1].lower() for p in paths)
    n_files    = len(paths)

    if n_files == 1 and extensions == {'.basis'}:
        return _try_lcmodel_basis(paths[0])

    if n_files == 1 and paths[0].lower().endswith('.nii.gz'):
        return _try_niftimrs(paths[0])

    ############## Folder of .nii.gz files (with optional .json sidecars) #############
    # spec2nii writes a .json sidecar alongside each .nii.gz — treat as NIfTI-MRS
    nii_files  = [p for p in paths if p.lower().endswith('.nii.gz')]
    json_files = [p for p in paths if p.lower().endswith('.json')]
    other      = [p for p in paths
                  if not p.lower().endswith('.nii.gz')
                  and not p.lower().endswith('.json')]
    if nii_files and not other:
        return _try_niftimrs_folder(nii_files)

    if n_files == 1 and extensions == {'.mat'}:
        return _try_single_mat(paths[0])

    if extensions == {'.json'}:
        return _try_fsLmrs(paths)

    # Filter out LCMBasisSet.txt before jMRUI/GAVA check
    txt_files = [p for p in paths
                 if p.lower().endswith('.txt')
                 and os.path.basename(p) != 'LCMBasisSet.txt']
    if txt_files and not [p for p in paths if not p.lower().endswith('.txt')
                          and os.path.basename(p) != 'LCMBasisSet.txt']:
        # Try GAVA first (single .txt file with semicolons + tabs)
        if len(txt_files) == 1:
            result = _try_gava(txt_files[0])
            if result['status'] != 'failed':
                return result
        return _try_jmrui(txt_files)
    if extensions == {'.txt'} and not any(
            os.path.basename(p) == 'LCMBasisSet.txt' for p in paths):
        if len(paths) == 1:
            result = _try_gava(paths[0])
            if result['status'] != 'failed':
                return result
        return _try_jmrui(paths)

    ############## Single .csv file — PyAMARES prior knowledge #############
    if n_files == 1 and extensions == {'.csv'}:
        return _try_pyamares(paths[0])

    ############## Folder with no-extension files — SpinWizard / JET #############
    # Case 1: folder path passed directly (Browse folder)
    # Only try SpinWizard if the folder contains extension-less files
    if n_files == 1 and os.path.isdir(paths[0]):
        folder_files = [
            f for f in os.listdir(paths[0])
            if not f.startswith('.') and not f.startswith('._')
            and os.path.isfile(os.path.join(paths[0], f))
        ]
        SPINWIZARD_IGNORE = {'LCMBasisSet.txt', 'BasisSetParameters.txt',
                              'BasisRFOffset.txt'}
        no_ext_in_folder = [f for f in folder_files
                            if os.path.splitext(f)[1] == ''
                            and f not in SPINWIZARD_IGNORE]
        if no_ext_in_folder:
            result = _try_spinwizard(paths)
            if result['status'] != 'failed':
                return result
    # Case 2: files passed individually (Browse files) — any no-extension files present
    SW_IGNORE = {'LCMBasisSet.txt', 'BasisSetParameters.txt', 'BasisRFOffset.txt'}
    no_ext_files = [p for p in paths
                    if os.path.splitext(p)[1] == '' and os.path.isfile(p)
                    and os.path.basename(p) not in SW_IGNORE]
    if no_ext_files:
        # Pass all paths so companion files can be found in the same folder
        result = _try_spinwizard(paths)
        if result['status'] != 'failed':
            return result

    ############## Single .xml file — MIDAS or VeSPA VIFF prior #############
    if n_files == 1 and extensions == {'.xml'}:
        result = _try_midas(paths[0])
        if result['status'] != 'failed':
            return result
        return _try_vespa(paths[0])

    ############## Single .txt file — GAVA or jMRUI #############
    # (already handled above for jMRUI folders; this is for single GAVA files)

    if extensions == {'.mat'}:
        return _try_mat_folder(paths)

    if extensions == {'.raw'}:
        return _try_raw_folder(paths)

    if len(extensions) > 1:
        return _try_mixed(paths)

    return _failed("File format not recognized.")


############## Detectors — each reads one file for real metadata #############

def _try_lcmodel_basis(path):
    try:
        from readers.read_lcmodel import read_lcmodel_basis
        results = read_lcmodel_basis(path)
        if not results:
            return _failed("No metabolites found in .BASIS file.")
        r = results[0]
        # Round sw to avoid floating point precision issues (e.g. 3999.9998)
        sw = round(r['sw'])
        return {
            'status':           'full',
            'format':           'LCModel .BASIS',
            'ambiguous_options': [],
            'data': {
                'sw':    sw,
                'sf':    r['sf'],
                'n':     r['n'],
                'Bo':    round(r['sf'] / 42.577, 3) if r['sf'] else None,
                'names': [x['name'] for x in results],
                'te':    r.get('te'),
                'seq':   r.get('seq'),
            },
            'missing': [] if r.get('te') else ['te'],
            'message': 'LCModel .BASIS file detected.',
        }
    except Exception as ex:
        return _failed(f"Could not read .BASIS file: {ex}")


def _try_niftimrs_folder(paths):
    """Folder of .nii.gz files — read first for metadata."""
    try:
        from readers.read_niftimrs import read_niftimrs
        # Prefer _conj files if present
        conj = [p for p in paths if '_conj' in os.path.basename(p)]
        first = conj[0] if conj else paths[0]
        results = read_niftimrs(first)
        if not results:
            return _failed("No data found in NIfTI-MRS files.")
        r = results[0]
        # Get all names from filenames
        names = [
            os.path.splitext(os.path.splitext(
                os.path.basename(p))[0])[0].replace('_conj', '')
            for p in paths
            if '_conj' in os.path.basename(p) or
               not any('_conj' in os.path.basename(q)
                       for q in paths)
        ]
        return {
            'status':           'full',
            'format':           'NIfTI-MRS',
            'ambiguous_options': [],
            'data': {
                'sw':    r['sw'],
                'sf':    r.get('sf'),
                'n':     r['n'],
                'Bo':    round(r['sf'] / 42.577, 3) if r.get('sf') else None,
                'names': names,
            },
            'missing': [],
            'message': 'NIfTI-MRS basis set folder detected.',
        }
    except Exception as ex:
        return _failed(f"Could not read NIfTI-MRS folder: {ex}")


def _try_niftimrs(path):
    try:
        from readers.read_niftimrs import read_niftimrs
        results = read_niftimrs(path)
        if not results:
            return _failed("No data found in NIfTI-MRS file.")
        r = results[0]
        return {
            'status':           'full',
            'format':           'NIfTI-MRS',
            'ambiguous_options': [],
            'data': {
                'sw':    r['sw'],
                'sf':    r.get('sf'),
                'n':     r['n'],
                'Bo':    round(r['sf'] / 42.577, 3) if r.get('sf') else None,
                'names': [x['name'] for x in results],
            },
            'missing': [],
            'message': 'NIfTI-MRS file detected.',
        }
    except Exception as ex:
        return _failed(f"Could not read NIfTI-MRS file: {ex}")


def _struct_keys(raw, keys):
    """The fields of a .mat file's one variable, when it holds just one struct (MRSCloud
    saves each metabolite as outA); none otherwise."""
    if len(keys) != 1:
        return set()
    val = raw[next(iter(keys))]
    if isinstance(val, dict):
        return set(val.keys())
    return set(getattr(getattr(val, 'dtype', None), 'names', None) or [])


def _try_single_mat(path):
    """
    Inspect a single .mat file's fields to determine format.
    Priority: Osprey → INSPECTOR → MARSS → FID-A → MRSCloud
    Handles both scipy (v5/v7) and mat73 (v7.3 HDF5) formats.
    """
    raw = None
    try:
        import scipy.io as sio
        raw = sio.loadmat(path, simplify_cells=True)
    except NotImplementedError:
        # v7.3 HDF5 file — use mat73
        try:
            import mat73
            raw = mat73.loadmat(path)
        except ImportError:
            return _failed(
                f"Cannot read {os.path.basename(path)}: MATLAB v7.3 file. "
                f"Install mat73: pip install mat73"
            )
        except Exception as ex:
            return _failed(f"Could not load .mat file: {ex}")
    except Exception as ex:
        return _failed(f"Could not load .mat file: {ex}")

    if raw is None:
        return _failed(f"Could not load .mat file: {os.path.basename(path)}")

    keys = set(k for k in raw if not k.startswith('_'))

    ############## ProFit: has 'basis' key with cell array fid + f0 field #############
    # Check before Osprey since both use a 'basis' key
    if 'basis' in keys:
        b = raw.get('basis', {})
        b_keys = set(b.keys()) if isinstance(b, dict) else set()
        import numpy as np
        fid_val = b.get('fid') if isinstance(b, dict) else None
        is_cell = fid_val is not None and np.asarray(fid_val).dtype == object
        if {'fid', 'met', 'bw', 'f0'}.issubset(b_keys) and is_cell:
            return _read_and_return_profit(path)

    ############## Osprey: has 'BASIS' key with fids, spectralwidth, name #############
    if 'BASIS' in keys or 'basis' in keys:
        return _read_and_return_osprey(path)

    ############## INSPECTOR: has 'lcmBasis' key with data, sw_h, sf #############
    if 'lcmBasis' in keys or 'lcmbasis' in keys:
        return _read_and_return_inspector(path)

    ############## INSPECTOR fallback: data, sw_h, sf at top level #############
    if {'data', 'sw_h', 'sf'}.issubset(keys):
        return _read_and_return_inspector(path)

    ############## ProFit: has 'basis' key with fid (cell array), met, bw, f0 #############
    if 'basis' in keys:
        b = raw.get('basis', {})
        b_keys = set(b.keys()) if isinstance(b, dict) else set()
        if {'fid', 'met', 'bw', 'f0'}.issubset(b_keys):
            return _read_and_return_profit(path)

    ############## ProFit fallback: top-level fid (cell), met, bw, f0 #############
    if {'fid', 'met', 'bw', 'f0'}.issubset(keys):
        fid_val = raw.get('fid')
        if hasattr(fid_val, 'dtype') and fid_val.dtype == object:
            return _read_and_return_profit(path)

    # ############# MARSS: has 'exptDat' with fid, sw_h, sf, nspecC #############
    if 'exptDat' in keys:
        expt = raw['exptDat']
        expt_keys = set(expt.keys()) if isinstance(expt, dict) else set(expt.dtype.names or [])
        if 'fid' in expt_keys and 'sw_h' in expt_keys:
            return _read_and_return_marss_single(path)

    ############## FID-A: has 'out' or top-level fids + txfrq + nucleus #############
    if 'out' in keys:
        out = raw['out']
        out_keys = set(out.keys()) if isinstance(out, dict) else set(out.dtype.names or [])
        if 'fids' in out_keys and 'txfrq' in out_keys:
            return _read_and_return_fida_single(path)

    ############## MRSCloud: fids + txfrq + Bo, top-level or in the one struct it saves (outA) #############
    if {'fids', 'txfrq', 'Bo'}.issubset(keys | _struct_keys(raw, keys)):
        return _read_and_return_mrscloud_single(path)

    # Unknown — return ambiguous with whatever we can infer
    return {
        'status':           'ambiguous',
        'format':           'Unknown .mat',
        'ambiguous_options': ['Osprey .mat', 'INSPECTOR .mat',
                              'MARSS .mat',  'FID-A .mat', 'MRSCloud .mat'],
        'data':             {},
        'missing':          [],
        'message':          'Single .mat file — could not identify format from fields.',
    }


def _try_mat_folder(paths):
    """
    Folder of .mat files — inspect the first file to narrow down.
    Could be MARSS, FID-A, or MRSCloud.
    """
    try:
        import scipy.io as sio
        raw = sio.loadmat(paths[0], simplify_cells=True)
    except NotImplementedError:
        try:
            import mat73
            raw = mat73.loadmat(paths[0])
        except Exception:
            raw = {}
    except Exception:
        raw = {}

    keys = set(k for k in raw if not k.startswith('_'))

    # MARSS: exptDat with fid, sw_h, sf, nspecC
    if 'exptDat' in keys:
        expt = raw['exptDat']
        expt_keys = set(expt.keys()) if isinstance(expt, dict) else set(expt.dtype.names or [])
        if {'fid', 'sw_h', 'sf'}.issubset(expt_keys):
            return _read_and_return_marss_folder(paths)

    # FID-A: out struct with fids + txfrq + nucleus
    if 'out' in keys:
        out = raw['out']
        out_keys = set(out.keys()) if isinstance(out, dict) else set(out.dtype.names or [])
        if {'fids', 'txfrq', 'nucleus'}.issubset(out_keys):
            return _read_and_return_fida_folder(paths)

    # MRSCloud: fids + txfrq + Bo + name, top-level or in the one struct it saves (outA)
    if {'fids', 'txfrq', 'Bo'}.issubset(keys | _struct_keys(raw, keys)):
        return _read_and_return_mrscloud_folder(paths)

    # Still ambiguous — return with partial data from first file
    return {
        'status':           'ambiguous',
        'format':           'MARSS .mat',
        'ambiguous_options': ['MARSS .mat', 'FID-A .mat', 'MRSCloud .mat'],
        'data':             {},
        'missing':          [],
        'message':          'Folder of .mat files — could be MARSS, FID-A, or MRSCloud.',
    }


def _try_raw_folder(paths):
    """
    Folder of .raw files — detect whether MARSS or LCModel format.
    MARSS .raw has no $SEQPAR block — sw and sf will be missing.
    LCModel .raw has $SEQPAR with HZPPPM, NUNFIL, DELTAT.
    """
    try:
        from readers.read_lcmodel import read_lcmodel_raw
        r = read_lcmodel_raw(paths[0])
        names = [os.path.splitext(os.path.basename(p))[0] for p in paths]
        is_marss = r.get('is_marss_raw', False)

        if is_marss:
            # MARSS .raw — no spectral params in file
            return {
                'status':           'full',
                'format':           'MARSS .raw',
                'ambiguous_options': [],
                'data': {
                    'sw':    None,
                    'sf':    None,
                    'n':     r['n'],
                    'Bo':    None,
                    'names': names,
                },
                'missing': ['sw', 'sf', 'te'],
                'message': 'MARSS .raw basis functions detected. '
                           'Spectral width and Larmor frequency not stored — '
                           'will be requested on next screen.',
            }
        else:
            # LCModel full .raw — has $SEQPAR
            return {
                'status':           'full',
                'format':           'LCModel .raw',
                'ambiguous_options': [],
                'data': {
                    'sw':    r['sw'],
                    'sf':    r['sf'],
                    'n':     r['n'],
                    'Bo':    round(r['sf'] / 42.577, 3) if r['sf'] else None,
                    'names': names,
                },
                'missing': ['te'],
                'message': 'LCModel .raw basis functions detected.',
            }
    except Exception as ex:
        return {
            'status':           'ambiguous',
            'format':           'LCModel .raw',
            'ambiguous_options': ['LCModel .raw files', 'MARSS .raw files'],
            'data':             {},
            'missing':          [],
            'message':          f'Folder of .raw files detected (metadata read failed: {ex}).',
        }


def _try_fsLmrs(paths):
    try:
        from readers.read_fsLmrs import read_fsLmrs_file
        r = read_fsLmrs_file(paths[0])
        names = [os.path.splitext(os.path.basename(p))[0] for p in paths]
        return {
            'status':           'full',
            'format':           'FSL-MRS .json',
            'ambiguous_options': [],
            'data': {
                'sw':    r['sw'],
                'sf':    r.get('sf'),
                'n':     r['n'],
                'Bo':    None,
                'names': names,
            },
            'missing': ['sf'] if not r.get('sf') else [],
            'message': 'FSL-MRS JSON basis folder detected.',
        }
    except Exception as ex:
        return _failed(f"Could not read FSL-MRS .json file: {ex}")


def _try_jmrui(paths):
    # Verify it's actually jMRUI format
    try:
        with open(paths[0], 'r') as f:
            first_line = f.readline().strip()
        if 'jMRUI' not in first_line:
            return _failed(
                f"First .txt file does not appear to be jMRUI format."
            )
        from readers.read_jmrui import read_jmrui_file
        r = read_jmrui_file(paths[0])
        names = [os.path.splitext(os.path.basename(p))[0] for p in paths]
        return {
            'status':           'full',
            'format':           'jMRUI .txt',
            'ambiguous_options': [],
            'data': {
                'sw':    r['sw'],
                'sf':    r['sf'],
                'n':     r['n'],
                'Bo':    r.get('Bo'),
                'names': names,
            },
            'missing': [],
            'message': 'jMRUI text basis folder detected.',
        }
    except Exception as ex:
        return _failed(f"Could not read jMRUI .txt file: {ex}")


def _read_and_return_profit(path):
    try:
        from readers.read_profit import read_profit
        results = read_profit(path)
        if not results:
            return _failed("No metabolites found in ProFit .mat file.")
        r = results[0]
        return {
            'status':           'full',
            'format':           'ProFit .mat',
            'ambiguous_options': [],
            'data': {
                'sw':    r['sw'],
                'sf':    r['sf'],
                'n':     r['n'],
                'Bo':    r.get('Bo'),
                'names': [x['name'] for x in results],
            },
            'missing': [],
            'message': 'ProFit basis set detected.',
        }
    except Exception as ex:
        return _failed(f"Could not read ProFit .mat: {ex}")


def _try_spinwizard(paths):
    """Try to read SpinWizard / JET basis files (no extension, two-column text)."""
    try:
        import numpy as np
        from readers.read_spinwizard import read_spinwizard_folder, read_spinwizard_file

        # Resolve: could be a folder path, a list of files, or a single file
        if len(paths) == 1 and os.path.isdir(paths[0]):
            # Folder passed directly
            folder = paths[0]
            results = read_spinwizard_folder(folder)
        elif len(paths) == 1 and os.path.isfile(paths[0]):
            # Single file
            test_path = paths[0]
            data = np.loadtxt(test_path)
            if data.ndim == 1:
                data = data.reshape(1, -1)
            if data.shape[1] < 2:
                return _failed("Not a two-column file")
            results = [read_spinwizard_file(test_path)]
        else:
            # Multiple files — exclude companion/list files
            SW_IGNORE = {'LCMBasisSet.txt', 'BasisSetParameters.txt',
                         'BasisRFOffset.txt'}
            file_paths = [p for p in paths
                          if os.path.isfile(p)
                          and os.path.basename(p) not in SW_IGNORE]
            if not file_paths:
                return _failed("No files found")
            data = np.loadtxt(file_paths[0])
            if data.ndim == 1:
                data = data.reshape(1, -1)
            if data.shape[1] < 2:
                return _failed("Not a two-column file")

            # Read companion files from the same folder
            from readers.read_spinwizard import _read_companion_files
            folder   = os.path.dirname(file_paths[0])
            companion= _read_companion_files(folder)
            sw_c     = companion['sw']
            sf_c     = companion['sf']
            pc_c     = companion['ppmCalib']

            results = [read_spinwizard_file(p, sw=sw_c, sf=sf_c,
                                            ppm_calib=pc_c)
                       for p in file_paths]

        if not results:
            return _failed("No SpinWizard files found")

        names     = [r['name'] for r in results]
        n         = results[0]['n']
        sw        = results[0].get('sw')
        sf        = results[0].get('sf')
        ppm_calib = results[0].get('ppmCalib')
        Bo        = round(sf / 42.577, 3) if sf else None

        # Status is full if sw and sf were found in companion files
        missing = []
        if not sw: missing.append('sw')
        if not sf: missing.append('sf')
        status  = 'full' if not missing else 'partial'

        return {
            'status'           : status,
            'format'           : 'SpinWizard',
            'ambiguous_options': [],
            'data'             : {
                'sw'       : sw,
                'sf'       : sf,
                'n'        : n,
                'Bo'       : Bo,
                'ppmCalib' : ppm_calib,
                'names'    : names,
            },
            'missing'          : missing,
            'message'          : f'SpinWizard / JET basis files detected. '
                                 f'{len(results)} metabolites found.'
                                 + (f' SW={sw:.0f} Hz, SF={sf:.2f} MHz.' if sf else
                                    ' Spectral width and Larmor frequency required.'),
        }
    except Exception as ex:
        return _failed(f"Not a SpinWizard file: {ex}")


def _try_midas(path):
    """Try to read a MIDAS FITT_Generic_XML prior file."""
    try:
        import xml.etree.ElementTree as ET
        tree = ET.parse(path)
        root = tree.getroot()
        if root.tag != 'FITT_Generic_XML':
            return _failed(f"Not a MIDAS file: root tag is '{root.tag}'")

        prior = root.find('PRIOR_METABOLITE_INFORMATION')
        if prior is None:
            return _failed("No PRIOR_METABOLITE_INFORMATION element found")

        params = prior.findall('param')
        if not params:
            return _failed("No metabolite params found")

        # Extract unique metabolite names
        names = []
        for p in params:
            value = p.get('value', '')
            parts = value.split('++')
            if parts and parts[0] not in names:
                names.append(parts[0])

        return {
            'status'           : 'partial',  # sw/sf not in file
            'format'           : 'MIDAS .xml',
            'ambiguous_options': [],
            'data'             : {
                'sw'   : None,
                'sf'   : None,
                'n'    : 2048,
                'Bo'   : None,
                'names': names,
            },
            'missing'          : ['sw', 'sf'],
            'message'          : f'MIDAS prior XML detected. '
                                 f'{len(names)} metabolites found.',
        }
    except Exception as ex:
        return _failed(f"Not a MIDAS file: {ex}")


def _try_gava(path):
    """Try to read a GAVA text prior file."""
    try:
        names = []
        data_lines = 0
        with open(path, 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(';'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 8:
                    name = parts[0].strip()
                    if name and name not in names:
                        names.append(name)
                    data_lines += 1

        if data_lines == 0:
            return _failed("No tab-separated data lines found")

        return {
            'status'           : 'partial',
            'format'           : 'GAVA .txt',
            'ambiguous_options': [],
            'data'             : {
                'sw'   : None,
                'sf'   : None,
                'n'    : 2048,
                'Bo'   : None,
                'names': names,
            },
            'missing'          : ['sw', 'sf'],
            'message'          : f'GAVA text prior detected. '
                                 f'{len(names)} metabolites found.',
        }
    except Exception as ex:
        return _failed(f"Not a GAVA file: {ex}")


def _try_vespa(path):
    """Try to read a VeSPA VIFF XML prior file."""
    try:
        from readers.read_vespa import read_vespa
        results = read_vespa(path)
        if not results:
            return _failed("No metabolites found in VeSPA XML file.")
        names = [r['name'] for r in results]
        r = results[0]
        return {
            'status'           : 'full',
            'format'           : 'VeSPA .xml',
            'ambiguous_options': [],
            'data': {
                'sw'   : r.get('sw'),
                'sf'   : r.get('sf'),
                'n'    : r.get('n'),
                'Bo'   : round(r['sf'] / 42.577, 3) if r.get('sf') else None,
                'names': names,
            },
            'missing': ['sf', 'sw'] if not r.get('sf') else [],
            'message': f'VeSPA Analysis prior file detected. '
                       f'{len(results)} metabolites found.',
        }
    except RuntimeError as ex:
        return _failed(str(ex))
    except Exception as ex:
        return _failed(f"Could not read VeSPA XML: {ex}")


def _try_pyamares(path):
    """Detect and read a PyAMARES prior knowledge CSV file."""
    try:
        from readers.read_pyamares import read_pyamares
        result = read_pyamares(path)
        n_peaks = len(result.get('names', []))
        return {
            'status':           'full',
            'format':           'PyAMARES .csv',
            'ambiguous_options': [],
            'data': {
                'names' : result['names'],
                'sw'    : None,
                'sf'    : None,
                'n'     : None,
                'Bo'    : None,
                # Pass through peak params for write_oxsa
                'peak_params': result,
            },
            'missing': [],
            'message': f'PyAMARES prior knowledge CSV detected ({n_peaks} peaks).',
        }
    except Exception as ex:
        return _failed(f"Could not read PyAMARES CSV: {ex}")


def _try_mixed(paths):
    supported   = [p for p in paths if _is_supported(p)]
    unsupported = [p for p in paths if not _is_supported(p)]
    if not supported:
        return _failed(f"None of the {len(paths)} files could be read.",
                       files_failed=unsupported)
    return {
        'status':           'partial',
        'format':           'Mixed',
        'ambiguous_options': [],
        'data':             {},
        'missing':          [],
        'message':          f"{len(unsupported)} file(s) could not be read and will be skipped.",
    }


############## Per-format read helpers #############

def _read_and_return_osprey(path):
    try:
        from readers.read_osprey import read_osprey
        results = read_osprey(path)
        if not results:
            return _failed("No metabolites found in Osprey .mat file.")
        r = results[0]
        return {
            'status':           'full',
            'format':           'Osprey .mat',
            'ambiguous_options': [],
            'data': {
                'sw':    r['sw'],
                'sf':    r['sf'],
                'n':     r['n'],
                'Bo':    r.get('Bo'),
                'names': [x['name'] for x in results],
            },
            'missing': [],
            'message': 'Osprey .mat BASIS file detected.',
        }
    except Exception as ex:
        return _failed(f"Could not read Osprey .mat: {ex}")


def _read_and_return_inspector(path):
    try:
        from readers.read_inspector import read_inspector
        results = read_inspector(path)
        if not results:
            return _failed("No metabolites found in INSPECTOR .mat file.")
        r = results[0]
        return {
            'status':           'full',
            'format':           'INSPECTOR .mat',
            'ambiguous_options': [],
            'data': {
                'sw':    r['sw'],
                'sf':    r['sf'],
                'n':     r['n'],
                'Bo':    round(r['sf'] / 42.577, 3) if r['sf'] else None,
                'names': [x['name'] for x in results],
            },
            'missing': [],
            'message': 'INSPECTOR .mat basis file detected.',
        }
    except Exception as ex:
        return _failed(f"Could not read INSPECTOR .mat: {ex}")


def _read_and_return_marss_single(path):
    try:
        from readers.read_marss import read_marss_file
        r = read_marss_file(path)
        return {
            'status':           'full',
            'format':           'MARSS .mat',
            'ambiguous_options': [],
            'data': {
                'sw':    r['sw'],
                'sf':    r['sf'],
                'n':     r['n'],
                'Bo':    round(r['sf'] / 42.577, 3),
                'names': [r['name']],
            },
            'missing': [],
            'message': 'Single MARSS .mat basis function detected.',
        }
    except Exception as ex:
        return _failed(f"Could not read MARSS .mat: {ex}")


def _read_and_return_marss_folder(paths):
    try:
        from readers.read_marss import read_marss_file
        r = read_marss_file(paths[0])
        names = [os.path.splitext(os.path.basename(p))[0] for p in paths]
        return {
            'status':           'full',
            'format':           'MARSS .mat',
            'ambiguous_options': ['MARSS .mat', 'FID-A .mat', 'MRSCloud .mat'],
            'data': {
                'sw':    r['sw'],
                'sf':    r['sf'],
                'n':     r['n'],
                'Bo':    round(r['sf'] / 42.577, 3),
                'names': names,
            },
            'missing': [],
            'message': 'Folder of MARSS .mat basis functions detected.',
        }
    except Exception as ex:
        return _failed(f"Could not read MARSS .mat folder: {ex}")


def _read_and_return_fida_single(path):
    try:
        from readers.read_fida import read_fida_file
        r = read_fida_file(path)
        return {
            'status':           'full',
            'format':           'FID-A .mat',
            'ambiguous_options': [],
            'data': {
                'sw':    r['sw'],
                'sf':    r['sf'],
                'n':     r['n'],
                'Bo':    r.get('Bo'),
                'names': [r['name']],
            },
            'missing': [] if r.get('te') else ['te'],
            'message': 'FID-A .mat basis function detected.',
        }
    except Exception as ex:
        return _failed(f"Could not read FID-A .mat: {ex}")


def _read_and_return_fida_folder(paths):
    try:
        from readers.read_fida import read_fida_file
        r = read_fida_file(paths[0])
        names = [os.path.splitext(os.path.basename(p))[0] for p in paths]
        return {
            'status':           'full',
            'format':           'FID-A .mat',
            'ambiguous_options': ['FID-A .mat', 'MARSS .mat', 'MRSCloud .mat'],
            'data': {
                'sw':    r['sw'],
                'sf':    r['sf'],
                'n':     r['n'],
                'Bo':    r.get('Bo'),
                'names': names,
            },
            'missing': [] if r.get('te') else ['te'],
            'message': 'Folder of FID-A .mat basis functions detected.',
        }
    except Exception as ex:
        return _failed(f"Could not read FID-A .mat folder: {ex}")


def _read_and_return_mrscloud_single(path):
    try:
        from readers.read_mrscloud import read_mrscloud_file
        r = read_mrscloud_file(path)
        return {
            'status':           'full',
            'format':           'MRSCloud .mat',
            'ambiguous_options': [],
            'data': {
                'sw':    r['sw'],
                'sf':    r['sf'],
                'n':     r['n'],
                'Bo':    r.get('Bo'),
                'names': [r['name']],
            },
            'missing': [],
            'message': 'MRSCloud .mat basis function detected.',
        }
    except Exception as ex:
        return _failed(f"Could not read MRSCloud .mat: {ex}")


def _read_and_return_mrscloud_folder(paths):
    try:
        from readers.read_mrscloud import read_mrscloud_file
        r = read_mrscloud_file(paths[0])
        names = [os.path.splitext(os.path.basename(p))[0] for p in paths]
        return {
            'status':           'full',
            'format':           'MRSCloud .mat',
            'ambiguous_options': ['MRSCloud .mat', 'FID-A .mat', 'MARSS .mat'],
            'data': {
                'sw':    r['sw'],
                'sf':    r['sf'],
                'n':     r['n'],
                'Bo':    r.get('Bo'),
                'names': names,
            },
            'missing': [],
            'message': 'Folder of MRSCloud .mat basis functions detected.',
        }
    except Exception as ex:
        return _failed(f"Could not read MRSCloud .mat folder: {ex}")


############## Utilities #############

def _is_supported(path):
    ext = os.path.splitext(path)[1].lower()
    return ext in {'.mat', '.basis', '.raw', '.json', '.txt', '.gz'}


def _failed(message, files_failed=None):
    return {
        'status':            'failed',
        'format':            'Unreadable',
        'ambiguous_options': [],
        'data':              {},
        'missing':           [],
        'message':           message,
        'files_failed':      files_failed or [],
        'compatible_tools':  [],
    }
