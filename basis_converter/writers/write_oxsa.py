"""
writers/write_oxsa.py
Basis Set Converter — OXSA writer

Generates a pre-configured run_oxsa_amares.m MATLAB script from a
list of core struct dicts. The script is ready to run in MATLAB with
OXSA installed — metabolite names, PPM positions, linewidth bounds,
and relative amplitude scales are all filled in from the basis set.

The user only needs to:
    1. Set oxsaRoot and subjectFile paths in the script
    2. Run the script in MATLAB

Entry point:
    write_oxsa(basis_list, outpath, seqte, metablist_path) → .m file
"""

import os
import numpy as np
from datetime import datetime

from writers.centre import at_centre


############## Linewidth defaults (Hz) #############
LW_INIT = {
    'MM': 30, 'MM09': 30,
    'Glu': 8, 'Gln': 8, 'GABA': 8, 'GSH': 8, 'Glc': 8,
    'mI': 8, 'ml': 8, 'Tau': 8, 'PE': 8, 'Asc': 8, 'Asp': 8,
    'Lac': 7,
}
LW_BOUNDS = {
    'MM': [5, 150], 'MM09': [5, 150],
    'Glu': [1, 60], 'Gln': [1, 60], 'GABA': [1, 60],
    'GSH': [1, 60], 'Glc': [1, 60], 'mI': [1, 60], 'ml': [1, 60],
    'Tau': [1, 60], 'PE': [1, 60], 'Asc': [1, 60], 'Asp': [1, 60],
}
DEFAULT_LW_INIT   = 6
DEFAULT_LW_BOUNDS = [1, 40]

############## Chemical shift search tolerance (ppm) #############
CHEM_TOL = {
    'MM': 0.25, 'MM09': 0.25,
    'Glu': 0.12, 'Gln': 0.12, 'GABA': 0.12, 'GSH': 0.12, 'Glc': 0.12,
    'mI': 0.10, 'ml': 0.10, 'Tau': 0.10, 'PE': 0.10,
    'Asc': 0.10, 'Asp': 0.10, 'Lac': 0.08,
}
DEFAULT_CHEM_TOL = 0.06

############## Relative amplitude scales (NAA = 1.0) #############
REL_SCALE = {
    'NAA': 1.00, 'NAAG': 0.12,
    'Cr391': 0.40, 'CrNo391': 0.40,
    'PCr393': 0.20, 'PCrNo393': 0.20,
    'Cho': 0.16, 'GPC': 0.10, 'PCh': 0.16,
    'mI': 0.30, 'ml': 0.30, 'sI': 0.08, 'sl': 0.08,
    'Glu': 0.35, 'Gln': 0.18, 'GABA': 0.10,
    'GSH': 0.08, 'Tau': 0.10, 'Lac': 0.05,
    'Asc': 0.06, 'Asp': 0.06, 'Gly': 0.06,
    'Glc': 0.06, 'PE': 0.06,
    'MM': 0.20, 'MM09': 0.15,
}
DEFAULT_REL_SCALE = 0.05

PPM_CALIB = 4.65


def write_oxsa(basis_list, outpath,
               seqte=26.0,
               metablist_path=None):
    """
    Write a pre-configured run_oxsa_amares.m script.

    Parameters
    ----------
    basis_list     : list of core struct dicts
    outpath        : str, full path for output .m file
    seqte          : float, echo time in ms
    metablist_path : str or None, path to metabList.mat
    """
    basis_list = at_centre(basis_list)     # peaks are read off the spectrum at 4.65
    outpath = os.path.abspath(outpath)
    os.makedirs(os.path.dirname(outpath) or '.', exist_ok=True)

    if not basis_list:
        raise RuntimeError("basis_list is empty")

    # Get metabolite names and dominant PPM positions
    from writers.write_vespa import (_load_metablist, _get_lit_ppms, _peaks_at_ppms,
                                     METAB_PPM_TABLE, PPM_CALIB as VC)

    ppm_table = _load_metablist(metablist_path) if metablist_path \
                else {k: [(p, 1) for p in v] for k, v in METAB_PPM_TABLE.items()}

    # Build per-metabolite info
    metab_info = []
    for core in basis_list:
        name      = core.get('name', 'unknown')
        lit_lines = _get_lit_ppms(name, ppm_table)
        # Dominant PPM: the literature line where the basis is strongest, or a fallback
        peaks = _peaks_at_ppms(np.asarray(core['fid']).ravel(), float(core.get('sf') or 123.26),
                               float(core.get('sw') or 4000.0), lit_lines,
                               ppm_calib=float(core.get('ppmCalib') or VC)) if lit_lines else []
        dom_ppm = max(peaks, key=lambda p: p[1])[0] if peaks else _strongest_ppm(core, VC)

        lw_init   = LW_INIT.get(name, DEFAULT_LW_INIT)
        lw_bounds = LW_BOUNDS.get(name, DEFAULT_LW_BOUNDS)
        chem_tol  = CHEM_TOL.get(name, DEFAULT_CHEM_TOL)
        rel_scale = REL_SCALE.get(name, DEFAULT_REL_SCALE)

        metab_info.append({
            'name'     : name,
            'dom_ppm'  : dom_ppm,
            'lw_init'  : lw_init,
            'lw_bounds': lw_bounds,
            'chem_tol' : chem_tol,
            'rel_scale': rel_scale,
        })

    # Get shared params from first metabolite
    first     = basis_list[0]
    sf        = float(first.get('sf') or 123.26)
    sw        = float(first.get('sw') or 4000.0)
    ppm_calib = float(first.get('ppmCalib') or VC)
    source    = first.get('source', 'unknown')
    timestamp = datetime.now().strftime('%Y-%m-%dT%H:%M:%S')

    ############## Build metabolite list string #############
    metab_names_str = "{\n"
    for m in metab_info:
        metab_names_str += f"    '{m['name']}', ...\n"
    metab_names_str += "}"

    ############## Build components struct entries #############
    comp_entries = []
    for m in metab_info:
        cs_hz      = round((m['dom_ppm'] - ppm_calib) * sf, 4)
        cs_tol_hz  = round(m['chem_tol'] * sf, 4)
        lw_lo      = m['lw_bounds'][0]
        lw_hi      = m['lw_bounds'][1]
        entry = f"""
    %% {m['name']}
    c = struct();
    c.metab              = '{m['name']}';
    c.componentName      = '{m['name']}';
    c.ppmExpected        = {m['dom_ppm']:.6f};
    c.ppmInit            = {m['dom_ppm']:.6f};
    c.chemShiftHz        = {cs_hz:.4f};
    c.chemShiftBoundsHz  = [{cs_hz - cs_tol_hz:.4f}, {cs_hz + cs_tol_hz:.4f}];
    c.linewidthInitHz    = {m['lw_init']};
    c.linewidthBoundsHz  = [{lw_lo}, {lw_hi}];
    c.amplitudeInit      = {m['rel_scale']:.4f};
    c.amplitudeBounds    = [0, {round(m['rel_scale'] * 5, 4):.4f}];
    c.phaseInitDeg       = 0;
    c.phaseBoundsDeg     = [-180, 180];
    components(end+1)    = c;"""
        comp_entries.append(entry)

    comp_block = '\n'.join(comp_entries)

    ############## Write MATLAB script #############
    script = f'''%% run_oxsa_amares.m
%  Pre-configured OXSA AMARES fitting script
%  Generated by MRS Basis Set Converter on {timestamp}
%  Source basis set: {os.path.basename(source)}
%  TE = {seqte} ms  |  SF = {sf:.4f} MHz  |  SW = {sw:.1f} Hz
%  Carrier PPM = {ppm_calib:.4f}
%
%  INSTRUCTIONS:
%    1. Set oxsaRoot to your OXSA installation folder
%    2. Set subjectFile to your MRS data file (.mat MARSS format)
%    3. Set outputDir for results
%    4. Run this script in MATLAB
%
%  Reference: Purvis et al., PLOS ONE 2017, doi:10.1371/journal.pone.0185356

clear; clc; close all;

%%%%%%%%%%%% USER SETTINGS %%%%%%%%%%%%
oxsaRoot    = '/path/to/OXSA';           % <-- SET THIS
subjectFile = '/path/to/your/data.mat';  % <-- SET THIS (MARSS .mat format)
outputDir   = './oxsa_results/';         % <-- SET THIS
showPlot    = true;

%%%%%%%%%%%% SPECTRAL PARAMETERS (from basis set) %%%%%%%%%%%%
sf         = {sf:.6f};    % Larmor frequency (MHz)
sw         = {sw:.1f};     % Spectral width (Hz)
carrierPpm = {ppm_calib:.4f};  % Water reference (ppm)
seqte      = {seqte:.1f};     % Echo time (ms)

%%%%%%%%%%%% FITTING OPTIONS %%%%%%%%%%%%
options.apodization = 1;
options.lineshape   = 'L';
options.TolFun      = 1e-12;
options.TolX        = 1e-12;
options.MaxIter     = 2000;
options.MaxFunEvals = 1e6;

%%%%%%%%%%%% OXSA STARTUP %%%%%%%%%%%%
if exist(fullfile(oxsaRoot, 'startup.m'), 'file')
    run(fullfile(oxsaRoot, 'startup.m'));
    rehash toolboxcache;
else
    error('Could not find OXSA startup.m at: %s', oxsaRoot);
end

%%%%%%%%%%%% LOAD SUBJECT DATA %%%%%%%%%%%%
S = load(subjectFile);
if ~isfield(S, 'exptDat')
    error('Subject MAT file must contain exptDat struct (MARSS format)');
end
exptDat = S.exptDat;

fid  = double(exptDat.fid(:));
npts = double(exptDat.nspecC);
fid  = fid(1:npts);

% Scale and phase correct
scaleFactor = max(abs(fid));
fid = fid / scaleFactor;
fid = fid .* exp(-1j * angle(fid(1)));

dt       = 1 / sw;
timeAxis = ((0:npts-1)\' * dt) + dt;
freqHz   = ((0:npts-1)\' - npts/2) * (sw / npts);
ppmAxis  = carrierPpm + freqHz / sf;
specReal = real(fftshift(fft(fid)));

exptParams.samples          = npts;
exptParams.imagingFrequency = sf;
exptParams.dwellTime        = dt;
exptParams.beginTime        = dt;
exptParams.timeAxis         = timeAxis;
exptParams.ppmAxis          = freqHz / sf;

fprintf('Loaded: SF=%.4f MHz  SW=%.1f Hz  N=%d\\n', sf, sw, npts);

if showPlot
    figure(\'Color\',\'w\',\'Name\',\'Subject spectrum\');
    plot(ppmAxis, specReal, \'k\', \'LineWidth\', 1);
    set(gca,\'XDir\',\'reverse\'); xlim([0.5 4.5]);
    xlabel(\'ppm\'); ylabel(\'Real\'); title(\'Subject spectrum\');
end

%%%%%%%%%%%% BUILD COMPONENTS %%%%%%%%%%%%
%  Pre-configured from basis set: {len(metab_info)} metabolites
%  Dominant PPM, linewidth bounds, and relative scales are set automatically.

components = struct( ...
    \'metab\', {{}}, \'componentName\', {{}}, \'ppmExpected\', {{}}, ...
    \'ppmInit\', {{}}, \'chemShiftHz\', {{}}, \'chemShiftBoundsHz\', {{}}, ...
    \'linewidthInitHz\', {{}}, \'linewidthBoundsHz\', {{}}, ...
    \'amplitudeInit\', {{}}, \'amplitudeBounds\', {{}}, ...
    \'phaseInitDeg\', {{}}, \'phaseBoundsDeg\', {{}} );

{comp_block}

fprintf(\'Built %d components\\n\', numel(components));

%%%%%%%%%%%% REFINE AMPLITUDE INIT FROM SPECTRUM %%%%%%%%%%%%
%  Estimate local peak amplitude from the actual spectrum for each component.
for i = 1:numel(components)
    ppm0 = components(i).ppmExpected;
    tol  = 0.10;
    mask = ppmAxis >= ppm0-tol & ppmAxis <= ppm0+tol;
    if any(mask)
        localMax = max(specReal(mask));
        if isfinite(localMax) && localMax > 0
            components(i).amplitudeInit = max(localMax / (npts/2), 0.001);
            components(i).amplitudeBounds(2) = max(components(i).amplitudeInit * 5, 0.01);
        end
    end
end

%%%%%%%%%%%% BUILD PRIOR KNOWLEDGE %%%%%%%%%%%%
fields.Bounds = {{'peakName','chemShift','linewidth','amplitude','phase', ...
                  'chemShiftDelta','amplitudeRatio'}};
fields.IV     = {{'peakName','chemShift','linewidth','amplitude','phase'}};
fields.PK     = {{'peakName','multiplet','chemShiftDelta','amplitudeRatio', ...
                  'G_linewidth','G_amplitude','G_phase','RelPhase', ...
                  'G_chemShiftDelta','refPeak'}};

nComp = numel(components);
values.boundsCellArray = cell(nComp, 7);
values.IVCellArray     = cell(nComp, 5);
values.PKCellArray     = cell(nComp, 10);

for i = 1:nComp
    c = components(i);
    values.boundsCellArray(i,:) = {{c.componentName, c.chemShiftBoundsHz, ...
        c.linewidthBoundsHz, c.amplitudeBounds, c.phaseBoundsDeg, [], []}};
    values.IVCellArray(i,:) = {{c.componentName, c.chemShiftHz, ...
        c.linewidthInitHz, c.amplitudeInit, c.phaseInitDeg}};
    values.PKCellArray(i,:) = {{c.componentName, [], [], [], [], [], [], 1, [], 0}};
end

pk = AMARES.priorKnowledge.preparePriorKnowledge(fields, values);
pk.notes = sprintf('Generated by basis_set_converter | TE=%.1f ms | SF=%.4f MHz', ...
                   seqte, sf);

%%%%%%%%%%%% RUN AMARES %%%%%%%%%%%%
fprintf(\'\\nRunning OXSA AMARES...\\n\');
[fitResults, fitStatus, figHandle, CRBResults] = ...
    AMARES.amaresFit(fid, exptParams, pk, showPlot, options);

%%%%%%%%%%%% PRINT RESULTS %%%%%%%%%%%%
chemShiftHz_fit = fitResults.chemShift(:);
ppm_fit         = carrierPpm + chemShiftHz_fit / sf;
lw_fit          = fitResults.linewidth(:);
amp_fit         = fitResults.amplitude(:) * scaleFactor;
phase_fit       = fitResults.phase(:);

fprintf(\'\\n=== AMARES Results ===\\n\');
fprintf(\'%-15s %8s %10s %12s %10s\\n\', \'Peak\', \'ppm\', \'LW(Hz)\', \'Amplitude\', \'Phase\');
fprintf(\'%s\\n\', repmat(\'-\',1,60));
for i = 1:nComp
    fprintf(\'%-15s %8.4f %10.3f %12.4g %10.2f\\n\', ...
        components(i).componentName, ppm_fit(i), lw_fit(i), amp_fit(i), phase_fit(i));
end

%%%%%%%%%%%% SAVE RESULTS %%%%%%%%%%%%
if ~exist(outputDir, \'dir\'), mkdir(outputDir); end

componentNames = string({{components.componentName}}\');
resultTable = table(componentNames, ppm_fit, lw_fit, amp_fit, phase_fit, ...
    \'VariableNames\', {{\'Component\',\'ppm\',\'linewidthHz\',\'amplitude\',\'phaseDeg\'}});

save(fullfile(outputDir, \'oxsa_fitResults.mat\'), ...
    \'fitResults\',\'fitStatus\',\'CRBResults\',\'pk\', ...
    \'exptParams\',\'components\',\'resultTable\',\'scaleFactor\');

writetable(resultTable, fullfile(outputDir, \'oxsa_results.csv\'));
fprintf(\'\\nResults saved to: %s\\n\', outputDir);
'''

    with open(outpath, 'w', encoding='utf-8') as f:
        f.write(script)

    print(f"  Written: {os.path.basename(outpath)}")
    print(f"  {len(metab_info)} metabolites configured")
    print(f"  SF={sf:.4f} MHz  SW={sw:.1f} Hz  TE={seqte} ms")
    print(f"\n  Done. OXSA script: {os.path.basename(outpath)}")


############## Command line ##############

if __name__ == '__main__':
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    if len(sys.argv) < 3:
        print("Usage: python write_oxsa.py <input_folder> <output.m> [te_ms] [metabList.mat]")
        sys.exit(1)

    from readers.read_marss import read_marss_file, read_marss_folder
    input_path     = sys.argv[1]
    outpath        = sys.argv[2]
    seqte          = float(sys.argv[3]) if len(sys.argv) > 3 else 26.0
    metablist_path = sys.argv[4] if len(sys.argv) > 4 else None

    basis_list = read_marss_folder(input_path) if os.path.isdir(input_path) \
                 else [read_marss_file(input_path)]
    print(f"Loaded {len(basis_list)} metabolites\n")
    write_oxsa(basis_list, outpath, seqte=seqte, metablist_path=metablist_path)


def _strongest_ppm(core, default_calib):
    """Where a basis function without literature lines is strongest, between 0.5 and 4.3 ppm."""
    fid = np.asarray(core['fid']).ravel()
    sf, sw = float(core.get('sf') or 123.26), float(core.get('sw') or 4000.0)
    ppm = float(core.get('ppmCalib') or default_calib) + np.fft.fftshift(np.fft.fftfreq(len(fid), 1.0 / sw)) / sf
    spec = np.abs(np.fft.fftshift(np.fft.fft(fid)))
    window = (ppm > 0.5) & (ppm < 4.3)
    return float(ppm[window][np.argmax(spec[window])])
