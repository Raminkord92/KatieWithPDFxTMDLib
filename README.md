# PDFxTMDLib usage in KaTie input files
KatieWithPDFxTMDLib is a clone of Katie Parton level event generator https://bitbucket.org/hameren/katie/src/master/.
This document explains how to use PDFxTMDLib (TMDs and collinear PDFs) from KaTie input files. Example inputs are in `diphoton/`:
- `diphoton/pdfxtmdKt` (TMD / kT-dependent)
- `diphoton/pdfxtmdCol` (collinear / cPDF)

## Setup (library paths)
Edit `settings.py` and set the PDFxTMDLib path before building:
```
# Path to the directory where libPDFxTMDLib.so is located.
# Leave empty if PDFxTMD is not available on this machine.
PDFXTMDpath = '/usr/local/lib'
```

Then generate the library:
```
./config.py lib
```

## Quick start
1) Build KaTie with PDFxTMDLib linked (see `./config.py lib` above).
2) In the input file, set `PDFxTMDSet` and choose either TMD (default) or cPDF mode.
3) Optionally use PDFxTMDLib for alphaQCD with `PDFxTMDCoupling`.

## Input keywords

### Required
- `PDFxTMDSet = <setName> <beamA> <beamB> <member>`
  - `setName`: PDFxTMDLib set name (string)
  - `beamA`, `beamB`: beam IDs (typically `2212` for proton)
  - `member`: set member index (integer, often `0`)

### Optional
- `PDFxTMDType = cpdf`
  - Use collinear PDFs from PDFxTMDLib.
  - If omitted, KaTie uses PDFxTMD TMDs (kT-dependent) by default.
- `PDFxTMDCoupling = yes|no|<setName>`
  - If `yes`, alphaQCD is taken from PDFxTMDLib. If a set name is provided,
    that set is used for the coupling.
- `PDFxTMDCouplingSet = <setName>`
  - Explicit set name for alphaQCD (overrides default behavior).
- `useLHAPDF = no|yes`
  - If `no`, KaTie will skip LHAPDF initialization and rely on PDFxTMDLib.

## Choosing TMD vs collinear (cPDF)

### TMD mode (kT-dependent)
Use off-shell beams and do not set `PDFxTMDType`:
```
PDFxTMDSet = PB-LO-HERAI+II-2020-set1 2212 2212 0
PDFxTMDCoupling = yes
useLHAPDF = no

offshell = 1 1
```

### Collinear mode (cPDF)
Use on-shell beams and set `PDFxTMDType = cpdf`:
```
PDFxTMDSet = MSTW2008lo68cl 2212 2212 0
PDFxTMDType = cpdf
PDFxTMDCoupling = yes
useLHAPDF = no

offshell = 0 0
```

## Notes and tips
- `PDFxTMDSet` is used for both beams by default. If you need different sets
  per beam, use the two-set form supported in your workflow (see local examples).
- If you set `useLHAPDF = no`, KaTie will not call LHAPDF init or LHAPDF PDFs.
- If you set `PDFxTMDType = cpdf`, KaTie will not initialize PDFxTMD TMDs.

## Troubleshooting
- **Error:** `LHAPDF::ReadError Info file not found for PDF set ...`
  - Ensure `useLHAPDF = no` when you only want PDFxTMDLib.
  - Regenerate the run directory after changing input options.
- **Still seeing LHAPDF calls in generated code**
  - Re-run `./run.sh prepare <input> <dir>` so the template is re-expanded.

## Example inputs
See:
- `diphoton/pdfxtmdKt` (TMD)
- `diphoton/pdfxtmdCol` (cPDF)
