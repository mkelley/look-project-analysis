# look-project-analysis
LCO Outbursting Objects Key Project Analysis Tools

1. Download BANZAI processed files to e91/{date}/
2. If any files need updated WCS, save the updates as FITS header formatted metadata to data/wcs/{label}/ with the extension .hdr.  The header file base name and the BANZAI processed base name must match, save any reduction level suffix (e.g., NEOExchange uses e92).
3. Run `file-summary.py` to find files to process.
phot.py -q &&\
stack.py -q &&\
color-and-ostat.py -q &&\
phot-to-json.py &&\
temporal-filter.py -q
