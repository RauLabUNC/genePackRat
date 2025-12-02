# locusPackRat 0.4.0

## New Features
- Overlap detection in `addRatTable()` with `overlap_mode` parameter for region-based linking
- `removeRatTable()` function for removing supplementary tables from projects
- Column abbreviation feature for Excel output prefixes (avoids column name collisions)
- Legacy genome support: hg19 and mm10 now available alongside hg38 and mm39
- `queryMouseMine()` now uses httr for more reliable API communication
- Automatic formatting for MouseMine and OpenTargets data in Excel sheet outputs

## Bug Fixes
- Row stripping in `makeGeneSheet()` only occurs with >= 2 rows present
- Fix for handling multiple regions in `generateLocusZoomPlot()`
- Rare error fix in OpenTargets query response handling

## Documentation
- Updated QTL mapping vignette
- Improved function documentation and column name descriptions

---

# locusPackRat 0.3.0

- Initial support for locus zoom plot generation via plotgardener
- Query functions for external databases (MouseMine, OpenTargets)
- Project-based workflow with persistent storage

---

# locusPackRat 0.2.0

- Project directory structure generation
- Vignette for basic workflow

---

# locusPackRat 0.1.0

- Initial package release
- Core functions: `initPackRat()`, `addRatTable()`, `makeGeneSheet()`
- Support for human and mouse genomes (hg38, mm39)
