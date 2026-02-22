
## emery (Version 0.6.2)

### Major changes

- New functions: random_start_binary(), unique_obs_summary()
- New accessor functions for MultiMethodMLEstimate objects: getNames(), getResults(), setFreqs()
- estimate_ML_binary() and estimate_ML_ordinal() now (optionally) accept a summary of unique observations and their 
frequencies which can result in significantly increased calculation speeds, especially when bootstrapping large data sets
- New slot added to MultiMethodMLEstimate S4 objects to support above

### Minor changes

- Automated testing added for some functions

## emery (Version 0.6.0)

### Major changes

- New functions: aggregate_boot_ML(), plot.boot_ML, bin_auc

### Minor changes

- Standardized outputs from estimate_ML for data types
- Added warnings to some functions
- Progress bar for boot_ML
- New color palette for graphs

## emery (Version 0.5.1)

### Major changes

- CRAN resubmission.

### Minor changes

- Corrected reference format in DESCRIPTION file
- Added description of return value to `show-MultiMethodMLEstimate-method`
- Other minor documentation text corrections

## emery (Version 0.5.0)

### Major changes

- Initial CRAN submission.

