# Ring Experiment Analysis Summary

Generated: 2025-07-20 21:43:32.051652

## Experiment Overview

- Network sizes tested: 50, 100
- Lambda values tested: 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95
- Perturbation types: const_high
- Total parameter combinations: 22

## Best Performing Parameters

### By F1 Score:

|   N|perturbation_type | lambda_rounded| f1_mean| mcc_mean| recall_mean| precision_mean|
|---:|:-----------------|--------------:|-------:|--------:|-----------:|--------------:|
|  50|const_high        |           0.01|   0.981|    0.980|       0.964|              1|
| 100|const_high        |           0.01|   0.989|    0.988|       0.978|              1|
