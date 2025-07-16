# Ring Experiment Analysis Summary

Generated: 2025-07-16 17:39:29.282367

## Experiment Overview

- Network sizes tested: 50, 100
- Lambda values tested: 0.01, 0.03, 0.05, 0.07, 0.09, 0.1, 0.15, 0.2, 0.229, 0.25, 0.265, 0.3, 0.307, 0.356, 0.4, 0.413, 0.468, 0.478, 0.5, 0.505, 0.545, 0.554, 0.588, 0.6, 0.634, 0.643, 0.65, 0.684, 0.7, 0.738, 0.745, 0.75, 0.796, 0.8, 0.85, 0.859, 0.9, 0.929, 0.95, 0.953, 0.963, 0.964, 0.971, 0.975, 0.976, 0.979, 0.981, 0.982, 0.984, 0.985, 0.987, 0.989, 0.99, 0.991, 0.992
- Perturbation types: const_high, const_low, lambda_double, lambda_half
- Total parameter combinations: 78

## Best Performing Parameters

### By F1 Score:

|   N|perturbation_type | lambda_rounded| f1_mean| mcc_mean| recall_mean| precision_mean|
|---:|:-----------------|--------------:|-------:|--------:|-----------:|--------------:|
|  50|const_high        |           0.99|   0.958|    0.956|       0.920|              1|
| 100|const_high        |           0.99|   0.962|    0.960|       0.926|              1|
