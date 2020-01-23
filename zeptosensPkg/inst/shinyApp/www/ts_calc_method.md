# Target Score Calculation Method

Target Score currently provided two method in calculation. One listed as _Line by Line_ and the other listed as _Pooled_.

## Line by Line
Calculation line by line limited the calculation by putting number of doses as one and calculate target score for each 
every line for the provided Perturbation Response File dataset.
## Pooled
Calculation Pooled provided the calculation by putting number of doses as number of rows which sum up target score for each 
every line treated as different doses for the provided Perturbation Response File dataset. (ts <- colSums(tsd))