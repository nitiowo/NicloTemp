# 06_phenotype

Summarizes and visualizes phenotype data across treatment groups.

## Scripts

| Script | Description |
|--------|-------------|
| `01_summarize_traits.R` | Aggregates cleaned phenotype data to treatment-group level (mean per temp × niclosamide × infection group). Outputs separate tables for uninfected and infected conditions, plus an individual-level long-format table |
| `02_trait_overview.R` | Visualizes group-level trait summaries as \ heatmaps and per-trait bar/line plots across treatment groups |

## Phenotype Categories

- **Uninfected snail traits:** 
- **Infected snail / parasite traits:** 

## Outputs

- `output/trait_summaries/trait_summaries_uninfected.csv`
- `output/trait_summaries/trait_summaries_infected.csv`
- `output/trait_summaries/trait_summaries_individual.csv`
- `output/figures/` — heatmaps and per-trait plots
