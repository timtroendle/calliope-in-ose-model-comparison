# Introduction

...

# Results

## Baseline scenario

Figure @fig:timeseries-baseline shows the generation and demand time series:

![Generation profiles.](../build/output/baseline/plot.png){ #fig:timeseries-baseline .class}

The capacities that have been installed per country are shown in Table @tbl:installed-capacities-baseline.

```table
---
caption: 'Installed capacities in GW. {#tbl:installed-capacities-baseline}'
include: ../build/output/baseline/capacity-publish.csv
markdown: True
---
```

And here, i.e. in Table @tbl:trade-baseline, the electricity traded between countries.

```table
---
caption: 'Traded electricity GWh. {#tbl:trade-baseline}'
include: ../build/output/baseline/trade.csv
markdown: True
---
```

## Low cost scenario

In the low cost scenario, battery cost are 50% of baseline cost. Table @tbl:installed-capacities-low-cost shows the difference in installed capacities to the baseline.

```table
---
caption: 'Installed capacities in GW, diff to baseline. {#tbl:installed-capacities-low-cost}'
include: ../build/output/capacity-diff.csv
markdown: True
alignment: LR
---
```

And eventually, here are the levelised costs in the two scenarios:

```table
---
caption: 'Levelised cost in [€ct/kWh]. {#tbl:levelised-cost}'
include: ../build/output/cost-diff.csv
markdown: True
alignment: LR
---
```

# Bibliography
