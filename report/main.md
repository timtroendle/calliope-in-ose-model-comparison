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

The storage capacities in Table @tbl:installed-storage-capacities-baseline.

```table
---
caption: 'Installed storage capacities in GWh. {#tbl:installed-storage-capacities-baseline}'
include: ../build/output/baseline/storage-capacity-publish.csv
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

In the low cost scenario, battery cost are 50% of baseline cost. Table @tbl:installed-capacities-low-cost shows the difference in installed capacities in Germany to the baseline.

```table
---
caption: 'Installed capacities in Germany [GW], diff to baseline. {#tbl:installed-capacities-low-cost}'
include: ../build/output/low-cost/capacity-diff.csv
markdown: True
alignment: LR
---
```

And the diff in storage capacities in Germany in Table @tbl:installed-storage-capacities-low-cost.

```table
---
caption: 'Installed storage capacities in Germany [GWh], diff to baseline. {#tbl:installed-storage-capacities-low-cost}'
include: ../build/output/low-cost/storage-capacity-diff.csv
markdown: True
alignment: LR
---
```

## Germany only baseline scenario

In the Germany only baseline scenario, Germany has no connections to partner countries. Table @tbl:installed-capacities-baseline-germany shows the difference in installed capacities in Germany to the baseline.

```table
---
caption: 'Installed capacities in Germany [GW], diff to baseline. {#tbl:installed-capacities-baseline-germany}'
include: ../build/output/baseline-germany/capacity-diff.csv
markdown: True
alignment: LR
---
```

And the diff in storage capacities in Germany in Table @tbl:installed-storage-capacities-baseline-germany.

```table
---
caption: 'Installed storage capacities in Germany [GWh], diff to baseline. {#tbl:installed-storage-capacities-baseline-germany}'
include: ../build/output/baseline-germany/storage-capacity-diff.csv
markdown: True
alignment: LR
---
```

## Germany only low cost scenario

In the Germany only low cost scenario, Germany has no connections to partner countries and battery costs are 50%. Table @tbl:installed-capacities-low-cost-germany shows the difference in installed capacities in Germany to the baseline.

```table
---
caption: 'Installed capacities in Germany [GW], diff to baseline. {#tbl:installed-capacities-low-cost-germany}'
include: ../build/output/low-cost-germany/capacity-diff.csv
markdown: True
alignment: LR
---
```

And the diff in storage capacities in Germany in Table @tbl:installed-storage-capacities-low-cost-germany.

```table
---
caption: 'Installed storage capacities in Germany [GWh], diff to baseline. {#tbl:installed-storage-capacities-low-cost-germany}'
include: ../build/output/low-cost-germany/storage-capacity-diff.csv
markdown: True
alignment: LR
---
```


## Levelised costs

Here are the levelised costs in all scenarios:

```table
---
caption: 'Levelised cost in [€/MWh]. {#tbl:levelised-cost}'
include: ../build/output/cost.csv
markdown: True
alignment: LR
---
```

# Bibliography
