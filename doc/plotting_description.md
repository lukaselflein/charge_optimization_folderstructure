# Plots

## Content
`averages.png`: Charge per atom. The means of the constrained/unconstrained and bader calculations are compared to the charges fitted to the averaged cost functions.
`avg_heatmap.png`: Distance between constrained and unconstrained charges, with different sigma and lnrhoref parameters. Charges are calculated from averaged cost functions.
`boxplot.png`: Like `averages.png`, but showing the full distributions instead of the means. The whiskers are the 100% confidence interval, i.e. the min/max values of the distribution.
`combined_box_swarm_bader.png`: Boxplot of the bader charge distribution, with each snapshot as a color-coded point.
`combined_box_swarm_constrained.png`: Like `combined_box_swarm_bader.png`, but for constrained horton charges.
`combined_box_swarm_unconstrained.png`: Like `combined_box_swarm_bader.png`, but for unconstrained horton charges.
`dispersion.png`: Broadness of the constrained charge distribution, depending on parameters.
`lnrho_vs_diff_sigma_hue.png`: like `avg_heatmap.png`, but with lineplots instead of a heatmap.
`pointplot.png`
`snapshot_heatmap.png`: Like `avg_heatmap.png`, but with charges calculated on individual snapshot costfunctions and then averaged.
