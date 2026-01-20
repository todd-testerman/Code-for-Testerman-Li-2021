Maaslin2\_NAFLD\_Metagenomes
================
Todd Testerman
6/16/2021

``` r
library(Maaslin2)

input_data_path_abun = read.delim("all_samples_merged_outputs/humann_pathabundance_relab_all_samples.tsv", sep = "\t", row.names = 1) 
input_data_path_abun = data.frame(t(input_data_path_abun), check.names = F)
input_data_path_abun

input_metadata_path_abun = read.delim("metadata_Maaslin.txt", sep = "\t", row.names = 1)
input_metadata_path_abun = data.frame(t(input_metadata_path_abun), check.names = F)
input_metadata_path_abun


input_data_gene_family = read.delim("all_samples_merged_outputs/humann_genefamilies_relab_all_samples.tsv", sep = "\t", row.names = 1) 
input_data_gene_family = data.frame(t(input_data_gene_family), check.names = F)
input_data_gene_family

input_metadata_gene_family = read.delim("metadata_Maaslin.txt", sep = "\t", row.names = 1)
input_metadata_gene_family = data.frame(t(input_metadata_gene_family), check.names = F)
input_metadata_gene_family
```

``` r
fit_data_path_abun = Maaslin2(
    input_data = input_data_path_abun, 
    input_metadata = input_metadata_path_abun, 
    output = "NAFLD_DA_Output_heatmap", 
    fixed_effects = "NAFLD", 
    plot_heatmap = T
)
```

``` r
fit_data_gene_family = Maaslin2(
    input_data = input_data_gene_family, 
    input_metadata = input_metadata_gene_family, 
    output = "NAFLD_Gene_Family_DA_Output", 
    fixed_effects = "NAFLD",
    cores = 16
)
```

``` r
sessionInfo()
```
