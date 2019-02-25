# R ggplot Function for Prettier Manhattan Plots
Author: Pagé Goddard

Source: https://github.com/pcgoddard/Burchardlab_Tutorials/wiki/GGplot2-Manhattan-Plot-Function#use

March 22 2017

## Contents
* [Motivation](#motivation)
* [Usage](#use)
* [Script](#script)
* [Example Plots](#examples)
* [BONUS: How I Choose Color Palettes](#colors)

Based on the code from [R Graph Gallery](https://www.r-graph-gallery.com/wp-content/uploads/2018/02/Manhattan_plot_in_R.html)

second ref: [Getting Genetics Done](http://www.gettinggeneticsdone.com/2010/01/gwas-manhattan-plots-and-qq-plots-using.html)

<a name="motivation"></a>
***Motivation***: *the qqman package in R has a manhattan plot function that is very efficient but limited in customization options. Of particular note, it is difficult to annotate manhattan plots (highlight SNPs or label SNPs below a certain p-value threshold) in a visually pleasing & legible manner. Here I provide a function to make manhattan plots in R with full ggplot customizability.*

<a name="use"></a>
### Usage
1. Format GWAS results as you would for qqman: `SNP` `CHR` `BP` `P` (tab-delim)
2. Define external significance and suggestive thresholds for data and plots (or edit `# add genome-wide sig and sugg lines` section in the function code)
2. Copy and paste the script [below](#script) into your environment then call the function:
```r
# Variables ====
mypalette <- c("#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B") # chr color palette
mysnps <- c("rs11801961","rs116558464","rs61703161") # snps to highlight
sig = 5e-8 # significant threshold line
sugg = 1e-6 # suggestive threshold line

# Define Function ====
# see below

# Run Function ====
gg.manhattan(df, threshold=1e-6, hlight=mysnps, col=mypalette, ylims=c(0,10), title="My Manhattan Plot")
```
* **`df`** = input dataframe: SNP CHR BP P (script is capital sensitive)
* **`threshold`** = numeric value or NA; label snps with pvalue < threshold
    - `ERROR: Aesthetics must be either length 1 or the same as the data (1): label, alpha, x, y` means no snps meet the threshold criteria
    - for no labels, use `NA`
* **`hlight`** = list of SNP names to highlight
    - adjust color inside function
    - for no SNP higlighting, use `NA`
* **`col`** = vector listing color/colors desired for manhattan plot
* **`title`** = plot title, with quotations (e.g. "My Manhattan Plot")
*  additional ggplot arguments (`ylim`, `cex`, etc. can be altered within the plot section of the function)

<a name="script"></a>
### Script
**copy and paste me!**
**customize me!**
```r
# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

gg.manhattan <- function(df, threshold, hlight, col, ylims, title){
  # format df
  df.tmp <- df %>% 
    
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%
    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) %>%
    mutate( is_annotate=ifelse(P < threshold, "yes", "no"))
  
  # get chromosome center positions for x-axis
  axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  ggplot(df.tmp, aes(x=BPcum, y=-log10(P))) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
    scale_color_manual(values = rep(col, 22 )) +

    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
    
    # add plot and axis titles
    ggtitle(paste0(title)) +
    labs(x = "Chromosome") +
    
    # add genome-wide sig and sugg lines
    geom_hline(yintercept = -log10(sig)) +
    geom_hline(yintercept = -log10(sugg), linetype="dashed") +
    
    # Add highlighted points
    #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +
    
    # Custom the theme:
    theme_bw(base_size = 22) +
    theme( 
      plot.title = element_text(hjust = 0.5),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}
```

<a name="examples"></a>
### Example Figures
```r
# with labeled snps

png("prefev1_loco_rsid_180322.png", width=1425, height=975)
print(gg.manhattan(loco_asm_fev, sugg, NA, c(1.5,9), reds.c, "Pre-treatment FEV1 % Predicted - MLMA (LOCO) GWAS Results\n"))
dev.off()
```
<img alt="gg.manhattan plot with labeled snps" src="https://github.com/pcgoddard/Burchardlab_Tutorials/blob/master/prefev1_loco_rsid_180322.png" width="800"/>


```r
# without labeled snps

png("prefev1_loco_180322.png", width=1425, height=975)
print(gg.manhattan(loco_asm_fev, NA, NA, c(1.5,9), reds.c, "Pre-treatment FEV1 % Predicted - MLMA (LOCO) GWAS Results\n"))
dev.off()
```
<img alt="gg.manhattan plot without labeled snps" src="https://github.com/pcgoddard/Burchardlab_Tutorials/blob/master/prefev1_loco_180322.png" width="800"/>

using qqman
```r
library(qqman)

png("prefev1_loco_qqman_rsid_180322.png", width=1425, height=975)
print(manhattan(loco_asm_fev, 
                suggestiveline = -log10(sugg), 
                genomewideline = -log10(sig), 
                col=reds.c, 
                cex=1.1, 
                ylim = c(1.5, 8),
                annotatePval = sugg, annotateTop = FALSE,
                main = "Pre-treatment FEV1 % Predicted - MLMA (LOCO) GWAS Results\nQQman Plot"))
dev.off()
```
<img alt="qqman plot without labeled snps" src="https://github.com/pcgoddard/Burchardlab_Tutorials/blob/master/prefev1_loco_qqman_rsid_180322.png" width="800"/>

<a name="colors"></a>
### BONUS: How I Choose Color Palettes
* [paletton.com](http://paletton.com/#uid=1000u0kllllaFw0g0qFqFg0w0aF) is my go-to resource for making color palettes.
* For the above figures, I listed the hex codes for all 5 shades of the selected color in a vector:
```R
reds <- c("#FF817E", "#E9534F", "#D92B26", "#AE1612", "#870300")
```
* you can adjust color brightness, saturation, and contrast by **shifting the dots** in the middle of the color wheel or toggling the **fine tune** options
* To create three palettes with different tones within one color family, I used the `adjacent 3 colors option` Below are 3 example palettes.
    * [3 blue-hue palettes](http://paletton.com/#uid=53G0p0ks-tYhHBXmJvnuHnX-0iU)
    * [3 red-hue palettes](http://paletton.com/#uid=5000k0kqprfgbzFl8tjsJlOxCgW)
    * [3 green-hue palettes](http://paletton.com/#uid=b2A2U2a0kqprfgbzFl8tjsJlOxCgW)
* to get a save-able html table, go to the **bottom R corner** --> `TABLES/EXPORT` --> `as HTML`

Here is a sample of a complex color palette scheme I have used for manhattan plots in a project that employed multiple phenotypes and multiple GWAS methods:
```R
# w = warmer tones
# n = neutral
# c = cooler tones

# Trait 1
reds.w <- c("#FFAD7E", "#E9874F", "#D96726", "#AE4A12", "#873100") # method 1
reds.n <- c("#FF817E", "#E9534F", "#D92B26", "#AE1612", "#870300") # method 2
reds.c <- c("#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B") # method 3

# Trait 2
blues.w <- c("#4EAFAF", "#2C9696", "#0F8F8F", "#057272", "#005A5A") # method 1
blues.n <- c("#5D82BB", "#3B64A5", "#1E4F9E", "#103B7E", "#082B64") # method 2
blues.c <- c("#6E65C2", "#4D43AE", "#3226A6", "#211785", "#150D69") # method 3

# Trait 3
greens.w <- c("#DDF67A", "#C2E04C", "#ADD024", "#88A711", "#678100") # method 1
greens.n <- c("#A7E672", "#85D047", "#6AC122", "#4F9B10", "#367800") # method 2
greens.c <- c("#A6C965", "#7B9F34", "#567714", "#375201", "#203000") # method 3
```

more sources for colorschemes : [link 1](http://www.sthda.com/english/wiki/colors-in-r),[link 2](https://www.r-bloggers.com/color-palettes-in-r/)
