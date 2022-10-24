Data description
================

``` r
library(Biostrings)
library(here)
library(tidyverse)
```

# brassicaceae_proteomes.rda

This object contains a list of proteomes for Brassicaceae species. FASTA
files containing raw proteomes were downloaded from BrassicaDB. Included
species:

| Species                         | Assembly                  | Reference            |
|:--------------------------------|:--------------------------|:---------------------|
| A.alpina_cv.Pajares             | GCA 000733195.1 v4        | Willing et al, 2015  |
| A.arabicum                      | Aarabicum v3.1            | Nguyen et al, 2019   |
| A.lyrata_cv.MN47                | Alyrata 384 v2.1          | Rawat et al, 2015    |
| A.thaliana                      | Athaliana 447 Araport11   | Cheng et al, 2017    |
| B.carinata_cv.zd-1              | GCA 016771965.1           | Song et al, 2021     |
| B.juncea_cv.AU213               | GCA 020002505.1           | Yang et al, 2021     |
| B.juncea_cv.T84-66              | GCA 020002515.1           | Yang et al, 2021     |
| B.napus_cv.ZS11                 | GCA 000686985.2           | Sun et al, 2017      |
| B.nigra_cv.NI100                | CGI Ni100 ONT Assembly v2 | Perumal et al, 2020  |
| B.nigra_cv.Sangam               | GCA 016432835.1           | Paritosh et al, 2020 |
| B.oleracea_cv.HDEM              | GCA 900416815.2           | Belser et al, 2018   |
| B.oleracea_cv.JZS               | GWHAASO00000000 v2.0      | Cai et al, 2020      |
| B.rapa_cv.Chiifu                | GCA 000309985.3           | Zhang et al, 2018    |
| B.rapa_cv.FPsc                  | Brapa FPsc 277 v1.3       | Nordberg et al, 2014 |
| B.retrofracta                   | GCA 015832515.1           | Kliver et al, 2018   |
| B.stricta                       | Bstricta 278 v1.2         | Lee et al, 2017      |
| C.grandiflora                   | Cgrandiflora 266 v1.1     | Slotte et al, 2013   |
| C.hirsuta                       | Chirsuta MPIPZ v1.0       | Gan et al, 2016      |
| C.rubella_cv.MonteGargano       | Crubella 474 v1.1         | Slotte et al, 2013   |
| C.sativa_cv.DH55                | GCA 000633955.1           | Kagale et al, 2014   |
| E.salsugineum                   | Esalsugineum 173 v1.0     | Yang et al, 2013     |
| R.raphanistrum_ssp.landra       | GWHANWL00000000 v1.0      | Zhang et al, 2021    |
| R.raphanistrum_ssp.raphanistrum | GWHANWM00000000 v1.0      | Zhang et al, 2021    |
| R.sativus_cv.XYB36-2            | Rapsa_Xiang v1.0          | Zhang et al, 2015    |
| S.parvula                       | Sparvula 574 v2.2         | Oh et al, 2014       |
| V.vinifera_cv.PN40024           | Vvinifera 457 v2.1        | Jaillon et al, 2007  |

``` r
# Load proteomes
library(Biostrings)
proteomes <- list.files(here("data", "Proteomes"), full.names = TRUE)
seqs <- lapply(proteomes, readAAStringSet)

# Define function to get primary transcripts only
get_primary_transcripts <- function(seq_stringset = NULL) {
    n <- names(seq_stringset)
    
    # For each transcript, find the encoding gene
    unique_seqs <- gsub("^.*gene=", "", n, ignore.case = TRUE)
    unique_seqs <- gsub("\\s.*$", "", unique_seqs)
    unique_seqs <- gsub("\\.t?[0-9]+(.p)?$", "", unique_seqs)
    
    df_clean <- data.frame(
        names = n, 
        clean_names = unique_seqs, 
        len = width(seq_stringset)
    ) %>% 
        group_by(clean_names) %>% 
        filter(len == max(len))
    
    df_clean <- df_clean[!duplicated(df_clean$clean_names), ]
    
    # Return the AAStringset datatype filtered by the remaining names
    return(seq_stringset[df_clean$names])
}

# Filter seqs to keep only primary transcripts
brassicaceae_proteomes <- lapply(seq_along(seqs), function(x) {
    clean_seq <- get_primary_transcripts(seqs[[x]])
    return(clean_seq)
})

# Save object
save(
    brassicaceae_proteomes,
    file = here("data", "brassicaceae_proteomes.rda"),
    compress = "xz"
)
```
