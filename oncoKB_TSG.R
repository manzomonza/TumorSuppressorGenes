## Onco KB list

oncokb <- readr::read_tsv("cancerGeneList.tsv") |>
  janitor::clean_names() |>
  dplyr::relocate(is_oncogene, is_tumor_suppressor_gene)

tsg <- oncokb |>
  dplyr::filter(is_tumor_suppressor_gene == "Yes")

### ENSEMBLDB

library(ensembldb)
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
## Evaluate whether we have protein annotation available
hasProteinData(edb)

protLength <- function(genestring){
  prts <- proteins(edb, filter = GeneNameFilter(genestring),
                   return.type = "AAStringSet")
  whichMaxProtL <- width(prts)[which.max(width(prts))]
  if(length(whichMaxProtL)>0){
    #print(genestring)
    df = data.frame(gene = genestring, maxLength = whichMaxProtL )
    return(df)
  }
}
genes = ensembldb::genes(edb)
##### CHECK FOR H1 genes
hist1_genes = grep("^HIST1", genes$gene_name, value = TRUE)

check_symbols = tsg$hugo_symbol
# Create new vector with all genes to check lengths for
check_symbols = unique(c(check_symbols,hist1_genes, "BCL7A"))
# Create data.frame with gene name and max protein length (protLength)
tsg <- dplyr::bind_rows(lapply(check_symbols, protLength))

# Convert to list
tsg_ls <- as.list(tsg$maxLength)
names(tsg_ls) <- tsg$gene
saveRDS(tsg_ls, "OncoKB_TSG_maxLength.RDS")

# tsg_ls <- readRDS("OncoKB_TSG_maxLength.RDS")
# tsg_ls[["PTEN"]]

