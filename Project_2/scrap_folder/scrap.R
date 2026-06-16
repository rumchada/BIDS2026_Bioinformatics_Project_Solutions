#codes that are completely unnecessary but still cool

library(TCGAbiolinks)

#stack TCGA
#file_ids <- c("abb04117-9caf-4882-a8fc-80450bdc9446",
#   "95a72aff-106e-4eb3-a897-1629a4b55c55",
#  "1f12acd4-1569-431c-b7b2-448638086ca4",
#  "a47f3dac-45bb-4f7d-ad28-a7ead099891e",
#  "babd0017-64a6-440f-8b82-e7a4de5cad8b")

#Check if those dataset exist within the TCGA-PRAD query
#query_results <- results[results$id, ]

#data_ids <- query_results$id
#data_ids

query <- GDCquery(
  project = "TCGA-PRAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
# Primary Tumor: Direct tumor samples
#Solid Tissue Normal: These patient samples are called "solid tissue normals" and are taken from normal tissues near the tumor
results <- getResults(query)
results
#Check if those dataset exist within the TCGA-PRAD query
data_ids <- results$cases

data_vector <- list()
meta_data <- list()


for(barcode in data_ids){
  
  query_uuid <- GDCquery(
    project = "TCGA-PRAD",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    barcode = barcode)
  
  
  GDCdownload(query_uuid)
  
  se <- GDCprepare(query_uuid)
  
  data_vector[[barcode]] <- se
  
}

#processing meta_data of each pull dataset
meta_data <- lapply(data_vector, function(x) {
  #convert into a datafram and transpose
  df <- t(as.data.frame(colData(x))) %>%
    as.data.frame() %>%
    #convert the rownames to columns
    tibble::rownames_to_column("variable") %>%
    as_tibble() %>%
    #the current meta-data is a list type
    mutate(across(everything(), as.character)) %>%
    #filter for the covaritates chosen
    filter(variable %in% covariates) %>%
    #across each row of the variable column
    #use subsitiute any numerical with "," with a "."
    mutate(across(-variable, ~ gsub(",", ".", .x))) %>%
    # if there is an NA data type convert that into a "none" character type
    mutate(across(-variable, ~ ifelse(is.na(.x), "none", .x)))
  
})
meta_data