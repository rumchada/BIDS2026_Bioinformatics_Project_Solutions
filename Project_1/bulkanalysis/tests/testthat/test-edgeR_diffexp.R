#source("../../R/edgeR_diffexp.R", chdir = TRUE)

library(bulkanalysis)
library(testthat)
library(RNAseqQC)

fake_ds1 <- readRDS("../test_data/fake_dds.rds")


dds_object <- fake_ds1
#Add the condition column name
condition_col = "condition"
#set a refreence group name
ref_group_name = "CONTROL"


#Performs inital differential expression on your DDS Object's raw counts dataframe
# The use give the condition column from the metadata to make the test parameters
# and the control group


#Input esting
test_that(
  "pvalue_input validations",{

    invalid_pvalues <- c(NA, "0.05")

    for (i in invalid_pvalues){
      # pvalue check
      expect_error(
        edgeR_diffexp(
          dds_object = dds_object,
          condition_col = condition_col,
          ref_group_name = ref_group_name,
          init_pval_cutoff = i
        )
      )}#end of loop

  })#end of test_that function

test_that(
  #condition column test
  "expect error for failed column name", {
    expect_error(
      edgeR_diffexp(
        dds_object = dds_object,
        condition_col = "fake_column",
        ref_group_name = ref_group_name,
        init_pval_cutoff = 0.5,
        adjust.method = "bonferroni"
      )
    )# end of internal test
  })


#expecting invalid p.adjust method
test_that(
  "valid p.adjust methods",{
    valid_p.adjust_method <- c("bonferroni","BH")
    for(i in valid_p.adjust_method)
      expect_no_error(
        edgeR_diffexp(
          dds_object = dds_object,
          condition_col = condition_col,
          ref_group_name = ref_group_name,
          init_pval_cutoff = 0.5,
          adjust.method = i
        )
      )
  })


