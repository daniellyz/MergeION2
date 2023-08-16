ref <- library2$metadata
ref2 <- ref
set.seed(34523)
ref2$ADDUCT <- sample(c("M+H","M+2H","M+Na","M+K","M+NH4", "M+"), length(ref$ADDUCT), replace = TRUE)


test <- process_metadata(ref=ref2, processing.algorithm = "Default",  polarity = c("Positive"), adductType = "M+2H")


test_that("subAdductType", {
  expect_equal(unique(test$ADDUCT), "M+2H")
})
