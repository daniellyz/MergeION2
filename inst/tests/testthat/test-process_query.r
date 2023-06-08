# test query using equal string
test <- process_query(library2, query = c("SCANS=6"))

test_that("prcessQuery", {
  expect_equal(test$SELECTED$sp, library2$sp[which(library2$metadata$SCANS == 6)]
  )
})

# test query using regex
testRegex <- process_query(library2, query = c("COMPOUND@xxx@^Lo"))

test_that("prcessQuery", {
  expect_equal(testRegex$SELECTED$sp, library2$sp[grepl("^Lo",library2$metadata$COMPOUND)]
  )
})

