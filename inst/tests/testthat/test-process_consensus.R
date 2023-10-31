libraryTest <- list(complete = library2, consensus = NULL, network = NULL)

test_that("consensus library", {

  aa <- capture.output(newLibrary <- process_consensus(input_library = libraryTest),   type = "message")
  expect_equivalent(aa, "Generating consensus library")
  expect_gte(nrow(newLibrary$consensus$metadata), 1)
  
  bb <- capture.output(
    newLibrary2 <- process_consensus(input_library = libraryTest, IDsUpdated = c("4", "6", "8", "12", "16")
                                     ),   type = "message")
  
  
  expect_equivalent(bb, "Missing consensus libary in input_library, generating consensus library for all compound present")
  expect_equal(newLibrary$consensus, newLibrary2$consensus)
  
  
  cc <- capture.output(
    newLibrary3 <- process_consensus(input_library = newLibrary2, IDsUpdated = c("4", "6", "8", "12", "16")
    ),   type = "message")
  
  expect_equivalent(cc, "Updating existing concensus library")
  
  expect_equal(newLibrary3$consensus$metadata, 
               newLibrary2$consensus$metadata[match(newLibrary3$consensus$metadata$SCANS, newLibrary2$consensus$metadata$SCANS),])
  
  
})