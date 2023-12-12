
test <- dataForPlotly(library2, id = 2) 

test_that("plotlyData", {
  expect_equal(as.character(unique(test$ID)), "2")
})

test_that("plotlyData", {
  expect_is(test$ID, "factor")
})


test_that("network visulisation", {
  
  skip("development")
  
  
})