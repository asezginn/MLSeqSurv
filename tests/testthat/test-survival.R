


# print(getwd())

test_that("NULL data doesn't work", {
  expect_error(MLSeqSurv::survival(NULL, "blackboost"), "MLSeqSurv object is null")
})

# test_that("NULL method doesn't work", {
#   expect_error(MLSeqSurv::survival(data_obj, NULL), "Method is not specified.")
# })
#
# test_that("Unavailable method doesn't work", {
#   expect_error(MLSeqSurv::survival(data_obj, "xd"), "Specified method is not available.")
# })

# test all methods
