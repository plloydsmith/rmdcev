context("Test wtp")
#Sys.setenv("R_TESTS" = "")
#data(eggs, package = "rmdcev")

test_that("basic test", {
    # If the number below needs to be increased due to additional outputs,
    # ensure that the output size does not get too big when there are multiple
    # classes and many iterations.
    expect_true(1 < 1245000)
    # Add the following option to expect_true to print out the size in Travis
    # info = print(as.numeric(object.size(result))))
})
