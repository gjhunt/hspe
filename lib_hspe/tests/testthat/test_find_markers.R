.runThisTest <- Sys.getenv("RunAllHPSETests") == "yes"

if (.runThisTest) {
    
    library("hspe")
    
    Y <- matrix(c(1, 0, 1, 0, 0, 1, 0, 1), byrow = TRUE, nrow = 4)
    pure_samples <- list(1, 3)
    
    test_that("ratio markers are found as expected", {
        ml <- find_markers(Y, pure_samples = pure_samples)
        expect_equal(ml$L[[1]][[1]], 1)
        expect_equal(ml$L[[2]][[1]], 2)
        expect_equal(ml$V[[1]][[1]], 2)
        expect_equal(ml$V[[2]][[1]], 2)
    })
    
    test_that("diff markers are found as expected", {
        ml <- find_markers(Y, pure_samples = pure_samples, marker_method = "diff")
        expect_equal(ml$L[[1]][[1]], 1)
        expect_equal(ml$L[[2]][[1]], 2)
        expect_equal(ml$V[[1]][[1]], 1)
        expect_equal(ml$V[[2]][[1]], 1)
    })
    
    test_that("regression markers are found as expected", {
        ml <- find_markers(Y, pure_samples = pure_samples, marker_method = "regression")
        expect_equal(ml$L[[1]][[1]], 1)
        expect_equal(ml$L[[2]][[1]], 2)
        expect_equal(ml$V[[1]][[1]], 1)
        expect_equal(ml$V[[2]][[1]], 1)
    })
    
}
