.runThisTest <- Sys.getenv("RunAllHPSETests") == "yes"

if (.runThisTest) {
    
    
    ## Basic
    truth <- shen_orr_ex$annotation$mixture
    pure_samples <- lapply(1:3, function(i) {
        which(truth[, i] == 1)
    })
    Y <- shen_orr_ex$data$log
    n_markers <- 20
    
    # Basic: Y and References and pure_samples
    test_that("basic hspe Y and refs and ps", {
        Y <- shen_orr_ex$data$log
        n_markers <- 20
        references <- Y[unlist(pure_samples), ]
        Y <- Y[-unlist(pure_samples), ]
        dt_out <- hspe(Y, references = references, pure_samples = pure_samples, n_markers = n_markers, 
            marker_method = "ratio", seed = 4261992)
        expect_equal_to_reference(dt_out, file = "basic_hspe.rds")
    })
    
    # Basic: Y and References
    test_that("basic hspe Y and refs", {
        Y <- shen_orr_ex$data$log
        n_markers <- 20
        comp <- readRDS("basic_hspe.rds")
        references <- Y[unlist(pure_samples), ]
        Y <- Y[-unlist(pure_samples), ]
        combined_refs <- t(sapply(pure_samples, function(x) colMeans(references[x, 
            ])))
        dt_out <- hspe(Y, references = combined_refs, n_markers = n_markers, markers = comp$markers, 
            seed = 4261992)
        expect_equal(dt_out[-5], comp[-5], tolerance = 1e-05)
    })
    
    # Basic: Y and pure_samples
    test_that("basic hspe Y and pure_samples", {
        Y <- shen_orr_ex$data$log
        n_markers <- 20
        comp <- readRDS("basic_hspe.rds")
        dt_out <- hspe(Y, pure_samples = pure_samples, n_markers = n_markers, markers = comp$markers, 
            seed = 4261992)
        expect_equal(dt_out[-5], comp[-5], tolerance = 1e-05)
        # expect_equal(dt_out$estimates[-unlist(pure_samples),],comp$estimates,tolerance=1E-5)
    })
    
    # Basic: n_markers
    test_that("basic hspe n_markers", {
        Y <- shen_orr_ex$data$log
        n_markers <- 20
        comp <- readRDS("basic_hspe.rds")
        references <- Y[unlist(pure_samples), ]
        Y <- Y[-unlist(pure_samples), ]
        combined_refs <- t(sapply(pure_samples, function(x) colMeans(references[x, 
            ])))
        dt_out <- hspe(Y, references = combined_refs, n_markers = NULL, markers = comp$markers, 
            seed = 4261992)
        expect_equal(dt_out, comp)
        dt_out <- hspe(Y, references = combined_refs, n_markers = c(10, 11, 12), 
            markers = comp$markers, seed = 4261992)
        expect_equal_to_reference(dt_out$estimates, "basic_hspe_markers.rds")
    })
    
    # Basic: markers
    test_that("basic hspe markers", {
        Y <- shen_orr_ex$data$log
        n_markers <- 20
        dt_out <- hspe(Y, pure_samples = pure_samples, n_markers = c(10, 11, 12), 
            marker_method = "regression", seed = 4261992)
        expect_equal_to_reference(dt_out$estimates, "basic_hspe_marker_reg.rds")
    })
    
}

