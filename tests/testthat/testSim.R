test_that("Simulation outputs are as expected", {
  ps <- 100
  out <- simpleSim(inc_mean_param = 6,
            inf_mean_param = 6,
            R_0 = 2.6,
            case_fatality = 0.02,
            pop_size = ps,
            inc_shape_param = 2,
            inf_shape_param = 2,
            seed = 12345)
  
  expect_lte(nrow(out$cases), ps)
  expect_gte(nrow(out$cases), 1)
  expect_equal(nrow(out$pop), ps)
  expect_true(is.data.frame(out$pop))
  expect_true(is.data.frame(out$cases))
  expect_true(all(c('id', 'doi', 'doo','dor', 'ill', 'died') %in% names(out$cases)))
  expect_true(all(out$cases$doi < out$cases$doo))
  expect_true(all(out$cases$doo < out$cases$dor))
  expect_false(any(out$cases$ill))
})