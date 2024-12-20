test_that("fixed results: calc_2bnab_BHtiter and calc_3bnab_BHtiter", {


  test_2mab_data = structure(list(
    a = c(400, 111, 11, 160, 582, 4245, 1335),
    b = c(0, 4610, 1139, 774, 1788, 5667, 120)
  ),
  row.names = c(NA, 6L),
  class = "data.frame")
  expect_equal(
    object = calc_2bnab_BHtiter(test_2mab_data$a, test_2mab_data$b),
    expected = c(100,1280.179, 298.0105, 327.9148, 886.0958, 3986.585, 452.2978),
    tolerance = 1e-2
  )

  test_3mab_data = structure(list(
    a = c(0, 333, 83, 191, 1200, 18, 47),
    b = c(0, 108, 15, 128, 167, 5534, 4636),
    c = c(400, 47, 405, 972, 3501, 28, 1024)
  ),
  row.names = c(NA, 6L),
  class = "data.frame")

  expect_equal(
    object = calc_3bnab_BHtiter(test_3mab_data$a, test_3mab_data$b, test_3mab_data$c),
    expected = c(100, 202.3578, 184.8139, 509.6948, 1915.56, 1439.631, 2051.053),
    tolerance = 1e-2
    )

})

test_that("calc_bnab_BHtiter when passed 1,2,3 active titers", {
  #pt80 is the default for these functions

  expect_equal(calc_3bnab_BHtiter(0, 0, 0), 0)
  expect_equal(calc_2bnab_BHtiter(0, 0), 0)
  expect_error(calc_3bnab_BHtiter(c(0,0)))


  expect_equal(calc_3bnab_BHtiter(1080*4, 0, 0), 1080)
  expect_equal(calc_3bnab_BHtiter(0, 1080*4, 0), 1080)
  expect_equal(calc_3bnab_BHtiter(0, 0, 1080*4), 1080)

  expect_equal(calc_3bnab_BHtiter(0, 1080, 900),
               calc_2bnab_BHtiter(1080, 900))
  expect_equal(calc_3bnab_BHtiter(1080, 900, 0),
               calc_2bnab_BHtiter(900, 1080))

})

test_that("benchmarking-check vectorization-2mab",{

  expect_error(calc_2bnab_BHtiter(c(16, 1080), c(1080)),
               "length(pt50_1) == length(pt50_2) is not TRUE",
               fixed = T)
  expect_error(calc_2bnab_BHtiter(c(16, 1080), c(1080,100,100)),
               "length(pt50_1) == length(pt50_2) is not TRUE",
               fixed = T)

  test_size = if(interactive()) 100 else 5
  set.seed(10)
  benchmark_data = data.frame(
    a = 10^runif(test_size, 1, 4),
    b = 10^runif(test_size, 1, 4)
  )

  rbenchmark::benchmark(
    implict_vector = { # this is currently fastest
      implicit = with(benchmark_data, calc_2bnab_BHtiter(a, b)
      )
    },
    explicit_vector_base = {
      explicit = mapply(calc_2bnab_BHtiter,
                        benchmark_data$a,
                        benchmark_data$b)
    },
    explicit_vector_purrr = {
      explicit_purrr = purrr::pmap_dbl(
        benchmark_data, function(a,b) calc_2bnab_BHtiter(a,b))
    }, replications = 10)


  expect_equal(implicit, explicit)
  expect_equal(implicit, explicit_purrr)

})


test_that("check vectorization-benchmarking-3mab", {

  # R recycles, so this first case is insidious without error checking
  expect_error(calc_3bnab_BHtiter(c(10, 0), c(16, 1080), c(1080)),
               "length(pt50_1) == length(pt50_2) & length(pt50_2) == length(pt50_3) is not TRUE",
               fixed = T)
  expect_error(calc_3bnab_BHtiter(c(10, 0), c(16, 1080), c(1080,100,100)),
               "length(pt50_1) == length(pt50_2) & length(pt50_2) == length(pt50_3) is not TRUE",
               fixed = T)

  expect_equal(
    object = calc_3bnab_BHtiter(c(10, 0), c(16, 1080), c(1080, 900)),
    expected = c(calc_3bnab_BHtiter(10, 16, 1080), calc_3bnab_BHtiter(0, 1080, 900))
  )

  test_size = if(interactive()) 100 else 5
  set.seed(10)
  benchmark_data = data.frame(
    a = 10^runif(test_size, 1, 4),
    b = 10^runif(test_size, 1, 4),
    c = 10^runif(test_size, 1, 4)
  )


  # implicit_vector           10   0.104    1.000     0.101    0.002          0         0
  rbenchmark::benchmark(
    implicit_vector = { # this is currently fastest
      implicit = with(benchmark_data,
                      calc_3bnab_BHtiter(a, b, c)
      )
      },
    explicit_vector_base = {
      explicit = mapply(calc_3bnab_BHtiter,
                        benchmark_data$a,
                        benchmark_data$b,
                        benchmark_data$c)
    },
    explicit_vector_purrr = {
      explicit_purrr = purrr::pmap_dbl(
        benchmark_data, function(a,b,c) calc_3bnab_BHtiter(a,b,c))
    }, replications = 10)


  expect_equal(implicit, explicit)
  expect_equal(implicit, explicit_purrr)


})

test_that(".calc_3bnab_BHPTxx error handling", {

  expected_check = calc_3bnab_BHtiter(10, 20, 300)
  expected_min = calc_3bnab_BHtiter(0.001, 0.001, 0.001)
  expected_na = calc_3bnab_BHtiter(NA, 0.001, 0.001)


  # basic checks
  expect_equal(.calc_3bnab_BHPTxx(data.frame(a = 10, b = 20, c = 300), .80, 0.01, T),
               expected_check)

  expect_equal(.calc_3bnab_BHPTxx(data.frame(a = 0, b = 0, c = 0), .80, 0.01, T),
               0.01)
  expect_equal(.calc_3bnab_BHPTxx(data.frame(a = 0, b = 0, c = 0), .80, 0.01, T),
               expected_min)

  # NA input
  expect_equal(.calc_3bnab_BHPTxx(data.frame(a = NA, b = 0, c = 0), .80, 0.01, T),
               NA_real_)
  expect_equal(.calc_3bnab_BHPTxx(data.frame(a = NA, b = 0, c = 0), .80, 0.01, F),
               NA_real_)

  expect_equal(calc_3bnab_BHtiter(NA, 0.001, 0.001), NA_real_)
  expect_equal(calc_3bnab_BHtiter(NA, 0.001, 0.001, NaN_as_NA = F), NA_real_)

  # NaN
  expect_equal(.calc_3bnab_BHPTxx(data.frame(a = NaN, b = 0, c = 0), .80, 0.01, T),
               NA_real_)
  expect_equal(calc_3bnab_BHtiter(NA, 0.001, 0.001), NA_real_)

  expect_error(.calc_3bnab_BHPTxx(data.frame(a = NaN, b = 0, c = 0), .80, 0.01, F),
               "BH PT50 NaN, input: NaN,0,0")
  expect_error(calc_3bnab_BHtiter(NaN, 0.001, 0.001, NaN_as_NA = F),
               "BH PT50 NaN, input: NaN,0.001,0.001")

  # other input
  expect_error(.calc_3bnab_BHPTxx(data.frame(a = "a", b = 0, c = 0), .80, 0.01, F),
               "error with uniroot solver (titer input?): non-numeric argument to binary operator",
               fixed = T)

  expect_error(calc_3bnab_BHtiter("a", 0.001, 0.001),
               "is.numeric(pt50_1) & is.numeric(pt50_2) & is.numeric(pt50_2) is not TRUE",
               fixed = T)

  expect_error(calc_3bnab_BHtiter(0.001, "b", 0.001),
               "is.numeric(pt50_1) & is.numeric(pt50_2) & is.numeric(pt50_2) is not TRUE",
               fixed = T)

  expect_error(calc_3bnab_BHtiter(0.001, 0.001, "c"),
               "is.numeric(pt50_1) & is.numeric(pt50_2) & is.numeric(pt50_2) is not TRUE",
               fixed = T)
})
