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


  rbenchmark::benchmark(
    implict_vector = { # this is currently fastest
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
