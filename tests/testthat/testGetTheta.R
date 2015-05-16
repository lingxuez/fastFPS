test_that("getTheta() works in p=30, ndim=1:5",{
  for (ndim in c(1:5)){
  v.test = sort(runif(30, min=0, max=10), decreasing=TRUE)
  theta = getTheta(v.test, ndim)
  v.test.p = pmin(pmax(v.test-theta$theta, 0), 1)
  expect_equal(sum(v.test.p), ndim)
  expect_equal(v.test.p[theta$last_act_i+1], 0)
  }
})
