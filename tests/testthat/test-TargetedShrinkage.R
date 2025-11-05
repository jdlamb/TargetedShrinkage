test_that("target works",{
  Sigma <- solve(matrix(ncol=2,c(4,2,2,3)))
  means <- c(1,2)
  expect_equal(target(means,Sigma),16/11)
})

test_that("JorionShrinkage works", {
  expect_equal(JorionShrinkage(1:3,diag(rep(1,3)),20),
               c(1.11111111,2,2.88888889))
  expect_equal(mean(JorionShrinkage(1:3,diag(rep(1,3)),20)),mean(1:3))
})

test_that("JamesSteinShrinkage works", {
  expect_equal(JamesSteinShrinkage(1:3,diag(rep(1,3))),c(1.5,2,2.5))
  expect_equal(mean(JamesSteinShrinkage(1:3,diag(rep(1,3)))),mean(1:3))
})

test_that("targetedShrinkage works", {
  expect_equal(c(1.15,2.05,2.95,3.85),
               targetedShrinkage(1:4,diag(rep(0.5,4))))
  expect_equal(mean(targetedShrinkage(1:4,diag(rep(1,4)))),mean(1:4))
  expect_equal(target(1:4,diag(c(1,1,1,4))),7/3.25)
  expect_equal(c(1.39473684210526,2.05263157894737,2.71052631578947,3.36842105263158),
               targetedShrinkage(1:4,diag(c(1,1,1,4)),debug=FALSE))
  expect_equal(c(0.3,0.1,-0.1,-0.3)+targetedShrinkage(1:4,diag(c(1,1,1,1))),
               JamesSteinShrinkage(1:4,diag(c(1,1,1,1))))
 V <- matrix(ncol=4,c(1,0,0,1,0,1,0,0,0,0,1,0,1,0,0,4))
 x <- 1:4
 expect_equal(targetedShrinkage(x,V,debug=FALSE),c(1.2,2,2.8,3.6))
 # test translation invariance
 expect_equal(targetedShrinkage(x+pi,V,debug=FALSE),c(1.2,2,2.8,3.6)+pi)
})

test_that("quadraticLoss works", {
  x <- -1:1
  y <- 2*x
  expect_equal(quadraticLoss(y,x),2)
  M <- matrix(ncol=2,c(1,0,0,0,1,0))
  theta <- c(0,0)
  observed <- quadraticLoss(M,theta)
  expected = 2/3
  testthat::expect_equal(observed,expected)
  # Another test
  x <- 1:3
  y <- x+1
  Sigma <- diag(rep(1,3))
  expected <- 3
  observed <- quadraticLoss(x,y)
  testthat::expect_equal(observed,expected)
  observed <- quadraticLoss(x,y,Sigma)
  testthat::expect_equal(observed,expected)
  observed <- quadraticLoss(2*x+3,2*y+3,2*Sigma)
  testthat::expect_equal(observed,2*expected)
  # Test case with unequal variances
  Sigma <- diag(c(1,1/2,1/3))
  observed <- quadraticLoss(x,y,Sigma)
  expected <- 6 # 1+2+3
  testthat::expect_equal(observed,expected)
  observed <- quadraticLoss(2*x+3,2*y+3,2*Sigma)
  # Test case with Sigma singular
  Sigma <- diag(rep(1,4))
  Sigma[4,1:3] <- 1/sqrt(3)
  Sigma[1:3,4] <- 1/sqrt(3)
  testthat::expect_equal(det(Sigma),0) # singular
  x <- 1:4
  y <- x+1
  observed <- quadraticLoss(x,y,Sigma)
  expected <- 1+sqrt(3/4)
  testthat::expect_equal(observed,expected)
  observed <- quadraticLoss(2*x+3,2*y+3,2*Sigma)
})

test_that("whiteningShrinkage works",{
  Sigma <- diag(c(2,2,2))
  means <- 1:3
  expect_equal(targetedShrinkage(means,Sigma,scale=1),whiteningShrinkage(means,Sigma))
  expect_equal(targetedShrinkage(means,Sigma,scale=1),JamesSteinShrinkage(means,Sigma))
  expect_equal(targetedShrinkage(means+pi,Sigma,scale=1),whiteningShrinkage(means+pi,Sigma))
})
