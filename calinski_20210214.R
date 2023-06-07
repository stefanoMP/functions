calinski <- function (hhc, gMax = NULL) {
  
  dist <- cophenetic(hhc)
  attr(dist, "method") <- "cophenetic"
  
  dist <- as.matrix(dist)
  n <- nrow(dist)

  if (is.null(gMax)) gMax = round(1 + 3.3 * log(n, 10)) 
  
  dist <- dist^2
  A <- -dist/2
  A_bar <- apply(A, 1, mean)
  totalSum <- sum(diag(A) - 2 * A_bar + mean(A))
  ans <- rep(0, gMax)
  for (g in 2:gMax) {
    cclust <- cutree(hhc, k = g)
    withinSum <- 0
    for (k in 1:g) {
      if (sum(cclust == k) == 1) 
        next
      A <- as.matrix(-dist/2)[cclust == k, cclust == k]
      A_bar <- apply(A, 1, mean)
      withinSum <- withinSum + sum(diag(A) - 2 * A_bar + mean(A))
    }
    betweenSum <- totalSum - withinSum
    betweenSum <- betweenSum/(g - 1)
    withinSum <- withinSum/(n - g)
    ans[g] <- betweenSum/withinSum
  }
  
  class(ans) <- "calinski"
  attr(ans, "noof_items") <- n
#  attr(ans, "heights") <- hhc$height
  attr(ans, "distance") <- attr(dist, "method")
  
  return(ans)
}

plot.calinski <- function(obj, add = FALSE, from = 1, to = floor(0.75 * attr(obj, "noof_items")), 
                          height = 0.6, shift = 0.25, max_height = max(attr(obj, "heights"))) {
  
    if(!add) {
      plot(unclass(obj), type = "l", col = "grey", main = "Calinski & Harabasz curve", xlab = "number of groups", ylab = "CH scores")
      height  <-  1.0
      shift <- 0.0
    }
  
    G <- length(obj)
    nums <- paste(1:G); nums[1] <- ""
    xx <- 1:G
  
    if(add) obj <- obj/max(obj) * height * max_height
    
    ccol <- rep("black", G)
    for (g in 2:(G - 1)) {
      check <- obj[g - 1] < obj[g] & obj[g + 1] < obj[g]
      if(check) ccol[g] <- "red"
    }
    
    if(add) {
      xx <- floor(seq(from = from, to = to, length.out = G))
      obj <- obj + shift * max_height
      lines(xx, obj, col = "gray30", lty = "longdash")
    }
    
    text(xx, obj, nums, col = ccol)
    invisible(NULL)
  }

print.calinski <- function (obj) {
  message("No. of groups explored: ", G <- length(obj))
  maxima <- c()
  for (g in 2:(G - 1)) {
    check <- obj[g - 1] < obj[g] & obj[g + 1] < obj[g]
    if(check) maxima <- c(maxima, g)
  }
#  message("Suggested no. of groups: ", maxima) # 20210214
  message("Suggested no. of groups: ", paste(maxima, collapse = ", ")) # 20210315
}
  
#hc <- hclust(dist(USArrests), "ward.D2")
#aCalinski <- calinski(hc); plot(aCalinski)

#plot(hc, hang = -1) 
#plot(aCalinski, add = TRUE, height = 0.5, shift = 0.35)
