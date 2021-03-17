wilcox_test <- function (x, y, alternative = c("greater", "two.sided", "less"), verb = FALSE) 
{
    alternative <- match.arg(alternative)
    notIsNA_x <- !is.na(x)
    notIsNA_y <- !is.na(y)
    x <- x[notIsNA_x]
    y <- y[notIsNA_y]
    if(((sum(notIsNA_x) != length(notIsNA_x)) | sum(notIsNA_y) != length(notIsNA_y)) & verb) 
        message("Some missing values in data.")
    
    nx <- length(x)
    ny <- length(y)
    tmp <- c(x, y)
    names(tmp) <- c(rep("x", nx), rep("y", ny))
    tmp <- rank(tmp)
    Tstat <- sum(tmp[which(names(tmp) == "y")])
    Ustat <- nx * ny + ny * (ny + 1)/2 - Tstat
    mu <- nx * ny/2
    sigma <- sqrt(mu * (nx + ny + 1)/6)
    zValue <- Ustat - mu
    correction <- switch(alternative, two.sided = sign(zValue) * 
        0.5, greater = 0.5, less = -0.5)
    zValue <- (zValue - correction)/sigma
    pValue <- switch(alternative, less = 1 - pnorm(zValue), greater = pnorm(-zValue), 
        two.sided = 2 * pnorm(-abs(zValue)))
    ans <- list(statistic = Ustat, p.value = pValue, alternative = alternative, 
        pod = Ustat/nx/ny)
    return(ans)
}
