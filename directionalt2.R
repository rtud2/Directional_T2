#' Two Sample Directional $T^2$ Tests
#'
#' Two Sample Directional $T^2$ Test. See Follmann (1996) and Glimm (2009) for more details.
#' @param x A data matrix 
#' @param y A second data matrix
#' @param direction Direction of Alternative Hypothesis. Default is "greater", which test $mu_x > mu_y$. Possible values: c("greater","less") 
#' @param method Directional $T^2$ method. Possible values: c("follmann", "glimm")
#' @return a list with 4 slots (additional slot for method = "glimm") including the mean difference, T^2 statistic, df1, df2, and p-value.
#' @export

directional_t2 <- function( x, y, method = "glimm", direction = "greater") {
    # check x and y are matrices
    if(ncol(x) != ncol(y)){
        stop(paste("dimension mismatch. ncol(x)=",ncol(x),", ncol(y)=",ncol(y)))
    }
    p <- ncol(x)

    # check directions
    if( direction == "greater" ){
        larger <- x
        smaller <- y
    } else if (direction == "less") {
        larger <- y
        smaller <- x
    } else {
        stop(paste0("direction = ", direction," is invalid"))
    }

    larger_mean = colMeans(larger)
    smaller_mean = colMeans(smaller)
    mean_diff = matrix(larger - smaller, nrow = 1)
    n1 <- nrow(larger)
    n2 <- nrow(smaller)

    larger_cov <- cov(larger)
    smaller_cov <- cov(smaller)
    pooled_cov <- ((n1 - 1) * larger_cov + (n2 - 1) * smaller_cov)/(n1 + n2 - 2)
    C <- (n1 + n2 - p - 1)/(p * (n1 + n2 -1)) * (n1 * n2)/(n1 + n2)

    if( method == "glimm" ){
        min_mu_obj <- function(mu, G, mean_diff, n1, n2, p, C){
        x_bar <- matrix(mu, nrow = 1)  - mean_diff 
        C * x_bar %*% solve(G) %*% t(x_bar)
        }

        ex_val <- optim(par = rep(0,p),
                        fn = min_mu_obj,
                        upper = rep(0,p),
                        method = "L-BFGS-B",
                        mean_diff = mean_diff,
                        G = pooled_cov,
                        n1 = n1, 
                        n2 = n2,
                        p = 3,
                        C = C)

        return(list(negative.orthant = ex_val$par,
                    mean.diff = mean_diff,
                    t.2 = ex_val$value,
                    df1 = p, 
                    df2 = n1 + n2 - p -1,
                    p.value = pf(ex_val$value, df1 = p, df2 =  n1 + n2 - p - 1, lower.tail = F)/2))

    } else if (method == "follmann"){
        my_stat <-  C * mean_diff %*% solve(pooled_cov) %*% t(mean_diff)
        return(list(mean.diff = mean_diff,
                    t.2 = my_stat,
                    df1 = p,
                    df2 =  n1 + n2 - p - 1,
                    p.value = pf(my_stat, df1 = p, df2 =  n1 + n2 - p - 1, lower.tail = F)/2))

    } else {
        stop(paste0("method = ", method," is invalid"))
    }
}
