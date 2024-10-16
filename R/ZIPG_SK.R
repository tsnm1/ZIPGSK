#' Simultaneous knockoff for multi-source count data
#'
#' @param W Count data
#' @param class_K Matrix data, representing the input of count data;
#' @param data_x Vector data, representing different sources of multi-source data; default is NULL, representing one source or pooling method;
#' @param M The sequencing depth of count data, which was the row sum of matrix data;
#' @param y Response variable, representing binary, continuous or mixed variables;
#' @param T_var True correlation variable, which represents the set of variables that are truly correlated with the response variable, defaults to NULL;
#' @param fdr the target FDR level, default 0.2;
#' @param offset value between 0 and 1;
#' @param test_statistic Statistics for a single source dataset, c('DE','GLM','RF')
#' @param filter_statistics The statistics of the dataset were calculated from multiple sources. Denote c(cumprod,max, sum) by 1,2,3;
#' @param test1 Test method, select when selecting 'DE', default is "wilcox.test", c("wilcox.test",'ks.test','mmd','distance_JS').
#'
#' @return
#' c_w: Results of test statistics for multiple single-source datasets;
#' c_w_b: The filter statisticss of the dataset were calculated from multiple sources;
#' S: The variable selection set;
#' FDRPower: The resulting FDR and Power are calculated after given T_var.
#' @export
#'
#' @examples
#'
#' data(data_pg_copula)
#' i <- 2 # or 1
#' data_K_j <- data_pg_copula[[i]]
#' count_K_j <- data_K_j # [data_K_j[,1]==c,]
#' n_x1 <- 3
#' data_x <- as.data.frame(count_K_j[, c(1:n_x1 + 1)])
#' W <- as.data.frame(count_K_j[, -c(1:(n_x1 + 2))])
#' M <- count_K_j[, n_x1 + 2]
#' n_w <- dim(W)[1]
#' n_data <- 2
#' y <- rep(rep(c(1, 2), n_data), rep(n_w / 2 / n_data, n_data * 2))
#' class_K <- rep(c(1:n_data), rep(n_w / n_data, n_data))
#' n_p  <-  c(40,50)
#' n_p_all <- c(400,500)
#' T_var <- 1:n_p[i]
#' name_data <- names(table(class_K))
#' fdr <- 0.2
#'
#' ZIPG_DE_S3 <- ZIPG_SK(
#'   W = W, class_K = class_K, data_x = data_x, M = M, y = y, T_var = T_var,
#'   fdr = fdr, test_statistic = "DE", filter_statistics = 1
#' )
#' ZIPG_DE_S3$S
#' ZIPG_DE_S3$FDRPower
#'
ZIPG_SK <- function(W = W, class_K = class_K, data_x = NULL, M = NULL, y = y, T_var = 1, fdr = 0.2, offset = 1,
                    test_statistic = "DE", filter_statistics = 3, test1 = "wilcox.test") {
  if (!require(ZIPG)) devtools::install_github("roulan2000/ZIPG")
  # if (!require(scDesign2)) devtools::install_github("JSB-UCLA/scDesign2")
  if (!require(knockoff)) install.packages(knockoff)
  library(ZIPG)
  # library(scDesign2)
  library(knockoff)
  if (is.null(class_K)) {
    class_K <- rep(1, dim(W)[1])
  }
  if (is.null(data_x)) { # sum(data_x) == 1
    data_x <- as.data.frame(W[, c(1:3)])
    data_x[data_x != 0] <- 0
  }
  if (is.null(M)) {
    M <- apply(W, 1, sum)
  }
  name_data <- names(table(class_K))
  model_K <- list()
  c <- 1
  for (k in name_data) {
    # c = 1
    W_k <- W[class_K == k, ]
    data_x_k <- data_x[class_K == k, ]
    M_k <- M[class_K == k]

    W_k <- apply(W_k, 2, function(col) {
      col_replace <- mean(col, na.rm = T)
      col[is.na(col)] <- col_replace
      if (any(is.na(col))) {
        col[is.na(col)] <- 0
      }
      return(col)
    })

    copula_result <- ZIPG_Estimate_3_1(data_x_k, W_k, M_k)
    model_K[[c]] <- copula_result
    c <- c + 1
  }

  c_w <- c()
  for (k1 in 1:length(name_data)) {
    sub <- class_K == name_data[k1]
    W_k <- W[sub, ]
    data_x_k <- data_x[sub, ]
    M_k <- M[sub]
    y_k <- y[sub]
    copula_result <- model_K[[k1]]

    W_k <- apply(W_k, 2, function(col) {
      col_replace <- mean(col, na.rm = T)
      col[is.na(col)] <- col_replace
      if (any(is.infinite(col))) {
        col[is.infinite(col)] <- max(col[!is.infinite(col)])
      }
      if (any(is.na(col))) {
        col[is.na(col)] <- 0
      }
      return(col)
    })

    W_k_1 <- simulate_count_copula_3(copula_result, data_x_k, M_k)

    W_k_1 <- apply(W_k_1, 2, function(col) {
      col_replace <- mean(col, na.rm = T)
      col[is.na(col)] <- col_replace
      if (any(is.na(col))) {
        col[is.na(col)] <- 0
      }
      return(col)
    })

    if (test_statistic == "DE") {
      c_w_k <- contrast_score_computation(W_k, W_k_1, y_k, test1)
    } else if (test_statistic == "GLM") {
      test_statistic1 <- stat.glmnet_coefdiff
      random_m <- matrix(runif(dim(W_k)[1] * dim(W_k)[2], min = 0, max = 1), nrow = dim(W_k)[1])
      W_k <- W_k + random_m
      W_k_1 <- W_k_1 + random_m
      c_w_k <- test_statistic1(W_k, W_k_1, y_k)
    } else if (test_statistic == "RF") {
      test_statistic1 <- stat.random_forest
      c_w_k <- test_statistic1(W_k, W_k_1, y_k)
    }
    c_w <- rbind(c_w, c_w_k)
  }
  # c_w_b <- switch(filter_statistics,
  #                 apply(c_w, 2, cumprod)[dim(c_w)[1],],
  #                 apply(c_w, 2, max),
  #                 apply(c_w, 2, sum))
  if (k1 == 1) {
    c_w_b <- c_w
  } else {
    c_w_b <- switch(filter_statistics,
      apply(c_w, 2, cumprod)[dim(c_w)[1], ],
      apply(c_w, 2, max),
      apply(c_w, 2, sum)
    )
  }
  result_fdr <- clipper_BC_2(c_w_b, fdr)
  S <- result_fdr$discovery

  # FDRPower = NULL
  # if(!is.null(T_var)){
  #   FDRPower <- FDR_Power(S,T_var)
  # }
  FDRPower <- FDR_Power(S, T_var)
  FDRPower <- c(result_fdr$FDR, FDRPower$FP)

  return(list(fdr = fdr, c_w = c_w, c_w_b = c_w_b, S = S, FDRPower = FDRPower))
}


#' An Aggregating version of ZIPG_SK
#'
#' @param W Same as ZIPG_SK;
#' @param class_K Matrix data, representing the input of count data;
#' @param data_x Vector data, representing different sources of multi-source data; default is NULL, representing one source or pooling method;
#' @param M The sequencing depth of count data, which was the row sum of matrix data;
#' @param y Response variable, representing binary, continuous or mixed variables;
#' @param T_var True correlation variable, which represents the set of variables that are truly correlated with the response variable, defaults to NULL;
#' @param fdr the target FDR level, default 0.2;
#' @param offset value between 0 and 1;
#' @param B Number of aggregations;
#' @param Bstat For aggregation, 1 represents the aggregated test statistic and 2 represents the aggregated fliter statistic;
#' @param test_statistic Statistics for a single source dataset, c('DE','GLM','RF')
#' @param filter_statistics The statistics of the dataset were calculated from multiple sources. Denote c(cumprod,max, sum) by 1,2,3;
#' @param test1 Test method
#' @param combine_1 When combine_1 aggregation mode is 1, 'simul' means that the statistic is calculated based on the proposed method, and 'fisher' means that the statistic is calculated based on the combination of p values.
#'
#' @return
#' S: The variable selection set;
#' c_w_pis: Same as c_w;
#' contrast_score: fliter statistics for the non-aggregated case (B=1);
#' FDRPower: The resulting FDR and Power are calculated after given T_var.
#' chis: The threshold of the chi-square distribution using fisher's method in the first aggregation case (Bstat = 1, combine_1 == 'fisher');
#' p_comb:
#' e_value: Combined E-values for multi-source datasets;
#' e_w_B: The e value of the multi-source dataset in the aggregated case, B times;
#' e_value: E-values for single-source datasets.
#' @export
#'
ZIPG_SK_B <- function(W = W, class_K = NULL, data_x = NULL, M = NULL, y = y, T_var = NULL, fdr = 0.2, offset = 1,
                      B = 1, Bstat = 2, test_statistic = "DE", filter_statistics = 3, test1 = "wilcox.test",
                      combine_1 = "simul") {
  if (!require(ZIPG)) devtools::install_github("roulan2000/ZIPG")
  # if (!require(scDesign2)) devtools::install_github("JSB-UCLA/scDesign2")
  if (!require(knockoff)) install.packages(knockoff)
  library(ZIPG)
  # library(scDesign2)
  library(knockoff)
  if (is.null(class_K)) {
    class_K <- rep(1, dim(W)[1])
  }
  if (is.null(data_x)) { # sum(data_x) == 1
    data_x <- as.data.frame(W[, c(1:3)])
    data_x[data_x != 0] <- 0
  }
  if (is.null(M)) {
    M <- apply(W, 1, sum)
  }
  name_data <- names(table(class_K))
  model_K <- list()
  c <- 1
  for (k in name_data) {
    # c = 1
    W_k <- W[class_K == k, ]
    data_x_k <- data_x[class_K == k, ]
    M_k <- M[class_K == k]

    W_k <- apply(W_k, 2, function(col) {
      col_replace <- mean(col, na.rm = T)
      col[is.na(col)] <- col_replace
      if (any(is.na(col))) {
        col[is.na(col)] <- 0
      }
      return(col)
    })

    copula_result <- ZIPG_Estimate_3_1(data_x_k, W_k, M_k)
    model_K[[c]] <- copula_result
    c <- c + 1
  }

  c_w <- e_w <- e_w_Bk <- c()
  for (k1 in 1:length(name_data)) {
    sub <- class_K == name_data[k1]
    W_k <- W[sub, ]
    data_x_k <- data_x[sub, ]
    M_k <- M[sub]
    y_k <- y[sub]
    copula_result <- model_K[[k1]]

    W_k <- apply(W_k, 2, function(col) {
      col_replace <- mean(col, na.rm = T)
      col[is.na(col)] <- col_replace
      if (any(is.infinite(col))) {
        col[is.infinite(col)] <- max(col[!is.infinite(col)])
      }
      if (any(is.na(col))) {
        col[is.na(col)] <- 0
      }
      return(col)
    })

    e_w_B <- c_w_B <- c()
    for (b in 1:B) {
      W_k_1 <- simulate_count_copula_3(copula_result, data_x_k, M_k)

      W_k_1 <- apply(W_k_1, 2, function(col) {
        col_replace <- mean(col, na.rm = T)
        col[is.na(col)] <- col_replace
        if (any(is.na(col))) {
          col[is.na(col)] <- 0
        }
        return(col)
      })

      if (test_statistic == "DE") {
        c_w_k <- contrast_score_computation(W_k, W_k_1, y_k, test1)
      } else if (test_statistic == "GLM") {
        test_statistic1 <- stat.glmnet_coefdiff
        random_m <- matrix(runif(dim(W_k)[1] * dim(W_k)[2], min = 0, max = 1), nrow = dim(W_k)[1])
        W_k <- W_k + random_m
        W_k_1 <- W_k_1 + random_m
        c_w_k <- test_statistic1(W_k, W_k_1, y_k)
      } else if (test_statistic == "RF") {
        test_statistic1 <- stat.random_forest
        c_w_k <- test_statistic1(W_k, W_k_1, y_k)
      }
      c_w_B <- rbind(c_w_B, c_w_k)

      gamma <- fdr / 2
      offset <- 1
      e_w_b <- ekn(c_w_k, gamma, offset)
      e_w_B <- rbind(e_w_B, e_w_b)
    }
    e_w_Bk <- rbind(e_w_Bk, e_w_B)
    c_w <- rbind(c_w, c_w_B)

    e_w_B_k <- apply(e_w_B, 2, mean)
    e_w <- rbind(e_w, e_w_B_k)
  }

  if (B == 1) {
    c_w_k <- switch(filter_statistics,
      apply(c_w, 2, cumprod)[dim(c_w)[1], ],
      apply(c_w, 2, max),
      apply(c_w, 2, sum)
    )
    result_fdr <- clipper_BC_2(c_w_k, fdr)
    S <- result_fdr$discovery

    # FDRPower = NULL
    # if(!is.null(T_var)){
    #   FDRPower <- FDR_Power(S,T_var)
    # }
    FDRPower <- FDR_Power(S, T_var)
    FDRPower <- c(result_fdr$FDR, FDRPower$FP)
    return(list(
      fdr = fdr, B = B, c_w_pis = c_w, contrast_score = c_w_k,
      S = S, FDRPower = FDRPower
    ))
  } else if (Bstat == 1) {
    # B test_statistics, mean; First k and then B; Finally K row E-value;
    if (combine_1 == "simul") {
      e_value_B <- switch(filter_statistics,
        apply(e_w, 2, cumprod)[dim(e_w)[1], ],
        apply(e_w, 2, max),
        apply(e_w, 2, sum)
      )
      result_fdr <- clipper_BC_2(e_value_B, fdr)
      S <- result_fdr$discovery

      # FDRPower = NULL
      # if(!is.null(T_var)){
      #   FDRPower <- FDR_Power(S,T_var)
      # }
      FDRPower <- FDR_Power(S, T_var)
      FDRPower <- c(result_fdr$FDR, FDRPower$FP)
      return(list(
        fdr = fdr, B = B, c_w_pis = c_w, # contrast_score = c_w_B,
        e_value = e_w, e_value_B = e_value_B, combine_1 = combine_1,
        S = S, FDRPower = FDRPower
      ))
    } else if (combine_1 == "fisher") {
      p_e_w <- e_w
      p_e_w[p_e_w == 0] <- 1
      p_comb <- apply(-2 * log(1 / p_e_w), 2, sum)
      chis <- qchisq(fdr, 2 * length(name_data), lower.tail = FALSE)
      S <- which(p_comb > chis)

      # FDRPower = NULL
      # if(!is.null(T_var)){
      #   FDRPower <- FDR_Power(S,T_var)
      # }
      FDRPower <- FDR_Power(S, T_var)
      FDRPower <- c(fdr, FDRPower$FP)
      return(list(
        fdr = fdr, B = B, c_w_pis = c_w, chis = chis,
        e_value = e_w, p_comb = p_comb, combine_1 = combine_1,
        S = S, FDRPower = FDRPower
      ))
    }
  } else if (Bstat == 2) {
    # B filter_statistics, mean; First b and then K; Finally B row E-value;
    e_w_B <- c()
    for (b in 1:B) {
      c_w_bK <- c_w[seq(b, length(name_data) * B + b - 1, B), ] # t:11:((n_k*T)+t-1)
      c_w_b <- switch(filter_statistics,
        apply(c_w_bK, 2, cumprod)[dim(c_w_bK)[1], ],
        apply(c_w_bK, 2, max),
        apply(c_w_bK, 2, sum)
      )
      gamma <- fdr / 2
      offset <- 1
      e_w_b <- ekn(c_w_b, gamma, offset)
      e_w_B <- rbind(e_w_B, e_w_b)
    }
    e_value_B <- apply(e_w_B, 2, mean)
    result_ebh <- ebh_2(e_value_B, fdr)
    S <- result_ebh$discovery

    FDRPower <- NULL
    # if(!is.null(T_var)){
    #   FDRPower <- FDR_Power(S,T_var)
    # }
    FDRPower <- FDR_Power(S, T_var)
    FDRPower <- c(result_ebh$FDR, FDRPower$FP)
    return(list(
      fdr = fdr, B = B, c_w_pis = c_w, e_w_B = e_w_B,
      e_value = e_value_B, S = S, FDRPower = FDRPower
    ))
  }
}


#' ZIPG_SK method based on the cross-validation case
#'
#' @param W Same as ZIPG_SK;
#' @param class_K Matrix data, representing the input of count data;
#' @param data_x Vector data, representing different sources of multi-source data; default is NULL, representing one source or pooling method;
#' @param M The sequencing depth of count data, which was the row sum of matrix data;
#' @param y Response variable, representing binary, continuous or mixed variables;
#' @param T_var True correlation variable, which represents the set of variables that are truly correlated with the response variable, defaults to NULL;
#' @param fdr the target FDR level, default 0.2;
#' @param offset value between 0 and 1;
#' @param cv Fold number of cross validation;
#' @param method Synthetic data generation methods, c('ZIPG','ZINB','knockoff');
#' @param B Number of aggregations;
#' @param Bstat For aggregation, 1 represents the aggregated test statistic and 2 represents the aggregated fliter statistic;
#' @param test_statistic Statistics for a single source dataset, c('DE','GLM','RF')
#' @param filter_statistics The statistics of the dataset were calculated from multiple sources. Denote c(cumprod,max, sum) by 1,2,3;
#' @param test1 Test method
#' @param combine_1 When combine_1 aggregation mode is 1, 'simul' means that the statistic is calculated based on the proposed method, and 'fisher' means that the statistic is calculated based on the combination of p values.
#' @param In Boolean value that represents the result of the computation based on the Intersection.
#'
#' @return Same as ZIPG_SK and ZIPG_SK_B;
#' @export
#'
ZIPG_SK_cv <- function(W = W, class_K = NULL, data_x = NULL, M = NULL, y = y, T_var = NULL, fdr = 0.2, offset = 1,
                       cv = 1, method = "ZIPG",
                       B = 1, Bstat = 2, test_statistic = "DE", filter_statistics = 3, test1 = "wilcox.test",
                       combine_1 = "simul", In = F) {
  if (!require(ZIPG)) devtools::install_github("roulan2000/ZIPG")
  # if (!require(scDesign2)) devtools::install_github("JSB-UCLA/scDesign2")
  if (!require(knockoff)) install.packages(knockoff)
  library(ZIPG)
  library(scDesign2)
  library(knockoff)

  if (is.null(class_K)) {
    class_K <- rep(1, dim(W)[1])
  }
  if (is.null(data_x)) { # sum(data_x) == 1
    data_x <- as.data.frame(W[, c(1:3)])
    data_x[data_x != 0] <- 0
  }
  if (is.null(M)) {
    M <- apply(W, 1, sum)
  }
  n_data <- names(table(class_K))
  index <- c()
  for (K in names(table(class_K))) {
    index1 <- sample(which(class_K == K))
    index1 <- matrix(index1[1:(floor(length(index1) / cv) * cv)], nrow = cv)
    index <- cbind(index, index1)
  }
  index <- as.matrix(index)
  S_cv <- list()
  for (i in c(1:cv)) {
    W_f <- W[-index[i, ], ]
    class_K_f <- class_K[-index[i, ]]
    data_x_f <- data_x[-index[i, ], ]
    M_f <- M[-index[i, ]]
    y_f <- y[-index[i, ]]

    if (B == 1) {
      result_fold <- ZIPG_SK_other(
        W = W_f, class_K = class_K_f, data_x = data_x_f, M = M_f, y = y_f,
        fdr = fdr, offset = offset, method = method,
        test_statistic = test_statistic, filter_statistics = filter_statistics,
        test1 = test1
      )
      # return(list(fdr = fdr, c_w = c_w, c_w_b = c_w_b, S = S))
    } else {
      result_fold <- ZIPG_SK_B_other(
        W = W_f, class_K = class_K_f, data_x = data_x_f, M = M_f, y = y_f,
        fdr = fdr, offset = offset, method = method,
        B = B, Bstat = Bstat, test_statistic = test_statistic, filter_statistics = filter_statistics,
        test1 = test1, combine_1 = combine_1
      )
    }
    S_cv[[i]] <- result_fold
  }
  # S_cv
  f2 <- ifelse(filter_statistics == 1, 3, 1)
  r_ZIPG_SK_cv <- list(fdr = fdr, cv = cv, S_cv = S_cv, K = length(n_data))
  r_cv <- result_cv(r_ZIPG_SK_cv, B = B, filter_statistics = f2, In = In, b_1 = T_var)

  return(list(r_ZIPG_SK_cv = r_ZIPG_SK_cv, r_cv = r_cv))
}

#' Title
#'
#' @param r_ZIPG_SK_cv
#' @param B
#' @param filter_statistics
#' @param In
#' @param b_1
#'
#' @return
#' SS1: Results based on the statistic of fliter in the method in the CV case
#' FP: Same as FDRPower;
#' SS_f1: results based on other fliter's statistics in the CV case.
#' @export
#'
result_cv <- function(r_ZIPG_SK_cv, B = 1, filter_statistics = NULL, In = F, b_1 = NULL) {
  cv <- r_ZIPG_SK_cv$cv
  fdr <- r_ZIPG_SK_cv$fdr
  K <- r_ZIPG_SK_cv$K
  FP <- SS <- SS_f <- SS_I <- list()
  SS1 <- SS_f1 <- SS_I1 <- c()
  S_cv <- r_ZIPG_SK_cv$S_cv
  for (i in 1:cv) {
    result_fold <- S_cv[[i]]
    if (B == 1) {
      SS[[i]] <- result_fold$S
      SS1 <- c(SS1, result_fold$S)
      FP[[i]] <- result_fold$FDRPower
      if (!is.null(filter_statistics)) {
        S_f <- filter_cw(c_w = result_fold$c_w, filter_statistics = filter_statistics, fdr = fdr)
        SS_f1 <- c(SS_f1, S_f$S)
      }
      if (I) {
        S_I <- inter_cw(c_w = result_fold$c_w, b_1 = b_1, fdr = fdr)
        SS_I1 <- c(SS_I1, S_I$S)
      }
    } else {
      combine_1 <- result_fold$combine_1
      Bstat <- ifelse(is.null(combine_1), 2, 1)
      SS1 <- c(SS1, result_fold$S)
      FP[[i]] <- result_fold$FDRPower
      if (!is.null(filter_statistics)) {
        S_f <- filter_cw_B(
          Bstat = Bstat, combine_1 = combine_1, c_w = result_fold$c_w_pis, e_w = result_fold$e_value,
          filter_statistics = filter_statistics, B = result_fold$B, fdr = fdr
        )
        # SS_f[[i]] <- S_f$S;
        SS_f1 <- c(SS_f1, S_f$S)
      }
    }
  }
  if (B == 1) {
    return(list(SS1 = SS1, FP = FP, SS_f1 = SS_f1, SS_I1 = SS_I1))
  } else {
    return(list(SS1 = SS1, FP = FP, SS_f1 = SS_f1, Bstat = Bstat, combine_1 = combine_1))
  }
}


#' Extension of ZIPG_SK to consider more synthetic data methods.
#'
#' @param W Same as ZIPG_SK;
#' @param class_K
#' @param data_x
#' @param M
#' @param y
#' @param T_var
#' @param fdr
#' @param offset
#' @param method Synthetic data generation methods, c('ZIPG','ZINB','knockoff');
#' @param test_statistic
#' @param filter_statistics
#' @param test1
#'
#' @return Same as ZIPG_SK;
#' @export
#'
#' @examples
#'
#' data(data_pg_copula)
#' i <- 1 # or 2
#' data_K_j <- data_pg_copula[[i]]
#' count_K_j <- data_K_j # [data_K_j[,1]==c,]
#' n_x1 <- 3
#' data_x <- as.data.frame(count_K_j[, c(1:n_x1 + 1)])
#' W <- as.data.frame(count_K_j[, -c(1:(n_x1 + 2))])
#' M <- count_K_j[, n_x1 + 2]
#' n_w <- dim(W)[1]
#' n_data <- 2
#' y <- rep(rep(c(1, 2), n_data), rep(n_w / 2 / n_data, n_data * 2))
#' class_K <- rep(c(1:n_data), rep(n_w / n_data, n_data))
#' n_p  <-  c(40,50)
#' n_p_all <- c(400,500)
#' T_var <- 1:n_p[i]
#' name_data <- names(table(class_K))
#' fdr <- 0.2
#'
#' ZIPG_DE_S3 <- ZIPG_SK_other(
#'   W = W, class_K = class_K, data_x = data_x, M = M, y = y, T_var = T_var,
#'   fdr = fdr, method = "ZIPG", test_statistic = "DE", filter_statistics = 1
#' )
ZIPG_SK_other <- function(W = W, class_K = NULL, data_x = NULL, M = NULL, y = y, T_var = NULL, fdr = 0.2, offset = 1,
                          method = "ZIPG", test_statistic = "DE", filter_statistics = 3, test1 = "wilcox.test") {
  if (!require(ZIPG)) devtools::install_github("roulan2000/ZIPG")
  if (!require(scDesign2)) devtools::install_github("JSB-UCLA/scDesign2")
  if (!require(knockoff)) install.packages(knockoff)
  library(ZIPG)
  library(scDesign2)
  library(knockoff)

  if (is.null(class_K)) {
    class_K <- rep(1, dim(W)[1])
  }
  if (is.null(data_x)) { # sum(data_x) == 1
    data_x <- as.data.frame(W[, c(1:3)])
    data_x[data_x != 0] <- 0
  }
  if (is.null(M)) {
    M <- apply(W, 1, sum)
  }
  name_data <- names(table(class_K))
  if (method == "ZIPG") {
    model_K <- list()
    c <- 1
    for (k in name_data) {
      W_k <- W[class_K == k, ]
      data_x_k <- data_x[class_K == k, ]
      M_k <- M[class_K == k]

      W_k <- apply(W_k, 2, function(col) {
        col_replace <- mean(col, na.rm = T)
        col[is.na(col)] <- col_replace
        if (any(is.na(col))) {
          col[is.na(col)] <- 0
        }
        return(col)
      })

      copula_result <- ZIPG_Estimate_3_1(data_x_k, W_k, M_k)
      model_K[[c]] <- copula_result
      c <- c + 1
    }
  }

  c_w <- c()
  for (k1 in 1:length(name_data)) {
    sub <- class_K == name_data[k1]
    W_k <- W[sub, ]
    data_x_k <- data_x[sub, ]
    M_k <- M[sub]
    y_k <- y[sub]

    W_k <- apply(W_k, 2, function(col) {
      col_replace <- mean(col, na.rm = T)
      col[is.na(col)] <- col_replace
      if (any(is.infinite(col))) {
        col[is.infinite(col)] <- max(col[!is.infinite(col)])
      }
      if (any(is.na(col))) {
        col[is.na(col)] <- 0
      }
      return(col)
    })

    if (method == "ZIPG") {
      copula_result <- model_K[[k1]]
      W_k_1 <- simulate_count_copula_3(copula_result, data_x_k, M_k)
    } else if (method == "knockoff") {
      W_k_1 <- create.second_order(W_k)
    } else if (method == "ZINB") {
      W_k_1 <- scDesign2_simulation(W_k, y_k)
    }

    W_k_1 <- apply(W_k_1, 2, function(col) {
      col_replace <- mean(col, na.rm = T)
      col[is.na(col)] <- col_replace
      if (any(is.na(col))) {
        col[is.na(col)] <- 0
      }
      return(col)
    })

    if (test_statistic == "DE") {
      c_w_k <- contrast_score_computation(W_k, W_k_1, y_k, test1)
    } else if (test_statistic == "GLM") {
      # test_statistic1 <- stat.glmnet_coefdiff
      random_m <- matrix(runif(dim(W_k)[1] * dim(W_k)[2], min = 0, max = 1), nrow = dim(W_k)[1])
      W_k_1 <- W_k_1 + random_m
      W_k <- W_k + random_m
      c_w_k <- stat.glmnet_coefdiff(W_k, W_k_1, y_k)
    } else if (test_statistic == "RF") {
      # test_statistic1 <- stat.random_forest
      c_w_k <- stat.random_forest(W_k, W_k_1, y_k)
    }
    c_w <- rbind(c_w, c_w_k)
  }
  # print(c_w)
  if (k1 == 1) {
    c_w_b <- c_w
  } else {
    c_w_b <- switch(filter_statistics,
      apply(c_w, 2, cumprod)[dim(c_w)[1], ],
      apply(c_w, 2, max),
      apply(c_w, 2, sum)
    )
  }

  result_fdr <- clipper_BC_2(c_w_b, fdr)
  S <- result_fdr$discovery

  # FDRPower = NULL
  # if(!is.null(T_var)){
  #   FDRPower <- FDR_Power(S,T_var)
  # }
  FDRPower <- FDR_Power(S, T_var)
  FDRPower <- c(result_fdr$FDR, FDRPower$FP)
  return(list(fdr = fdr, c_w = c_w, c_w_b = c_w_b, S = S, FDRPower = FDRPower))
}

#' Extension of ZIPG_SK_B to consider more synthetic data methods.
#'
#' @param W Same as ZIPG_SK_B;
#' @param class_K
#' @param data_x
#' @param M
#' @param y
#' @param T_var
#' @param fdr
#' @param offset
#' @param method Synthetic data generation methods, c('ZIPG','ZINB','knockoff');
#' @param B
#' @param Bstat
#' @param test_statistic
#' @param filter_statistics
#' @param test1
#' @param combine_1
#'
#' @return Same as ZIPG_SK_B;
#' @export
#'
#' @examples
#'
#' data(data_pg_copula)
#' i <- 1 # or 2
#' data_K_j <- data_pg_copula[[i]]
#' count_K_j <- data_K_j # [data_K_j[,1]==c,]
#' n_x1 <- 3
#' data_x <- as.data.frame(count_K_j[, c(1:n_x1 + 1)])
#' W <- as.data.frame(count_K_j[, -c(1:(n_x1 + 2))])
#' M <- count_K_j[, n_x1 + 2]
#' n_w <- dim(W)[1]
#' n_data <- 2
#' y <- rep(rep(c(1, 2), n_data), rep(n_w / 2 / n_data, n_data * 2))
#' class_K <- rep(c(1:n_data), rep(n_w / n_data, n_data))
#' n_p  <-  c(40,50)
#' n_p_all <- c(400,500)
#' T_var <- 1:n_p[i]
#' name_data <- names(table(class_K))
#' fdr <- 0.2
#' B <- 10
#'
#' ZIPG_DE_S3_B2 <- ZIPG_SK_B_other(
#'   W = W, class_K = class_K, data_x = data_x, M = M, y = y, T_var = T_var, fdr = fdr, method = "ZIPG",
#'   B = B, Bstat = 1, test_statistic = "DE", filter_statistics = 1, combine_1 = "simul"
#' )
#' ZIPG_DE_S3$S
#' ZIPG_DE_S3$FDRPower
#'
ZIPG_SK_B_other <- function(W = W, class_K = NULL, data_x = NULL, M = NULL, y = y, T_var = NULL, fdr = 0.2, offset = 1,
                            method = "ZIPG", B = 1, Bstat = 2, test_statistic = "DE", filter_statistics = 3, test1 = "wilcox.test",
                            combine_1 = "simul") {
  if (!require(ZIPG)) devtools::install_github("roulan2000/ZIPG")
  if (!require(scDesign2)) devtools::install_github("JSB-UCLA/scDesign2")
  if (!require(knockoff)) install.packages("knockoff")
  library(ZIPG)
  library(scDesign2)
  library(knockoff)

  if (is.null(class_K)) {
    class_K <- rep(1, dim(W)[1])
  }
  if (is.null(data_x)) { # sum(data_x) == 1
    data_x <- as.data.frame(W[, c(1:3)])
    data_x[data_x != 0] <- 0
  }
  if (is.null(M)) {
    M <- apply(W, 1, sum)
  }
  name_data <- names(table(class_K))
  if (method == "ZIPG") {
    model_K <- list()
    c <- 1
    for (k in name_data) {
      # print(c)
      W_k <- W[class_K == k, ]
      data_x_k <- data_x[class_K == k, ]
      M_k <- M[class_K == k]

      W_k <- apply(W_k, 2, function(col) {
        col_replace <- mean(col, na.rm = T)
        col[is.na(col)] <- col_replace
        if (any(is.na(col))) {
          col[is.na(col)] <- 0
        }
        return(col)
      })

      copula_result <- ZIPG_Estimate_3_1(data_x_k, W_k, M_k)
      model_K[[c]] <- copula_result
      c <- c + 1
    }
  }

  c_w <- e_w <- e_w_Bk <- c()
  for (k1 in 1:length(name_data)) {
    sub <- class_K == name_data[k1]
    W_k <- W[sub, ]
    data_x_k <- data_x[sub, ]
    M_k <- M[sub]
    y_k <- y[sub]
    # copula_result = model_K[[k1]]

    W_k <- apply(W_k, 2, function(col) {
      col_replace <- mean(col, na.rm = T)
      col[is.na(col)] <- col_replace
      if (any(is.infinite(col))) {
        col[is.infinite(col)] <- max(col[!is.infinite(col)])
      }
      if (any(is.na(col))) {
        col[is.na(col)] <- 0
      }
      return(col)
    })

    e_w_B <- c_w_B <- c()
    for (b in 1:B) {
      # W_k_1 <- simulate_count_copula_3(copula_result,data_x_k,M_k)
      # cat('k1',k1,'b',b,'\n')
      if (method == "ZIPG") {
        copula_result <- model_K[[k1]]
        W_k_1 <- simulate_count_copula_3(copula_result, data_x_k, M_k)
      } else if (method == "knockoff") {
        W_k_1 <- create.second_order(W_k)
      } else if (method == "ZINB") {
        W_k_1 <- scDesign2_simulation(W_k, y_k)
      }

      W_k_1 <- apply(W_k_1, 2, function(col) {
        col_replace <- mean(col, na.rm = T)
        col[is.na(col)] <- col_replace
        if (any(is.na(col))) {
          col[is.na(col)] <- 0
        }
        return(col)
      })

      if (test_statistic == "DE") {
        c_w_k <- contrast_score_computation(W_k, W_k_1, y_k, test1)
      } else if (test_statistic == "GLM") {
        # test_statistic1 <- stat.glmnet_coefdiff
        random_m <- matrix(runif(dim(W_k)[1] * dim(W_k)[2], min = 0, max = 1), nrow = dim(W_k)[1])
        W_k_1 <- W_k_1 + random_m
        W_k <- W_k + random_m
        c_w_k <- stat.glmnet_coefdiff(W_k, W_k_1, y_k)
      } else if (test_statistic == "RF") {
        test_statistic1 <- stat.random_forest
        c_w_k <- test_statistic1(W_k, W_k_1, y_k)
      }
      c_w_B <- rbind(c_w_B, c_w_k)

      gamma <- fdr / 2
      offset <- 1
      e_w_b <- ekn(c_w_k, gamma, offset)
      e_w_B <- rbind(e_w_B, e_w_b)
    }
    e_w_Bk <- rbind(e_w_Bk, e_w_B)
    c_w <- rbind(c_w, c_w_B)

    e_w_B_k <- apply(e_w_B, 2, mean)
    e_w <- rbind(e_w, e_w_B_k)
  }

  if (B == 1) {
    # c_w_k <- switch(filter_statistics,
    #                 apply(c_w, 2, cumprod)[dim(c_w)[1],],
    #                 apply(c_w, 2, max),
    #                 apply(c_w, 2, sum))
    if (k1 == 1) {
      c_w_k <- c_w
    } else {
      c_w_k <- switch(filter_statistics,
        apply(c_w, 2, cumprod)[dim(c_w)[1], ],
        apply(c_w, 2, max),
        apply(c_w, 2, sum)
      )
    }
    result_fdr <- clipper_BC_2(c_w_k, fdr)
    S <- result_fdr$discovery

    # FDRPower = NULL
    # if(!is.null(T_var)){
    #   FDRPower <- FDR_Power(S,T_var)
    # }
    FDRPower <- FDR_Power(S, T_var)
    FDRPower <- c(result_fdr$FDR, FDRPower$FP)
    return(list(
      fdr = fdr, B = B, c_w_pis = c_w, contrast_score = c_w_k,
      S = S, FDRPower = FDRPower
    ))
  } else if (Bstat == 1) {
    # B test_statistics, mean; First k and then B; Finally K row E-value;
    if (combine_1 == "simul") {
      e_value_B <- switch(filter_statistics,
        apply(e_w, 2, cumprod)[dim(e_w)[1], ],
        apply(e_w, 2, max),
        apply(e_w, 2, sum)
      )
      result_fdr <- clipper_BC_2(e_value_B, fdr)
      S <- result_fdr$discovery

      # FDRPower = NULL
      # if(!is.null(T_var)){
      #   FDRPower <- FDR_Power(S,T_var)
      # }
      FDRPower <- FDR_Power(S, T_var)
      FDRPower <- c(result_fdr$FDR, FDRPower$FP)
      return(list(
        fdr = fdr, B = B, c_w_pis = c_w, # contrast_score = c_w_B,
        e_value = e_w, e_value_B = e_value_B, combine_1 = combine_1,
        S = S, FDRPower = FDRPower
      ))
    } else if (combine_1 == "fisher") {
      p_e_w <- e_w
      p_e_w[p_e_w == 0] <- 1
      p_comb <- apply(-2 * log(1 / p_e_w), 2, sum)
      chis <- qchisq(fdr, 2 * length(name_data), lower.tail = FALSE)
      S <- which(p_comb > chis)

      # FDRPower = NULL
      # if(!is.null(T_var)){
      #   FDRPower <- FDR_Power(S,T_var)
      # }
      FDRPower <- FDR_Power(S, T_var)
      FDRPower <- c(fdr, FDRPower$FP)
      return(list(
        fdr = fdr, B = B, c_w_pis = c_w, chis = chis,
        e_value = e_w, p_comb = p_comb, combine_1 = combine_1,
        S = S, FDRPower = FDRPower
      ))
    }
  } else if (Bstat == 2) {
    # B filter_statistics, mean; First b and then K; Finally B row E-value;
    e_w_B <- c()
    for (b in 1:B) {
      c_w_bK <- c_w[seq(b, length(name_data) * B + b - 1, B), ] # t:11:((n_k*T)+t-1)
      c_w_b <- switch(filter_statistics,
        apply(c_w_bK, 2, cumprod)[dim(c_w_bK)[1], ],
        apply(c_w_bK, 2, max),
        apply(c_w_bK, 2, sum)
      )
      gamma <- fdr / 2
      offset <- 1
      e_w_b <- ekn(c_w_b, gamma, offset)
      e_w_B <- rbind(e_w_B, e_w_b)
    }
    e_value_B <- apply(e_w_B, 2, mean)
    result_ebh <- ebh_2(e_value_B, fdr)
    S <- result_ebh$discovery
    # FDRPower = NULL
    # if(!is.null(T_var)){
    #   FDRPower <- FDR_Power(S,T_var)
    # }
    FDRPower <- FDR_Power(S, T_var)
    FDRPower <- c(result_ebh$FDR, FDRPower$FP)
    return(list(
      fdr = fdr, B = B, c_w_pis = c_w, e_w_B = e_w_B, e_value = e_value_B,
      S = S, FDRPower = FDRPower
    ))
  }
}


# devtools::document()
