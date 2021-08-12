mxcov_a <- function(X, M, COV=NULL, d){                      # this function is for M~X estimation
  s_alpha <- matrix(0, 3, d)
  for(j in 1:d){
    if (is.null(COV)) {
      MX <- data.frame(M = M[, j], X = X)
    } else {
      MX <- data.frame(M = M[, j], X = X, COV = COV)
    }
    fit <- lm(M ~ ., data = MX)
    s_alpha[1,j] <- summary(fit)$coef[2,1]                   # coefficients of alpha
    s_alpha[2,j] <- (summary(fit)$coef[2,2])^2               # var of alpha
    s_alpha[3,j] <- summary(fit)$coef[2,4]                   # p-value of alpha
  }
  colnames(s_alpha) = colnames(M)[1:d]
  return(s_alpha=s_alpha)
}  

mainfun <- function(X, Y, M, COV = NULL, SISn = NULL, verbose = TRUE){
  
  n <- nrow(M); p <- ncol(M); censor = 1-mean(Y[,2])
  
  ### sis ###
  if(verbose) message("Step 1: Screening...", "     (", Sys.time(), ")")
  if (is.null(SISn)){SISn <- ceiling(2*n/log(n))}
  cor.mx <- c()
  if (is.null(COV)) {
    for (i in 1:p){cor.mx[i] <- abs(cor.test(M[,i], X)$estimate)}
  } else {
    for (i in 1:p){cor.mx[i] <- abs(pcor.test(M[,i], X, COV)$estimate)}
  }
  names(cor.mx) <- colnames(M)
  sort.mx <- sort(cor.mx, decreasing = T)
  subM1_ID <- names(head(sort.mx, SISn))
  subM1 <- M[,subM1_ID]
  M1_ID <- as.numeric(gsub("M","" ,subM1_ID))
  if(verbose) message("        ", length(subM1_ID), " mediators were selected in step 1")
  
  ### penalty ###
  if(verbose) message("Step 2: Regularization...", "     (", Sys.time(), ")")
  if (is.null(COV)){
    XCOV <- X
    MXCOV <- cbind(subM1, X)
  } else {
    XCOV <- cbind(X, COV)
    MXCOV <- cbind(subM1, X, COV)
  }
  b <- ncol(subM1); a <- ncol(XCOV)
  scadfun <- function(surv, X, weights){coef(ahaz(surv, X, univariate=TRUE))}
  pen <- ahazpen(Y, MXCOV, penalty = sscad.control(init.sol = scadfun), 
                 penalty.wgt = c(rep(1, b), rep(0, a)))
  cvpen <- tune.ahazpen(Y, MXCOV, penalty = sscad.control(init.sol = scadfun), 
                        tune = cv.control(rep = 5), 
                        penalty.wgt = c(rep(1,b),rep(0,a)))
  best.lambda <- which(pen$lambda == cvpen$lambda.min)
  betas <- pen$beta[1:b,]
  pen_ID <- which(betas[,best.lambda] != 0)
  subM2_ID <- colnames(subM1)[pen_ID]
  subM2 <- subM1[,subM2_ID]
  M2_ID <- as.numeric(gsub("M","" ,subM2_ID))
  
  if (length(M2_ID) == 0) {
    if(verbose) message ("Oops!None mediator is selected in step 2.")
    result00 <- data.frame("alpha_coef" = c(0), "alpha_var" = c(0), "alpha_pval" = c(0), 
                           "beta_coef" = c(0), "beta_var" = c(0), "beta_pval" = c(0),
                           "ab_true" = c(0), "ab_coef" = c(0), "ab_var" = c(0), "ab_pval" = c(1), 
                           "conf_low" = c(0), "conf_up" = c(0), 
                           "sob_pval_fdr" = c(1), "sob_pval_bon" = c(1), 
                           "join_pval_fdr" = c(0), "join_pval_bon" = c(0),
                           "sob_pval_by" = c(0), "join_pval_by" = c(0))
    rownames(result00) <- paste0("M",p+1)
    return(list(censor = censor, result = result00, subM1_ID = subM1_ID, subM2_ID = c()))
  }
  
  if(verbose) message("        mediator(s) ", paste(subM2_ID, collapse = " "), " is(are) selected in step 2")
  if(verbose) message("Step 3: Indirect effect test ...", "     (", Sys.time(), ")")
  ### significant test ###
  # true coef
  alpha_true <- alpha[M2_ID];  beta_true <- beta[M2_ID]
  ab_true <- alpha_true * beta_true
  names(ab_true) <- subM2_ID
  
  #beta estimated from ahaz
  if(is.null(COV)){
    subMXC <- as.matrix(cbind(subM2, X))
  } else {
    subMXC <- as.matrix(cbind(subM2, X, COV))
  }
  ahaz.fit <- ahaz(Y, subMXC)
  ahaz.coef <- summary(ahaz.fit)$coef
  beta_coef <- (ahaz.coef[,1])[1:length(M2_ID)]
  beta_var <- ((ahaz.coef[,2])^2)[1:length(M2_ID)]
  beta_pval <- (ahaz.coef[,4])[1:length(M2_ID)]
  
  # alpha estimated from linear regression
  alpha.fit <- mxcov_a(X, subM1, COV, length(M1_ID))                  
  alpha_coef <- (alpha.fit[1,subM2_ID])
  alpha_var <- (alpha.fit[2,subM2_ID])
  alpha_pval <- (alpha.fit[3,subM2_ID])
  
  # p value calculation
  ab_coef <- alpha_coef * beta_coef                                  # estimation of alpha * beta
  ab_var <- (alpha_coef^2) * (beta_var) + (beta_coef^2) * (alpha_var)
  
  # confidence interval
  conf_low <- ab_coef - 1.96 * sqrt(ab_var);  conf_up <- ab_coef + 1.96 * sqrt(ab_var)
  
  # sobel test for alpha and beta
  s.test <- (abs(ab_coef))/(sqrt(ab_var))                            # z-score for sobel test
  sob_pval <- 2 * (1-pnorm(s.test))                                  # p-value of sobel test
  sob_pval_fdr <- p.adjust(sob_pval, 'fdr', length(M2_ID))      
  sob_pval_fdr[sob_pval_fdr > 1] <- 1
  sob_pval_bon <- p.adjust(sob_pval, 'bonferroni', length(M2_ID))           
  sob_pval_bon[sob_pval_bon > 1] <- 1
  sob_pval_by <- p.adjust(sob_pval, 'BY', length(M2_ID))           
  sob_pval_by[sob_pval_by > 1] <- 1
  
  # joint test for alpha and beta
  alpha_pval_fdr <- p.adjust(alpha_pval, 'fdr', length(M2_ID))     
  alpha_pval_fdr[alpha_pval_fdr > 1] <- 1
  beta_pval_fdr <- p.adjust(beta_pval, 'fdr', length(M2_ID))
  beta_pval_fdr[beta_pval_fdr > 1] <- 1
  pval_bind_fdr <- rbind(alpha_pval_fdr, beta_pval_fdr)
  join_pval_fdr <- apply(pval_bind_fdr, 2, max)
  
  alpha_pval_bon <- p.adjust(alpha_pval, 'bonferroni', length(M2_ID))     
  alpha_pval_bon[alpha_pval_bon > 1] <- 1
  beta_pval_bon <- p.adjust(beta_pval, 'bonferroni', length(M2_ID))
  beta_pval_bon[beta_pval_bon > 1] <- 1
  pval_bind_bon <- rbind(alpha_pval_bon, beta_pval_bon)
  join_pval_bon <- apply(pval_bind_bon, 2, max)
  
  
  alpha_pval_by <- p.adjust(alpha_pval, 'BY', length(M2_ID))     
  alpha_pval_by[alpha_pval_by > 1] <- 1
  beta_pval_by <- p.adjust(beta_pval, 'BY', length(M2_ID))
  beta_pval_by[beta_pval_by > 1] <- 1
  pval_bind_by <- rbind(alpha_pval_by, beta_pval_by)
  join_pval_by <- apply(pval_bind_by, 2, max)
  
  sob_fdr <- sob_pval_fdr[which(sob_pval_fdr < 0.05)]
  
  if(verbose) message("Done!", "     (", Sys.time(), ")")
  if(verbose) message("        mediator(s) ", paste(names(sob_fdr), collapse = " "), " is(are) significant!")
  
  result <- data.frame(alpha_coef, alpha_var, alpha_pval, beta_coef, beta_var, beta_pval,
                       ab_true, ab_coef, ab_var, sob_pval, conf_low, conf_up, 
                       sob_pval_fdr, sob_pval_bon, join_pval_fdr, join_pval_bon, sob_pval_by, join_pval_by)
  
  return(list(censor = censor, result = result, subM1_ID = subM1_ID, 
              subM2_ID = subM2_ID, subM2 = subM2, sob_fdr = sob_fdr))
}