# This function takes the x and y matrices in a regression and runs a penalized elastic net regression,
# using cross-validation to choose the hyperparameters of the regression
# ------- @ INPUTS
# xMat: a numeric matrix of independent variables, where each row corresponds to a point in time
# yMat: a 1-column numeric matrix of dependent variables, where each row corresponds to a point in time
# bigK: an integer representing the number of groups to segment the data in for cross-validation
# bigJ: an integer representing the number of times to redo the group segmentation
# alphaSearchGrid: a vector representing the alphas searched over - use 1 for LASSO or 0 for ridge regression
# verbose: a boolean to indicate whether to print progress while running
# ------- @ OUPUTS
# cvTable: a tibble of combinations of lambda, alphas attempted during hyperparameters seslection,
# and their resulting mean cross-validated errors  
# cvChart: a chart showing the cvTable results
# elasticNetReg: a glmnet result obejct
# coef: a tibble showing the shrunk and regularized coefficient estimates
# ------- @ DEPENDENCIES
# glmnet
# tidyverse
# foreach
# doParallel

getElasticNet = function(xMat, yMat, bigK = 10, bigJ = 100, alphaSearchGrid = seq(0, 1, 0.1), verbose = FALSE) {
  
  if (nrow(xMat) != nrow(yMat)) stop('xMat and yMat dimensions not equal')
  if (ncol(yMat) > 1) stop('yMat can only have one column')
  
  kFolds = lapply(1:bigJ, function(k) sample(rep(seq(bigK), length = nrow(xMat))))
  
  cvTable0 =
    lapply(1:length(kFolds), function(.fold) {
      
      if (verbose == TRUE) message(.fold)
      
      # Segment data into 5 random folds for cross-validation testing
      # Folds should be held constant across alphas; see page 5 of glmnet docs
      kFolds = lapply(1:bigJ, function(k) sample(rep(seq(bigK), length = nrow(xMat))))
      
      lapply(alphaSearchGrid, function(.alpha) {
        # Iterate over (alpha, lambda, foldgroup) combinations
        cv =
          glmnet::cv.glmnet(
            xMat,
            yMat,
            family = 'gaussian',
            alpha = .alpha,
            foldid = kFolds[[.fold]]
          )
        
        tibble(
          foldgroup = .fold,
          alpha = .alpha,
          lambda = cv$lambda,
          cvError := cv$cvm
          #lambda1Se = (lambda == cv$lambda.1se),
          #lambda1Min = (lambda == cv$lambda.min)
        ) %>% 
          return(.)
      }) %>% dplyr::bind_rows(.)
    }) %>% dplyr::bind_rows(.)
  
  
  # Averages across foldgroups
  cvTable =
    dplyr::distinct(cvTable0[, c('alpha', 'lambda')]) %>%
    purrr::transpose(.) %>%
    lapply(., function(pair)
      dplyr::filter(cvTable0, alpha == pair$alpha & lambda == pair$lambda) %>%
        {tibble(
          alpha = pair$alpha,
          lambda = pair$lambda,
          meanCvError = mean(.$cvError),
          foldGroupsCount = nrow(.)
        )}
    ) %>%
    dplyr::bind_rows(.) %>%
    dplyr::filter(., foldGroupsCount >= bigJ/2) %>%
    # Get min at each alpha
    dplyr::group_by(., alpha) %>%
    dplyr::mutate(., isMinGivenAlpha = (meanCvError == min(meanCvError))) %>%
    dplyr::ungroup(.) %>%
    # Get (alpha, lambda) which minimizes
    dplyr::mutate(., isMin = (meanCvError == min(meanCvError))) #%>%
  
  cvChartTable = dplyr::mutate(cvTable, logLambda = log(lambda))
  
  cvChart =
    ggplot(cvChartTable) +
    geom_line(aes_string(x = 'logLambda', y = 'meanCvError', color = 'alpha', group = 'alpha')) +
    geom_point(data = dplyr::filter(cvChartTable, isMinGivenAlpha == TRUE), color = 'green',
               aes_string(x = 'logLambda', y = 'meanCvError')) +
    geom_point(data = dplyr::filter(cvChartTable, isMin == TRUE), color = 'red', shape = 4, size = 3,
               aes_string(x = 'logLambda', y = 'meanCvError')) +
    labs(
      title = 'Elastic Net Mean Cross-Validation Error By Lambda at Each Alpha',
      x = 'log(Lambda)',
      y = 'Mean CV Error',
      caption = 'Green: Minimum (Alpha, Lambda) at each Alpha, Red: Minimum (Alpha, Lambda) Overall'
    )
  
  optim = dplyr::filter(cvTable, isMin == TRUE)
  
  message(paste0('Optim alpha = ',optim$alpha,', lambda = ',optim$lambda))
  
  elasticNetReg =
    glmnet::glmnet(
      x = xMat,
      y = yMat,
      lambda = optim$lambda,
      alpha = optim$alpha,
      family = 'gaussian'
    )
  
  coef =
    glmnet::coef.glmnet(elasticNetReg) %>%
    as.matrix(.) %>%
    as.data.frame(.) %>%
    rownames_to_column(., 'coefname') %>%
    {
      if (!is.null(colnames(yMat))) dplyr::rename(., !!colnames(yMat) := 's0')
      else .
      } %>%
    as_tibble(.) %>%
    dplyr::mutate(., coefname = replace(coefname, coefname == '(Intercept)', 'constant')) %>%
    replace(., is.na(.), 0) %>%
    # Move constant row to front
    {dplyr::bind_rows(
      dplyr::filter(., coefname == 'constant'),
      dplyr::filter(., coefname != 'constant')
    )} %>%
    dplyr::filter(., .[[2]] != 0) # Remove shrunk coefficients
  
  # Get fitted values
  xMatFit =
    as_tibble(xMat) %>%
    dplyr::bind_cols(constant = rep(1, nrow(.)), .) %>%
    dplyr::select(., coef$coefname) %>%
    as.matrix(.)
  
  coefMatFit =
    coef %>%
    dplyr::select(., -coefname) %>%
    as.matrix(.)
  
  fit =
    xMatFit %*% coefMatFit %>%
    as_tibble(.)
  
  resids =
    tibble(hist = as_tibble(yMat)[[1]], fit = fit[[1]]) %>%
    dplyr::mutate(., !!colnames(fit) := hist - fit) %>%
    dplyr::select(., -hist, -fit)
  
  obs = nrow(resids)
  sse = as.matrix(resids) %>% .^{2} %>% sum(.)
  mse = sse/obs
  rmse = mse^.5
  mae = sum(abs(as.matrix(resids)))/obs
  
  coefList = (coef %>% dplyr::filter(., .[,2] != 0) %>% .[['coefname']])
  
  if (verbose == TRUE) message('Coefficients shrunk to: ',paste0(coefList, collapse = ', '))
  
  list(
    cvTable = cvTable,
    cvChart = cvChart,
    optim = optim,
    elasticNetReg = elasticNetReg,
    coef = coef,
    fit = fit,
    resids = resids,
    obs = obs,
    sse = sse,
    mse = mse,
    rmse = rmse,
    mae = mae,
    xMat = xMat,
    yMat = yMat,
    bigK = bigK,
    bigJ = bigJ,
    alphaSearchGrid = alphaSearchGrid
    ) %>%
    return(.)
}
