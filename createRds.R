################### Load Libraries & Inputs ###################
library(tidyverse)
library(glmnet)
library(foreach)
library(doParallel)

# chimpsRds = readRDS('R:/regulatory/Stress Testing 2019/Data/2018Q4/chimpsData_2018Q4.rds')
oldModels =
  list.files('R:/Projects/ArimanatorRunner/output/Stress Test 2019', pattern = '.rds', full.names = TRUE, ignore.case = TRUE) %>%
  setNames(., basename(.)) %>%
  lapply(., function(rds)
    readRDS(rds)
  )

fn = new.env()
list.files('D:/Ye/chimps/R/functions', pattern = '[.]R$', full.name = TRUE, recursive = TRUE) %>%
  lapply(., function(R)
    source(R)
    )
  

################### Create Dependent & Independent Matrices ###################
local({
  
  newModels =
    
    lapply(oldModels, function(.model) {

      # Seperate indep variables into a list of matrices, each matrix corresponding to one indep variable
      yDf =
        .model$dv$data$df %>%
          dplyr::mutate(
            .,
            Value =
              {
                if (is.null(.model$dv$data$params$valDiff)) .$Value
                else if (.model$dv$data$params$valDiff$diff == 1 & .model$dv$data$params$valDiff$log == TRUE)
                  log(.$Value/dplyr::lag(.$Value, 1))
                else if (.model$dv$data$params$valDiff$diff == 1) .$Value - dplyr::lag(.$Value, 1)
                else if (.model$dv$data$params$valDiff$log == TRUE) log(.$Value)
                else .$Value
              }
          ) %>%
          na.omit(.) %>%
          setNames(., c('strDate', .model$dv$name))

      
      # Create single dep variable matrix with all CHIMPS variables and their lags (up to 4)
      xDf0 =
        chimpsRds$base$df %>%
        dplyr::select(., -strDate) %>%
        # Conduct transformations
        {lapply(colnames(.), function(varname) {
          if (chimpsRds[[1]]$params[[varname]]$diff == 1 & chimpsRds[[1]]$params[[varname]]$log == TRUE) {
            tibble(!!varname := log(.[[varname]]/dplyr::lag(.[[varname]], 1)))
          } else if (chimpsRds[[1]]$params[[varname]]$diff == 1) {
            tibble(!!varname := .[[varname]] - dplyr::lag(.[[varname]], 1))
          } else if (chimpsRds[[1]]$params[[varname]]$log == 1) {
            tibble(!!varname := log(.[[varname]]))
          } else {
            tibble(!!varname := .[[varname]]) 
          }
        })} %>%
        dplyr::bind_cols(.) %>%
        # Create lagged columns
        dplyr::bind_cols(
          .,
          lapply(1:4, function(l)
            lapply(., function(col) dplyr::lag(col, l)) %>%
              dplyr::bind_cols(.) %>%
              setNames(., paste0(colnames(.), '.l', l))
          )
        ) %>%
        # Filter by dates
        dplyr::bind_cols(strDate = chimpsRds[[1]]$df$strDate, .) #%>%
        # Remove any columns with NAs
        # purrr::keep(., function(col) length(col[is.na(col)]) == 0)
      
    })
  
})



################### Fit Elastic Net ###################
local({
  # Create parallel clusters
  cl = makeCluster(20)
  registerDoParallel(cl)
  
  # Let each core do one elastic net regression
  elasticNetResults =
    foreach(yDf = yDfList, .inorder = FALSE, .packages = c('glmnet', 'tidyverse')) %dopar% {
      # lapply(yDfList[1:5], function(yDf) {
      
      xDf1 =
        dplyr::select(yDf, -strDate) %>%
        {lapply(1:4, function(l)
          lapply(., function(col) dplyr::lag(col, l)) %>%
            dplyr::bind_cols(.) %>%
            setNames(., paste0(colnames(.), '.l', l))
        )} %>%
        dplyr::bind_cols(strDate = yDf$strDate, .) %>%
        na.omit(.)
  
      xDf =
        dplyr::inner_join(
          xDf1,
          xDf0 %>%
            dplyr::filter(., strDate %in% xDf1$strDate) %>%
            purrr::keep(., function(col) length(col[is.na(col)]) == 0),
          by = 'strDate')
      
      xMat =
        xDf %>% dplyr::select(., -strDate) %>% as.matrix(.)
      
      yMat =
        yDf %>% dplyr::filter(., strDate %in% xDf$strDate) %>% dplyr::select(., -strDate) %>% as.matrix(.)
      
      getElasticNet(xMat = xMat, yMat = yMat, bigK = 10, bigJ = 50, alphaSearchGrid = c(1)) %>%
        c(., list(dates = xDf$strDate)) %>%
        return(.)
      # })
    } %>%
    setNames(., names(yDfList))
  
  # Close clusters
  stopCluster(cl)
  
  elasticNetResults <<- elasticNetResults
})



################### Get Fit Error Diagnostics ###################
local({
  
  fitDiagnostics =
    lapply(names(elasticNetResults) %>% setNames(., .), function(varname) {
      elasticNetResult = elasticNetResults[[varname]]
      tibble(
        Variable = varname,
        'Possible Covariates' = ncol(elasticNetResult$xMat) + 1,
        'Shrunk Covariates' = nrow(elasticNetResult$coef),
        Alpha = elasticNetResult$optim$alpha,
        Lambda = elasticNetResult$optim$lambda,
        Obs = elasticNetResult$obs,
        SSE = elasticNetResult$sse,
        MSE = elasticNetResult$mse,
        RMSE = elasticNetResult$rmse,
        MAE = elasticNetResult$mae
      ) %>%
        return(.)
    }) %>%
    dplyr::bind_rows(.)
  
  fitDiagnostics <<- fitDiagnostics
})



################### Get Forecasts ###################
local({
  
  pred =
    lapply((elasticNetResults), function(elasticNetRes) {

      predDf =
        as_tibble(elasticNetRes$yMat) %>%
        dplyr::bind_cols(strDate = elasticNetRes$dates, .) %>%
        tail(., 4) %>%
        dplyr::bind_rows(
          .,
          tibble(
            strDate =
              tail(.$strDate, 1) %>%
              {paste0(
                stringr::str_sub(., 1, 4),
                stringr::str_pad(as.numeric(stringr::str_sub(., -1)) * 3, 2, pad = '0'),
                '01'
              )} %>%
              lubridate::ymd(.) %>%
              {lubridate::ceiling_date(., 'month') - lubridate::days(1)} %>%
              {seq(from = . + lubridate::days(1), by = '+3 month', length.out = 20 + 1)[2:21] - lubridate::days(1)} %>%
              as.Date(., origin = lubridate::origin) %>%
              {paste0(lubridate::year(.), '.', lubridate::quarter(.))}
          )
        )
      
      for (t in 1:(nrow(predDf) - 4)) {
        
        rowIndex = t + 4 # Row index of the new date
        
        # Get lagged covariate
        xLag =
          predDf[(rowIndex - 1):(rowIndex - 4), 2] %>%
          t(.) %>%
          as_tibble(.) %>%
          setNames(., paste0(colnames(predDf)[[2]],'.l',1:4)) %>%
          dplyr::bind_cols(constant = 1, .)
        
        xExog =
          dplyr::filter(xDf0, strDate == predDf[[rowIndex, 'strDate']]) %>%
          dplyr::select(., -strDate)
        
        x =
          dplyr::bind_cols(xLag, xExog) %>%
          dplyr::select(., elasticNetRes$coef$coefname)
        
        if (all(elasticNetRes$coef$coefname == colnames(x)) != TRUE)
          stop('Error - b matrix and x matrix colnames not identical!')
        
        coefMat = as.matrix(elasticNetRes$coef[[2]])
        xMat = as.matrix(x)
        
        y = (t(coefMat) %*% t(xMat))
        
        # Add forecast date and append to prediction matrix
        predDf[[rowIndex, 2]] = y
      }
      
      predDf %>% .[5:nrow(.),] %>% return(.)
      
    })
  
  
  pred <<- pred
})


################### Fit Holdout ###################
local({
  
  # Create parallel clusters
  cl = makeCluster(20)
  registerDoParallel(cl)
  
  # Let each core do one elastic net regression
  holdoutElasticNetResult =
    foreach(elasticNetRes = elasticNetResults, .inorder = FALSE, .packages = c('glmnet', 'tidyverse')) %dopar% {
      
      xMat = elasticNetRes$xMat %>% as_tibble(.) %>% .[1:(nrow(.) - 4), ] %>% as.matrix(.)
      yMat = elasticNetRes$yMat %>% as_tibble(.) %>% .[1:(nrow(.) - 4), ] %>% as.matrix(.)
      
      getElasticNet(xMat = xMat, yMat = yMat, bigK = 10, bigJ = 50, alphaSearchGrid = c(1)) %>%
        c(., list(dates = elasticNetRes$dates %>% .[1:(length(.) - 4)])) %>%
        return(.)
      
    } %>%
    setNames(., names(yDfList))
  
  # Close clusters
  stopCluster(cl)
  
  holdoutElasticNetResult <<- holdoutElasticNetResult
})



################### Predict Holdout ###################
local({
  holdoutPred =
    lapply((holdoutElasticNetResult), function(holdoutRes) {
      
      predDf =
        as_tibble(holdoutRes$yMat) %>%
        dplyr::bind_cols(strDate = holdoutRes$dates, .) %>%
        tail(., 4) %>%
        dplyr::bind_rows(
          .,
          tibble(
            strDate =
              tail(.$strDate, 1) %>%
              {paste0(
                stringr::str_sub(., 1, 4),
                stringr::str_pad(as.numeric(stringr::str_sub(., -1)) * 3, 2, pad = '0'),
                '01'
                )} %>%
              lubridate::ymd(.) %>%
              {lubridate::ceiling_date(., 'month') - lubridate::days(1)} %>%
              {seq(from = . + lubridate::days(1), by = '+3 month', length.out = 20 + 1)[2:21] - lubridate::days(1)} %>%
              as.Date(., origin = lubridate::origin) %>%
              {paste0(lubridate::year(.), '.', lubridate::quarter(.))}
            )
          )
      
      for (t in 1:(nrow(predDf) - 4)) {
        
        rowIndex = t + 4 # Row index of the new date
        
        # Get lagged covariate
        xLag =
          predDf[(rowIndex - 1):(rowIndex - 4), 2] %>%
          t(.) %>%
          as_tibble(.) %>%
          setNames(., paste0(colnames(predDf)[[2]],'.l',1:4)) %>%
          dplyr::bind_cols(constant = 1, .)
        
        xExog =
          dplyr::filter(xDf0, strDate == predDf[[rowIndex, 'strDate']]) %>%
          dplyr::select(., -strDate)
        
        x =
          dplyr::bind_cols(xLag, xExog) %>%
          dplyr::select(., holdoutRes$coef$coefname)
        
        if (all(holdoutRes$coef$coefname == colnames(x)) != TRUE)
          stop('Error - b matrix and x matrix colnames not identical!')
        
        coefMat = as.matrix(holdoutRes$coef[[2]])
        xMat = as.matrix(x)
        
        y = (t(coefMat) %*% t(xMat))
  
        # Add forecast date and append to prediction matrix
        predDf[[rowIndex, 2]] = y
      }
      
      predDf %>% .[5:nrow(.),] %>% return(.)
      
    })
  
  holdoutPred <<- holdoutPred
})



################### Calculate Holdout Diagnostics ###################
local({
  
  holdoutDiagnostics =
    lapply(names(holdoutForecasts) %>% setNames(., .), function(name) {
      
      df =
        dplyr::inner_join(holdoutForecasts[[name]], yDfList[[name]], by = 'strDate') %>%
        setNames(., c('strDate', 'forecast', 'historical')) %>%
        dplyr::mutate(., resid = historical - forecast)
      
      resids = dplyr::select(df, resid)
      
      obs = nrow(resids)
      sse = as.matrix(resids) %>% .^{2} %>% sum(.)
      mse = sse/obs
      rmse = mse^.5
      mae = sum(abs(as.matrix(resids)))/obs
      
      list(df = df, resids = resids, obs = obs, sse = sse, mse = mse, rmse = rmse, mae = mae)
      })
  
  holdoutDiagnostics <<- holdoutDiagnostics
})









################### Export ###################
rdsModels =
  lapply(names(arimaRdsList) %>% setNames(., .), function(.rdsname) {
    
    newModel = arimaRdsList[[.rdsName]]
    newModel$dv$diagnostics$plot_ccf =
      dplyr::bind_cols(
        elasticNetResults[[.rdsName]]$yMat %>% as_tibble(.),
        elasticNetResults[[.rdsName]]$xMat %>%
          as_tibble(.) %>%
          dplyr::select(., elasticNetResults[[.rdsName]]$coef$coefname %>% .[. != 'constant'])
      ) %>% plot_ccf(., model$dv$name)
    
    # list(
    #   SNVRegulatory1.5::plot_ccf()
    #   # elasticNetResults = elasticNetResults[[.rdsname]],
    #   # elasticNetHoldoutResults = elasticNetHoldoutResults[[.rdsname]],
    #   # holdoutForecasts = holdoutForecasts[[.rdsname]],
    #   # holdoutDiagnostics = holdoutDiagnostics[[.rdsname]]
    #   
    # )
    return(newModel)
    
  })




