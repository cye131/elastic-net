################### Load Libraries & Inputs ###################
library(tidyverse)
library(glmnet)
library(foreach)
library(doParallel)

chimpsRdsInput = 'R:/regulatory/Stress Testing 2019/Data/2018Q4/chimpsData_2018Q4.rds'
dependentXlsxInput = 'R:/regulatory/Stress Testing 2019/Data/2018Q4/CHIMPs Credit NewDev 4_17_19.xlsx'

chimpsRds = readRDS(chimpsRdsInput)



################### Create Dependent & Independent Matrices ###################
# Get indep variables
yDfFull =
  readxl::read_xlsx(dependentXlsxInput, sheet = 'To MDG', skip = 0) %>%
  {dplyr::bind_cols(
    strDate = paste0('20', stringr::str_sub(.[[1]], 3, 5), '.', stringr::str_sub(.[[1]], 1, 1)),
    .[2:ncol(.)]
  )} %>%
  na.omit(.)

# Extract dates
historicalDates = yDfFull$strDate

# Seperate indep variables into a list of matrices, each matrix corresponding to one indep variable
yDfList =
  dplyr::filter(yDfFull, strDate %in% historicalDates) %>%
  dplyr::select(., -strDate) %>%
  {lapply(colnames(.) %>% setNames(., .), function(colname) dplyr::select(yDfFull, strDate, colname))}

  #dplyr::select(., -strDate) %>%
  # {lapply(colnames(.) %>% setNames(., .), function(colname) .[, colname] %>% as.matrix(.))}

# Create single dep variable matrix with all CHIMPS variables and their lags (up to 4)
xDf0 =
  chimpsRds[[1]]$df %>%
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
  dplyr::bind_cols(strDate = chimpsRds[[1]]$df$strDate, .) %>%
  dplyr::filter(., strDate %in% historicalDates) %>%
  # Remove any columns with NAs
  purrr::keep(., function(col) length(col[is.na(col)]) == 0)



################### Get Elastic Net Results ###################
# Create parallel clusters
cl = makeCluster(20)
registerDoParallel(cl)

# Let each core do one elastic net regression
elasticNetResultsList =
  foreach(yDf = yDfList, .inorder = FALSE, .packages = c('glmnet', 'tidyverse')) %dopar% {
  # lapply(yDfList[1:5], function(yDf) {
    
    xDf1 =
      dplyr::select(yDf, -strDate) %>%
      {lapply(1:4, function(l)
        lapply(., function(col) dplyr::lag(col, l)) %>%
          dplyr::bind_cols(.) %>%
          setNames(., paste0(colnames(.), '.l', l))
        )} %>%
      dplyr::bind_cols(strDate = yDf$strDate, .)
  
    xDf =
      list(xDf1, xDf0) %>%
      purrr::reduce(., function(df1, df2) dplyr::inner_join(df1, df2, by = 'strDate')) %>%
      na.omit(.)

    xMat =
      xDf %>% dplyr::select(., -strDate) %>% as.matrix(.)

    yMat =
      yDf %>% dplyr::filter(., strDate %in% xDf$strDate) %>% dplyr::select(., -strDate) %>% as.matrix(.)

    getElasticNet(xMat = xMat, yMat = yMat, bigK = 10, bigJ = 50, alphaSearchGrid = 1) %>%
      return(.)
  # })
  } %>%
  setNames(., names(yDfList))

# Close clusters
stopCluster(cl)



################### Get Fit Error Diagnostics ###################
errorDiagnostics =
  lapply(names(elasticNetResultsList) %>% setNames(., .), function(varname) {
    elasticNetResult = elasticNetResultsList[[varname]]
    tibble(
      Variable = varname,
      'Possible Covariates' = ncol(elasticNetResult$xMat) + 1,
      'Shrunk Covariates' = nrow(elasticNetResult$coef),
      Obs = elasticNetResult$obs,
      SSE = elasticNetResult$sse,
      MSE = elasticNetResult$mse,
      RMSE = elasticNetResult$rmse,
      MAE = elasticNetResult$mae
      ) %>%
      return(.)
  }) %>%
  dplyr::bind_rows(.)
  







