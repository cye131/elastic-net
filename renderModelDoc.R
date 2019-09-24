library('SNVRegulatory1.5')
options(rsconnect.http = 'wininet')
Sys.setenv(http_proxy = 'http://proxy.snv.net:8080')
Sys.setenv(https_proxy = 'http://proxy.snv.net:8080')
loadPkg <- function(name) {
  if (!require(name, character.only = T)) {
    install.packages(name)
    require(name, character.only = T)
  }
}

loadPkg('gridExtra')
loadPkg('tidyverse')
loadPkg('lubridate')
loadPkg('forecast')
loadPkg('zoo')
loadPkg('rmarkdown')
loadPkg('officer')
loadPkg('mschart')
loadPkg('flextable')
loadPkg('devEMF')

if (!require('captioner')) {
  loadPkg('httr')
  loadPkg('devtools')

  with_config(use_proxy(url = 'http://proxy.snv.net', port = 8080),
              install_github('adletaw/captioner'))

}


source('R:/regulatory/Stress Testing 2019/MRM/Model Documents/r/convertArimanatorRDS.R')

basePath <- 'D:/Ye/lasso/test/' #'R:/regulatory/Stress Testing 2019/MRM/Model Documents/'

#iter <- '8'

max_lags <- 5

#rds <- readRDS('R:/Projects/ArimanatorRunner/output/Stress Test 2019/model_19040_8_20190809_Fixed_Assets_NCOR.RDS')
chimpsData <- readRDS('R:/Projects/ArimanatorRunner/data/chimpsData_2018Q4.rds')
chimpsData <- chimpsData %>%
  lapply(.,
         function(x){
           for (varName in names(x$params)) {
             if (varName == 'tyield_30y') {
               x$params[[varName]] <- NULL
               next()
             } else {
               for (i in seq(1, max_lags)) {
                 newVarName <- paste0(varName, '_lag', i)
                 x$params[[newVarName]] <- x$params[[varName]]
                 x$params[[newVarName]]$lag <- i
               }
             }
           }
           return(data_ts_transform(x))
         })

for (rdsPath in dir('D:/Ye/lasso/test/inputs', '.rds', ignore.case = T, full.names = T)) {
  # rdsPath <- dir('./Stress Test 2019/CHIMPS/2019 Segmentation/', '.rds', ignore.case = T, full.names = T)[2]

  # names(chimpsData) <- c('base', names(chimpsData)[2:4], 'adverse', names(chimpsData)[6], 'severe')
  modelObj <- readRDS(rdsPath)
  # modelObj <- convertArimanatorRDS(rds, iter = iter, chimpsData = chimpsData, savePath = NULL)
  outputPath <- paste0(basePath, 'output/MRM', modelObj$id,'-',format(Sys.time(), '%Y%m%d_%H%M%S'),'.docx')

  render(input = paste0(basePath, 'template/ModelDoc.Rmd'),
         output_file = outputPath,
         params = list(rds = modelObj, outpath = outputPath))

  captions <- read_rds(str_replace(string = outputPath, '.docx', '.rds'))

  # uniARIMANormality <- tempfile(fileext = '.png')
  # ggsave(filename = uniARIMANormality,
  #        plot = marrangeGrob(modelObj$uni$diagnostics$plot_normal,  nrow = 1, ncol = 2, top = NULL))

  CCFPlots <- tempfile(fileext = '.png')
  ggsave(filename = CCFPlots, width = 6, height = 8,
         plot = marrangeGrob(modelObj$dv$diagnostics$plot_ccf, nrow = 4, ncol = 2, top = NULL))

  multiARIMANormality <- tempfile(fileext = '.png')
  ggsave(filename = multiARIMANormality, width = 6, height = 4,
         plot = marrangeGrob(modelObj$holdout$diagnostics$plot_normal,  nrow = 1, ncol = 2, top = NULL))

  FinalNormality <- tempfile(fileext = '.png')
  ggsave(filename = FinalNormality, width = 6, height = 4,
         plot = marrangeGrob(modelObj$final$diagnostics$plot_normal,  nrow = 1, ncol = 2, top = NULL))

  doc <- read_docx(outputPath) %>%
    cursor_reach(., 'uni_d0_plot') %>%
    body_add_gg(modelObj$dv$diagnostics$diff_0$plot_uni, style = 'Date', width = 6, height = 4) %>%
    # shortcuts$slip_in_plotref(depth = 1) %>%
    cursor_reach(., 'uni_d0_acf') %>%
    body_add_gg(modelObj$dv$diagnostics$diff_0$plot_acf, style = 'Date', width = 6, height = 4) %>%
    # shortcuts$slip_in_plotref(depth = 1) %>%
    cursor_reach(., 'uni_d1_plot') %>%
    body_add_gg(modelObj$dv$diagnostics$diff_1$plot_uni, style = 'Date', width = 6, height = 4) %>%
    # shortcuts$slip_in_plotref(depth = 1) %>%
    cursor_reach(., 'uni_d1_acf') %>%
    body_add_gg(modelObj$dv$diagnostics$diff_1$plot_acf, style = 'Date', width = 6, height = 4) %>%
    # shortcuts$slip_in_plotref(depth = 1) %>%
    cursor_reach(., 'CCFPlots') %>%
    body_add_img(CCFPlots, style = 'Date', width = 6, height = 8) %>%
    # shortcuts$slip_in_plotref(depth = 1) %>%
    # cursor_reach(., 'uniARIMANormality') %>%
    # body_add_img(uniARIMANormality, style = 'Date', width = 6, height = 4) %>%
    # shortcuts$slip_in_plotref(depth = 1) %>%
    cursor_reach(., 'multiARIMAACFPlots') %>%
    body_add_gg(modelObj$holdout$diagnostics$plot_acf, style = 'Date', width = 6, height = 4) %>%
    # shortcuts$slip_in_plotref(depth = 1) %>%
    cursor_reach(., 'multiARIMANormality') %>%
    body_add_img(multiARIMANormality, style = 'Date', width = 6, height = 4) %>%
    # shortcuts$slip_in_plotref(depth = 1) %>%
    cursor_reach(., 'multiARIMAInTimePlot') %>%
    body_add_gg(modelObj$holdout$diagnostics$plot_intime, style = 'Date', width = 6, height = 4) %>%
    # shortcuts$slip_in_plotref(depth = 1) %>%
    cursor_reach(., 'multiARIMAHoldoutPlot') %>%
    body_add_gg(modelObj$holdout$diagnostics$plot_holdout, style = 'Date', width = 6, height = 4) %>%
    # shortcuts$slip_in_plotref(depth = 1) %>%
    cursor_reach(., 'multiARIMAHoldoutFcasts') %>%
    body_add_gg(modelObj$holdout$diagnostics$plot_scenario, style = 'Date', width = 6, height = 4) %>%
    # shortcuts$slip_in_plotref(depth = 1) %>%
    cursor_reach(., 'FinalACFPlots') %>%
    body_add_gg(modelObj$final$diagnostics$plot_acf, style = 'Date', width = 6, height = 4) %>%
    # shortcuts$slip_in_plotref(depth = 1) %>%
    cursor_reach(., 'FinalNormality') %>%
    body_add_img(FinalNormality, style = 'Date', width = 6, height = 4) %>%
    # shortcuts$slip_in_plotref(depth = 1) %>%
    cursor_reach(., 'FinalInTimeFcasts') %>%
    body_add_gg(modelObj$final$diagnostics$plot_intime, style = 'Date', width = 6, height = 4) %>%
    # shortcuts$slip_in_plotref(depth = 1) %>%
    cursor_reach(., 'FinalScenFcasts') %>%
    body_add_gg(modelObj$final$diagnostics$plot_scenario, style = 'Date', width = 6, height = 4) %>%
    # shortcuts$slip_in_plotref(depth = 1) %>%
    cursor_reach(., 'FinalShockFcasts') %>%
    body_add_gg(modelObj$final$diagnostics$plot_shock, style = 'Date', width = 6, height = 4) %>%
    # shortcuts$slip_in_plotref(depth = 1) %>%
    body_replace_all_text(paste0(names(captions), collapse = '|'),'')

  print(doc, target = outputPath)
}

