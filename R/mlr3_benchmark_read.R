## DLBC 2 fold
library(mlr3)
library(mlr3proba)
library(mlr3verse)
library(mlr3tuningspaces)
library(reticulate)
library(survivalmodels)
library(plyr)
library(dplyr)
library(survival)
library(survminer)
library(edgeR)
library(DESeq2)
library(foreach)
library(doMC)
library(data.table)
library(future)
library(future.apply)
# library(caret)

# options(future.globals.maxSize= 4294967296)
# options(expressions = 5e5)
# progressr::handlers(global = TRUE)

# path_to_python <- "/home/rstudio/.local/share/r-miniconda/envs/r-reticulate/bin/python"
# use_python(path_to_python)
# install_pycox(
#   method = "auto",
#   conda = "auto",
#   pip = TRUE,
#   install_torch = TRUE
# )
#install_keras(pip = TRUE, install_tensorflow = TRUE)
#install.packages("devtools")
# devtools::install_github("rstudio/keras")
# devtools::install_github("rstudio/tensorflow")
# library(keras)
# library(tensorflow)
# library(pseudo)



# t_harmonised <- read.csv("cancerData/DLBC_ProteinCoding.csv",header=TRUE,sep=",")

lgr::get_logger("mlr3")$set_threshold("warn")


plan(multisession)

for (index in 1:30){
  
  tvsd_train2 <- read.csv(paste("Results_CompSurv/DLBC_ProteinCoding/preFeatSel/train_set", index, ".csv", sep = ""))
  tvsd_test2 <- read.csv(paste("Results_CompSurv/DLBC_ProteinCoding/preFeatSel/test_set", index, ".csv", sep = ""))
  
  tdatam_fs <- TaskSurv$new("tvsd_train2", tvsd_train2[,-1], time = "time", event = "status")
  task_fs <- tdatam_fs
  
  #### CREATE LEARNERS ####
  learners = list(
    # lrn("surv.akritas", id = "akritas", lambda = 1, reverse = FALSE), # 0.5 regardless of parameters. does not work
    lrn("surv.blackboost", id = "blackboost", mstop = 100), # almost perfect, error on 1/25 #
    lrn("surv.cforest", id = "cforest", ntree = 100), # error on 4/25, good enough
    lrn("surv.coxboost", id = "coxboost"), # error on 4/25 good enough
    # # lrn("surv.coxtime", id = "coxtime", optimizer = "adam", early_stopping = TRUE), # needs to be re tested, takes very long
    lrn("surv.ctree", id = "ctree"), # error on 8/25, okayish but needs further testing 229 ft
    # # lrn("surv.deephit", id = "deephit"),
    # lrn("surv.deepsurv", id = "deepsurv"), # works on local
    # # lrn("surv.dnnsurv", id = "dnnsurv", optimizer = "adam"), # returns 0.5 cindex for 24 results out of 25, may need to tune more impactful parameters
    # # lrn("surv.flexible", id = "flexible"), # doesn't work, future error / C stack usage
    # lrn("surv.gamboost", id = "gamboost"), # doesn't work, future error ++ /C stack usage
    lrn("surv.gbm", id = "gbm", bag.fraction = 0.9), # error on 4/25
    lrn("surv.glmboost", id = "glmboost"), # error on 4/25, it takes very long for this to be complete
    lrn("surv.glmnet", id = "lasso", alpha = 1,s = 0.01), # error on 4/25
    lrn("surv.glmnet", id = "ridge", alpha = 0,s = 0.01), # error on 6/25
    lrn("surv.glmnet", id = "elasticnet", alpha = 0.5,s = 0.01), # error on 4/25
    # lrn("surv.loghaz", id= "loghaz", optimizer = "adamw", frac = 0.3) # takes very long and outputs garbage results/ 352 fts
    # lrn("surv.mboost", id= "mboost", baselearner = "bols", center = TRUE), # error ++ future error/ c stack usage too close to the limit
    lrn("surv.obliqueRSF", id= "obliqueRSF"), # error on 4/25
    # # lrn("surv.parametric", id= "parametric", dist = "exponential"), # corrputed double linked list 
    # # lrn("surv.pchazard", id= "pchazard"), # runs on local
    lrn("surv.penalized", id= "penalized", lambda1=10, lambda2=10), # should work now but tuning is irrelevant
    lrn("surv.ranger", id= "ranger"), # should be working
    lrn("surv.rfsrc", id= "rfsrc"), # 5/25
    lrn("surv.rpart", id= "rpart"), # 2/25, cindex is lower compared to other learners but still above 0.75 range
    lrn("surv.svm", id= "svm", type = "hybrid", diff.meth = "makediff3", kernel = "lin_kernel", gamma.mu = c(100,1000)), # no tunable params, still works fine
    lrn("surv.xgboost", id = "xgboost_gbtree", booster = "gbtree"), # 2/25
    lrn("surv.xgboost", id = "xgboost_gblinear", booster = "gblinear"), # 7/25
    lrn("surv.xgboost", id = "xgboost_dart", booster = "dart") # 6/25
  )
  
  # from glmnet C++ code (error code -30083); Numerical error at 83th lambda value; solutions for larger values of lambda returned
  # This happened PipeOp lasso's $train() 
  # from glmnet C++ code (error code -77); Convergence for 77th lambda value not reached after maxit=100000 iterations; solutions for larger lambdas returned
  # This happened PipeOp lasso's $train() 
  
  #### CREATE TUNE SPACES ####
  
  tune_spaces <- list(
    
    # tune_ps_akritas = ParamSet$new(
    #   list(
    #     ParamDbl$new("lambda", lower = 0, upper = 1)
    #   )),
    # 
    tune_ps_blackboost = ParamSet$new(
      list(
        ParamFct$new("family", levels = c("gehan", "cindex")), # may need to implement coxph as a separate learner
        ParamInt$new("mstop", lower = 10, upper = 1000),
        ParamDbl$new("nu", lower = 0, upper = 0.1),
        ParamInt$new("mtry", lower = 1, upper = 12)
        
        # # ParamDbl$new("alpha", lower = 0, upper = 1),
        # ParamDbl$new("abseps", lower = 0, upper = 10),
        # ParamInt$new("maxdepth", lower = 0, upper = 10)
      )),
    
    tune_ps_cforest = ParamSet$new(
      list(
        ParamInt$new("ntree", lower = 250, upper = 2500),
        ParamInt$new("mtry", lower = 1, upper = 12)
        # ParamDbl$new("alpha", lower = 0, upper = 1),
        # ParamDbl$new("abseps", lower = 0, upper = 10),
        # ParamInt$new("maxdepth", lower = 0, upper = 10)
      )),
    #
    tune_ps_coxboost = ParamSet$new(
      list(
        # ParamDbl$new("penalty", lower = , upper = ), # according to the thesis, we need to put "optimCoxBoostPenalty here but this is a numeric parameter
        ParamInt$new("stepno", lower = 500, upper = 1500)
      )),
    # # 
    # # tune_ps_coxtime = ParamSet$new(
    # #   list(
    # #     ParamDbl$new("frac", lower = 0.1, upper = 1),
    # #     ParamLgl$new("standardize_time", default = TRUE),
    # #     # ParamUty$new("num_nodes")
    # #     ParamDbl$new("dropout", lower = 0.1, upper = 1),
    # #     ParamDbl$new("weight_decay", lower = 0.1, upper =  0.5),
    # #     ParamDbl$new("learning_rate", lower = 0.1, upper = 1),
    # #     ParamInt$new("epochs", lower = 50, upper = 150),
    # #     ParamLgl$new("early_stopping", default = TRUE)
    # #   )),
    # # 
    tune_ps_ctree = ParamSet$new(
      list(
        ParamDbl$new("alpha", lower = 0, upper = 1),
        ParamDbl$new("abseps", lower = 0, upper = 10),
        ParamInt$new("maxdepth", lower = 1, upper = 16)
      )),
    # # 
    # # tune_ps_dnnsurv = ParamSet$new(
    # #   list(
    # #     ParamDbl$new("validation_split", lower = 0.1, upper = 1),
    # #     ParamDbl$new("decay", lower = 0, upper =  0.5),
    # #     ParamInt$new("epochs", lower = 50, upper = 150),
    # #     ParamDbl$new("lr", lower = 0, upper = 1)
    # #     # ParamLgl$new("early_stopping", default = TRUE)
    # #   )),
    # # 
    # # # tune_ps_flexible <- ParamSet$new(
    # # #   list(
    # # #     ParamInt$new("k", lower = 1,upper = 7)
    # # #   )),
    # # 
    # # 
    # # tune_ps_gamboost = ParamSet$new(
    # #   list(
    # #     ParamFct$new("family", levels = c("gehan", "cindex")), # may need to implement coxph as another learner
    # #     ParamInt$new("mstop", lower = 10, upper = 1000),
    # #     # ParamDbl$new("nu", lower = 0, upper = 0.1),
    # #     ParamFct$new("baselearner", levels = c("bols","btree"))
    # #   )),
    # # 
    tune_ps_gbm = ParamSet$new(
      list(
        # ParamDbl$new("shrinkage", lower = 0, upper = 2),
        ParamInt$new("interaction.depth", lower = 1, upper = 16)
      )),
    # # 
    tune_ps_glmboost = ParamSet$new(
      list(
        ParamFct$new("family", levels = c("gehan", "cindex")), # may need to implement coxph as another learner
        ParamInt$new("mstop", lower = 10, upper = 1000),
        ParamDbl$new("nu", lower = 0, upper = 0.1)
        # ParamDbl$new("sigma", lower = 0, upper = 1)
      )),
    # # 
    tune_ps_lasso = ParamSet$new(
      list(
        ParamDbl$new("lambda.min.ratio", lower = 0, upper = 1)
        # ParamInt$new("s", lower = 0, upper = 1)
      )),
    #
    tune_ps_ridge = ParamSet$new(
      list(
        ParamDbl$new("lambda.min.ratio", lower = 0, upper = 1)
        # ParamInt$new("s", lower = 0, upper = 1)
      )),
    #
    tune_ps_elasticnet = ParamSet$new(
      list(
        ParamDbl$new("lambda.min.ratio", lower = 0, upper = 1)
        # ParamInt$new("s", lower = 0, upper = 1)
      )),
    # 
    # tune_ps_kaplan = ParamSet$new(
    #   list(
    #   )),
    
    # tune_ps_loghaz = ParamSet$new(
    #   list(
    #     # ParamDbl$new("frac", lower = 0.1, upper = 1),
    #     # ParamLgl$new("standardize_time", default = TRUE),
    #     # ParamUty$new("num_nodes")
    #     ParamDbl$new("dropout", lower = 0.1, upper = 1),
    #     ParamDbl$new("weight_decay", lower = 0.1, upper =  0.5),
    #     ParamDbl$new("learning_rate", lower = 0.1, upper = 1),
    #     ParamInt$new("epochs", lower = 50, upper = 150)
    #     # ParamLgl$new("early_stopping", default = TRUE)
    #   ))
    # 
    # tune_ps_mboost = ParamSet$new(
    #   list(
    #         # ParamFct$new("family", levels = c("gehan", "cindex")), # may need to implement coxph as another learner
    #         # ParamInt$new("mstop", lower = 10, upper = 2500),
    #         # ParamDbl$new("nu", lower = 0, upper = 0.1),
    #         ParamFct$new("baselearner", levels = "bols"),
    #         ParamLgl$new("center", default = TRUE)
    # 
    # )),
    # 
    
    # 
    # tune_ps_nelson = ParamSet$new(
    #   list(
    #   )),
    
    tune_ps_obliqueRSF = ParamSet$new(
      list(
        ParamDbl$new("alpha", lower = 0, upper = 1),
        ParamDbl$new("gamma", lower = 0, upper = 1)
      )),
    # # 
    # 
    # # tune_ps_parametric = ParamSet$new(
    # #   list(
    # #     # ParamFct$new("dist", levels = c("exponential","gaussian","lognormal","loglogistic"))
    # #   )),
    # 
    # # tune_ps_pchazard = ParamSet$new(
    # #   list(
    # #     
    # #     
    # # 
    # #   )),
    # 
    tune_ps_penalized = ParamSet$new(
      list(
        ParamDbl$new("epsilon", lower = 0, upper = 1)
        # ParamDbl$new("startbeta", lower = 1, upper = 16),
        # ParamDbl$new("startgamma", lower = 1, upper = 16),
        # ParamInt$new("steps", lower = 1, upper = 16)
        # ParamUty$new("lambda1", default = 5),
        # ParamUty$new("lambda2", default = 5)
      )),
    # # 
    tune_ps_ranger = ParamSet$new(
      list(
        ParamFct$new("splitrule", levels = "C"),
        ParamInt$new("num.trees", lower = 250, upper = 1000),
        ParamInt$new("mtry", lower = 1, upper = 12),
        ParamInt$new("min.node.size", lower = 1, upper = 20)
      )),
    #
    tune_ps_rfsrc = ParamSet$new(
      list(
        # ParamFct$new("splitrule",) #default is already logrank, may create another model for the bs.gradient version. see the thesis
        ParamInt$new("ntree", lower = 250, upper = 2500),
        ParamInt$new("mtry", lower = 1, upper = 12),
        ParamInt$new("nodesize", lower = 1, upper = 20)
      )),
    #
    tune_ps_rpart = ParamSet$new( # no params
      list(
        ParamInt$new("minbucket", lower = 1, upper = 20),
        ParamInt$new("maxdepth", lower = 2, upper = 30)
      )),
    #
    tune_ps_svm = ParamSet$new(
      list(
        # ParamFct$new("type", levels = "hybrid"),
        # ParamFct$new("diff.meth", levels = "makediff3"),
        # # ParamUty$new("gamma.mu", default = ),
        # ParamFct$new("kernel", levels = c("lin_kernel"))
        ParamInt$new("sigf", lower = 2, upper = 12),
        ParamInt$new("maxiter", lower = 20, upper = 50),
        ParamDbl$new("margin", lower = 0.01, upper = 0.1),
        ParamDbl$new("bound", lower = 5, upper = 15)
        
        # ParamDbl$new("eig.tol", lower = 0, upper = 1),
        # ParamDbl$new("conv.tol", lower = 0, upper = 1),
        # ParamDbl$new("posd.tol", lower = 0, upper = 1)
      )),
    #
    tune_ps_gbtree = ParamSet$new(
      list(
        ParamDbl$new("alpha", lower = 0, upper = 1),
        ParamDbl$new("eta", lower = 0, upper = 1),
        ParamDbl$new("gamma", lower = 0, upper = 1),
        ParamDbl$new("lambda", lower = 0, upper = 2),
        # ParamDbl$new("lambda_bias", lower = 0, upper = 2),
        ParamInt$new("nrounds", lower = 1, upper = 16) # important parameter
      )),
    #
    tune_ps_gblinear = ParamSet$new(
      list(
        ParamDbl$new("alpha", lower = 0, upper = 1),
        ParamDbl$new("eta", lower = 0, upper = 1),
        # ParamDbl$new("gamma", lower = 0, upper = 1), not used apparently
        ParamDbl$new("lambda", lower = 0, upper = 2),
        # ParamDbl$new("lambda_bias", lower = 0, upper = 2),
        ParamInt$new("nrounds", lower = 1, upper = 16) # important parameter
      )),
    
    tune_ps_dart = ParamSet$new(
      list(
        ParamDbl$new("alpha", lower = 0, upper = 1),
        ParamDbl$new("eta", lower = 0, upper = 1),
        ParamDbl$new("gamma", lower = 0, upper = 1),
        ParamDbl$new("lambda", lower = 0, upper = 2),
        # ParamDbl$new("lambda_bias", lower = 0, upper = 2),
        ParamInt$new("nrounds", lower = 1, upper = 16) # important parameter
      ))
    
  )
  
  #### FEATURE SELECTION ####
  measures = msrs(c("surv.cindex"))
  successful_fselection <- TRUE
  feature_count_df <- data_frame()
  task_list <- list()
  test_data_list <- list()
  fserrors_DLBC <- c()
  
  for (i in 1:18){
    
    tryCatch(
      expr = {
        
        evals20 = trm("evals", n_evals = 5)
        
        instance = FSelectInstanceSingleCrit$new(
          task = task_fs,
          learner = learners[[i]],
          resampling=rsmp("repeated_cv", repeats = 10, folds = 2),#rsmp("holdout"),
          measure = msr("surv.cindex"),
          terminator = evals20
        )
        
        fselector = fs("random_search")
        fselector$optimize(instance)
        
        features <- as.vector(instance$result_feature_set)
        feature_count_df <- rbind(feature_count_df, list(learners[[i]]$id, length(features)))
        gc()
        
        tvsd_train3 <- cbind(tvsd_train2[,1], subset(tvsd_train2, select=features), tvsd_train2$time, tvsd_train2$status)
        tvsd_test3 <- cbind(tvsd_test2[,1], subset(tvsd_test2, select=features), tvsd_test2$time, tvsd_test2$status)
        names(tvsd_train3)[names(tvsd_train3) == "tvsd_train2$time"] <- "time"
        names(tvsd_train3)[names(tvsd_train3) == "tvsd_train2$status"] <- "status"
        names(tvsd_test3)[names(tvsd_test3) == "tvsd_test2$time"] <- "time"
        names(tvsd_test3)[names(tvsd_test3) == "tvsd_test2$status"] <- "status"
        
        fwrite(tvsd_train3, file = paste("Results_CompSurv/DLBC_ProteinCoding/selectedfeatures/train_set", index, learners[[i]]$id, ".csv", sep = ""))
        fwrite(tvsd_test3, file = paste("Results_CompSurv/DLBC_ProteinCoding/selectedfeatures/test_set", index, learners[[i]]$id, ".csv", sep = ""))
        
        tdatam3 <- TaskSurv$new("tvsd_train3", tvsd_train3[,-1], time = "time", event = "status")
        task3 <- tdatam3
        task_list <- append(task_list, list(task3))
        test_data_list <- append(test_data_list, list(tvsd_test3))
        
      },
      error = function(e) {
        print(e)
        successful_fselection <<- FALSE
        fserrors_DLBC <<- c(fserrors_DLBC, e)
        # if feature selection goes wrong we do not attempt to train learners as it will fail too.
      }
    )
    
    
    
    
  }
  print(fserrors_DLBC)
  # feature_count_df
  
  # loghaz error
  # Rcpp::exception in py_call_impl(callable, dots$args, dots$keywords): ValueError: zero-size array to reduction operation maximum which has no identity
  
  #### TUNING ####
  if (successful_fselection){
    
    print("Started tuning the learners!")
    tuned_lrnr_results <- data_frame()
    go_next <- TRUE
    j <- 1
    while (j <= 18){
      go_next <- TRUE
      print(paste("Tuning the learner", learners[[j]]$id))
      at=AutoTuner$new(learner=learners[[j]],
                       resampling=rsmp("repeated_cv", repeats = 10, folds = 2),
                       measure=msr("surv.cindex"),
                       terminator=trm("evals", n_evals = 5),
                       tuner=tnr("random_search"),
                       search_space = tune_spaces[[j]]
      )
      timer_start_tune <- Sys.time()
      tryCatch(
        expr = {
          plan(multisession)
          at$train(task_list[[j]])
          print("Training complete")
          pred_result <- at$predict_newdata(newdata = test_data_list[[j]])$score(measures)
          timer_stop_tune <- Sys.time()
          print("Prediction Complete!")
          print(pred_result)
          print(timer_stop_tune - timer_start_tune)
          tuned_lrnr_results <- rbind(tuned_lrnr_results, list(learners[[j]]$id, pred_result, difftime(timer_stop_tune, timer_start_tune, units = "mins")))
          gc()
          # tuned_lrnr_results
        },
        error = function(e){
          # may need to do a super assignment so that the iteration is reattempted in the case of an error
          # problem with this approach is that a function may get stuck in an error loop
          # would be nice to detect if a certain error occured so we could restart based on that
          
          tuned_lrnr_results <<- rbind(tuned_lrnr_results, list(learners[[j]]$id, -1, -1))
          go_next <- FALSE 
          print(e)
        }
      )
      
      if (go_next){
        j <- j+1
      }
      
      
    }
    
    colnames(tuned_lrnr_results) <- c("Learner", "C index", "Runtime (minutes)")
    tuned_lrnr_results
    fwrite(tuned_lrnr_results, file = paste("Results_CompSurv/DLBC_ProteinCoding/cindex_time/cindex_time", index ,".csv",sep = ""))
    
    
  }
  
}

# Reading the results
# final_results_df <- data_frame()
# for (i in 1:30) {
#   
#   temp_result <- fread(paste("Results_CompSurv/DLBC_ProteinCoding/cindex_time/cindex_time", i ,".csv",sep = ""))
#   final_results_df <- rbind(final_results_df, temp_result)
# }
# final_results_df <- final_results_df[order(final_results_df$Learner),]
# colnames(final_results_df) <- c("learner", "cindex", "runtime")
# remainder_results_df <- fread("Results_CompSurv/DLBC_ProteinCoding/remainder.csv")
# final_results_df <- rbind(final_results_df,remainder_results_df)
# View(final_results_df)
# 
# # final_results_df %>%
# #   filter(learner == "blackboost") %>%
# #     summarise(x_cindex = mean(cindex))
# 
# 
# mean_cindex <- c()
# res_df <- data_frame()
# for (i in 1:18){
#   
#   snp_index <- final_results_df %>%
#     filter(learner == learners[[i]]$id) %>%
#     summarise(m_cindex = mean(cindex))
#   
#   res_df <- rbind(res_df, list(learners[[i]]$id, snp_index))
#   
# }
# 
# fwrite(res_df, file = "Results_CompSurv/DLBC_ProteinCoding/results_mean_df.csv")
# 
# 
# View(res_df)
# 
# # View(res_df)
# # mean_cindex
# # View(mean_cindex)
# # final_results_df %>%
# #   group_by(learner) %>%
# #   summarise(x_cindex = mean(cindex))
# # 
# # grouped_df <- final_results_df %>% group_by(learner)
# # 
# # grouped_df %>% summarise(x_cindex = median(cindex))
# # 
# # final_results_df %>%
# #   group_by(learner) 
# # 
# # class(final_results_df$runtime)
# # View(final_results_df)
# dim(final_results_df)
# 
# 
# total_runtime <- sum(final_results_df$runtime)
# total_runtime
# fwrite(final_results_df, file = "Results_CompSurv/final_cindex_results.csv")
