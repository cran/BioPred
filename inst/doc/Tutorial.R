## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----include = FALSE----------------------------------------------------------
library(knitr)
library (kableExtra)
knitr::opts_chunk$set(warning = FALSE, message = FALSE) 

## ----echo=TRUE----------------------------------------------------------------
library(BioPred)

kable_styling(kable(x=head(tutorial_data), format="html", caption = "The first 6 subjects"),
              bootstrap_options="striped",font_size=16)

rownames(tutorial_data)=NULL

## ----echo=TRUE----------------------------------------------------------------
X_feature=tutorial_data[,c("x1", "x2" ,"x3" ,"x4","x5","x6","x7","x8","x9","x10")]
y_label=tutorial_data$y.con
# Convert treatment variable to (1 -1) format (1 for treatment, -1 for control)
trt=ifelse(tutorial_data$treatment==1, 1, -1)
true_label=tutorial_data$subgroup_label
# Estimate the propensity score using logistic regression
data_logit=tutorial_data[,c("treatment","x1", "x2" ,"x3" ,"x4" , "x5", "x6" , "x7","x8","x9","x10")]
logit_model <- glm(treatment~ ., data = data_logit, family = binomial)
pi <- predict(logit_model, type = "response")
# Train the XGBoostSub_con model using the specified parameters
model <- XGBoostSub_con(X_feature, y_label, trt ,pi,Loss_type = "A_learning", params = list(learning_rate = 0.005, max_depth = 1, lambda = 5, tree_method = 'hist'), nrounds = 200, disable_default_eval_metric = 0, verbose = FALSE)

## ----echo=TRUE----------------------------------------------------------------
cat("Evaluation loss:\n")
print(model$evaluation_log)

## ----echo=TRUE----------------------------------------------------------------
eval_metric_test <- eval_metric_con(model, X_feature, y_label, pi, trt, Loss_type = "A_learning")
cat("Testing Evaluation Result:\n")
print(eval_metric_test$value)


## ----echo=TRUE----------------------------------------------------------------
biomarker_imp=predictive_biomarker_imp(model)
kable_styling(kable(x=biomarker_imp, format="html", caption = "Biomarker importance table"),
              bootstrap_options="striped",font_size=16)

## ----echo=TRUE----------------------------------------------------------------
subgroup_results=get_subgroup_results(model, X_feature, subgroup_label=true_label, cutoff = 0.5)
kable_styling(kable(x=head(subgroup_results$assignment), format="html", caption = "The first 6 subjects subgroup assigentment"),
              bootstrap_options="striped",font_size=16)
cat("Prediction accuracy of true subgroup label:\n")
cat(subgroup_results$acc)



## ----echo=TRUE,fig.height=6, fig.width=6--------------------------------------
cdf_plot (xvar="x2", data=tutorial_data, y.int=5, xlim=NULL, xvar.display="Biomarker_X2", group=NULL)

## ----echo=TRUE,fig.height=6, fig.width=8--------------------------------------
scat_cont_plot(
  yvar='y.con',
  xvar='x2',
  data=tutorial_data,
  ybreaks = NULL,
  xbreaks = NULL,
  yvar.display = 'y',
  xvar.display = "Biomarker_X2"
)

## ----echo=TRUE,fig.height=6, fig.width=6--------------------------------------
cutoff_result_con=fixcut_con(yvar='y.con', xvar="x2", dir=">", cutoffs=c(-0.1,-0.3,-0.5,0.1,0.3,0.5), data=tutorial_data, method="t.test", yvar.display='y.con', xvar.display='biomarker x2', vert.x=F)

cutoff_result_con$fig


## ----echo=TRUE----------------------------------------------------------------
res=cut_perf (yvar="y.con", censorvar=NULL, xvar="x2", cutoff=c(0.5), dir="<=", xvars.adj=NULL, data=tutorial_data, type='c', yvar.display='y.con', xvar.display="biomarker x2")
kable_styling(kable(x=res$comp.res.display[,-2:-4], format="html", caption = "Performace at optimal cutoff"),
              bootstrap_options="striped",font_size=16)


## ----echo=TRUE,fig.height=10, fig.width=15------------------------------------
tutorial_data$biogroup=ifelse(tutorial_data$x2<=0.5,'biomarker_positive','biomarker_negative')

res = subgrp_perf_pred(yvar="y.time", censorvar="y.event", grpvar="biogroup", grpname=c("biomarker_positive",'biomarker_negative'),trtvar="treatment_categorical", trtname=c("Placebo" , "Treatment"), xvars.adj=NULL,data=tutorial_data, type="s")
kable_styling(kable(x=res$group.res.display, format="html", caption = "BioSubgroup Stat"),
              bootstrap_options="striped",font_size=16)
res$fig


## ----echo=TRUE----------------------------------------------------------------
y_label=tutorial_data$y.bin
model <- XGBoostSub_bin(X_feature, y_label, trt ,pi,Loss_type = "A_learning", params = list(learning_rate = 0.01, max_depth = 1, lambda = 5, tree_method = 'hist'), nrounds = 300, disable_default_eval_metric = 0, verbose = FALSE)

## ----echo=TRUE----------------------------------------------------------------
biomarker_imp=predictive_biomarker_imp(model)
kable_styling(kable(x=biomarker_imp, format="html", caption = "Biomarker importance table"),
              bootstrap_options="striped",font_size=16)

## ----echo=TRUE----------------------------------------------------------------
subgroup_results=get_subgroup_results(model, X_feature, subgroup_label=true_label, cutoff = 0.5)
kable_styling(kable(x=head(subgroup_results$assignment), format="html", caption = "The first 6 subjects subgroup assigentment"),
              bootstrap_options="striped",font_size=16)
cat("Prediction accuracy of true subgroup label:\n")
cat(subgroup_results$acc)

## ----echo=TRUE,fig.height=6, fig.width=6--------------------------------------
roc_bin_plot (yvar='y.bin', xvars="x10", dirs="auto", data=tutorial_data, yvar.display='y.bin', xvars.display="Biomarker x10")

## ----echo=TRUE,fig.height=6, fig.width=6--------------------------------------
roc_bin_plot (yvar='y.bin', xvars="x2", dirs="auto", data=tutorial_data, yvar.display='y.bin', xvars.display="Biomarker x2")

## ----echo=TRUE,fig.height=6, fig.width=9--------------------------------------
cutoff_result_bin=fixcut_bin (yvar='y.bin', xvar="x10", dir=">", cutoffs=c(-0.1,-0.3,-0.5,0.1,0.3,0.5), data=tutorial_data, method="Fisher", yvar.display='Response_Y', xvar.display='Biomarker_X10', vert.x=F)
cutoff_result_bin$fig

## ----echo=TRUE----------------------------------------------------------------
res=cut_perf (yvar="y.con", censorvar=NULL, xvar="x10", cutoff=c(0.5), dir=">", xvars.adj=NULL, data=tutorial_data, type='c', yvar.display='y.con', xvar.display="biomarker x2")
kable_styling(kable(x=res$comp.res.display[,-2:-4], format="html", caption = "caption"),
              bootstrap_options="striped",font_size=16)


## ----echo=TRUE,fig.height=10, fig.width=15------------------------------------
tutorial_data$biogroup=ifelse(tutorial_data$x10>0.5,'biomarker_positive','biomarker_negative')

res = subgrp_perf_pred(yvar="y.time", censorvar="y.event", grpvar="biogroup", grpname=c("biomarker_positive",'biomarker_negative'),trtvar="treatment_categorical", trtname=c("Placebo" , "Treatment"), xvars.adj=NULL,data=tutorial_data, type="s")
kable_styling(kable(x=res$group.res.display, format="html", caption = "BioSubgroup Stat"),
              bootstrap_options="striped",font_size=16)
res$fig


## ----echo=TRUE----------------------------------------------------------------
# Prepare the data for training the XGBoostSub_sur
y_label=tutorial_data$y.time
delta=tutorial_data$y.event

# Train the XGBoostSub_sur model using the specified parameters
model <- XGBoostSub_sur(X_feature, y_label, trt ,pi,delta,Loss_type = "A_learning", params=list(learning_rate = 0.005, max_depth = 1,lambda = 9, gamma=1, min_child_weight=1,
max_bin=128, tree_method = 'exact',subsample=0.8), nrounds = 250, disable_default_eval_metric = 1, verbose = FALSE)

## ----echo=TRUE----------------------------------------------------------------
biomarker_imp=predictive_biomarker_imp(model)
kable_styling(kable(x=biomarker_imp, format="html", caption = "Biomarker importance table"),
              bootstrap_options="striped",font_size=16)

## ----echo=TRUE----------------------------------------------------------------
subgroup_results=get_subgroup_results(model, X_feature, subgroup_label=true_label, cutoff = 0.5)
kable_styling(kable(x=head(subgroup_results$assignment), format="html", caption = "The first 6 subjects subgroup assigentment"),
              bootstrap_options="striped",font_size=16)
cat("Prediction accuracy of true subgroup label:\n")
cat(subgroup_results$acc)

## ----echo=TRUE,fig.height=10, fig.width=10------------------------------------
res = cat_summary(yvar="risk_category", yname=c("High Risk","Intermediate Risk" ,"Low Risk"), xvars="treatment_categorical", xname.list=list(c("Placebo" , "Treatment")), data=tutorial_data)

kable_styling(kable(x=res$cont.display, format="html", caption = "Contingency Table"),
              bootstrap_options="striped",font_size=16)

## ----echo=TRUE,fig.height=10, fig.width=15------------------------------------
res = subgrp_perf_pred(yvar="y.time", censorvar="y.event", grpvar="risk_category", grpname=c("High Risk","Intermediate Risk" ,"Low Risk"),trtvar="treatment_categorical", trtname=c("Placebo" , "Treatment"), xvars.adj=NULL,data=tutorial_data, type="s")
kable_styling(kable(x=res$group.res.display, format="html", caption = "Subgroup Stat"),
              bootstrap_options="striped",font_size=16)
res$fig

