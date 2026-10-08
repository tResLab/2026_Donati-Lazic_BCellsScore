library(tidyverse)
library(survminer)
libary(survival)
library(glmnet)


###---Input dataset has been deposited in Array Express at the accession number E-MTAB-15890, as the matrix file 

data<- as.data.frame(matrix_E_MTAB_15890)

###Data are present from the row number three

data<- data[3:773,] 

#Traspone dataset to obtain genes in columns and patients in rows

data_t<- as.data.frame(t(data))

colnames(data_t)<- data_t[1,]

data_t<- data_t[2:771,]

data_t<- data_t %>% mutate_all(as.numeric)

data_t$Sample_ID<- rownames(data_t)

rownames(data_t)<- NULL



###---the vector containing the 18 protective genes from which the score was derive

genes_Bcells<- as.vector(c("TNFRSF13B", "POU2AF1", "TNFRSF17", "CR2", "CD180", "FCER2", "BLK",
                           "CD79A", "CD22", "CD19", "CD79B", "MS4A1", "PAX5", "TNFRSF13C", "TLR9", "MEF2C", "HLA-DOB", "TLR10"))


###---Define a dataset including only the Sample_ID and the 18 protective genes from which the score was derive

data_score<- select(data_t, c("Sample_ID", one_of(genes_Bcells)))




###---To extend the score to new individual patients in clinical practice, we derived it from a matrix of Z-score-transformed gene expression values

train_z_matrix <- scale(data_score[, genes_Bcells])

# -- Save and Freeze the mean (mu) and standard deviation (sd) for each gene obtained from the training cohort data
#    These will be used to scale the gene expression data of any new patient

train_means <- attr(train_z_matrix, "scaled:center")
train_sds   <- attr(train_z_matrix, "scaled:scale")

save(train_means, train_sds, file = "scaling_params_train.RData")




###--- Clinical data can be obtained from the Metadata Template deposited in Array Express at the accession number E-MTAB-15890
#      The clinical data that are relevant to derive this score are time PFS ("characteristics: PFS_time(months)") and the PFS event ("characteristics: PFS_event"))

Metadata<- select(Metadata, c("Sample name", "characteristics: PFS_time(months)", "characteristics: PFS_event")

Metadata$'Sample name'<- Metadata$Sample_ID

# -- Check that Metadata and data_score had the same patients' order

all(data_score$Sample_ID == Metadata$Sample_ID) ##TRUE




###--- Data needed to employ the Elastic Net penalized regression model, running the glmnet package 

# Data matrix

x<-train_z_matrix 

# Survival object

y <- Surv(data_complete$PFS, data_complete$`evento PFS`)

# Sequence of alpha to test

alpha_values <- seq(0, 1, by = 0.05)

# Fixed lambda values (100 values from 10^2 to 10^-4)

common_lambda_grid <- 10^seq(2, -4, length.out = 100)

# Dataframe to save result

results_z <- data.frame(alpha = numeric(), lambda.min = numeric(), cvm.min = numeric())



###---Begin tuning loop to select the optimal alpha hyperparameter
#     Loop with strict cross-validation fold rules
#     Fix the CV folds so that each alpha evaluates EXACTLY the same patients

set.seed(42)
fold_ids <- sample(rep(seq(10), length.out = nrow(x)))

for (a in alpha_values) {
  
  cv_out <- cv.glmnet(
    x = x, 
    y = y, 
    family = "cox", 
    alpha = a, 
    lambda = common_lambda_grid, 
    foldid = fold_ids,           
    standardize = FALSE,         
    grouped = FALSE              
  )
  
  # Extract the minimum deviance obtained for this specific alpha
  min_cvm <- min(cv_out$cvm)
  opt_lambda <- cv_out$lambda[which.min(cv_out$cvm)]
  
  results_z <- rbind(results_z, data.frame(
    alpha = a,
    lambda.min = opt_lambda,
    cvm.min = min_cvm
  ))
}


# Plot of the Partial Likelihood Deviance according to Alpha
plot(results_z$alpha, results_z$cvm.min, 
     type = "b", pch = 19, col = "red",
     xlab = "Alpha (0 = Ridge, 1 = Lasso)", 
     ylab = "Partial Likelihood Deviance (CV Error)", 
     main = "Curve to select optimal alpha",
     ylim = c(4.5,4.7))


##For the final selection of the Elastic Net model, alpha=0.1 was chosen. 
# At parity of Cross-Validation Error compared to alpha=0.05,  alpha 0.1 is associated with a lower lambda value.
# This reduces the shrinkage bias on the non-zero coefficients while selecting only 8 predictive variables instead of the 13 identified at alpha = 0.05. 
# This 38% reduction in model complexity, achieved without any loss in predictive accuracy, ensures greater interpretability and superior robustness against overfitting.





###--- MODEL FOR COEFFICIENTS CALCULATION WITH FIXED ALPHA VALUES

set.seed(123)

# Cross-validation for Elastic Net (alpha 0.1) 
cvfit_FINAL <- cv.glmnet(x, y, family = "cox", alpha = 0.1) 
plot(cvfit_FINAL)

# Best lambda
best_lambda <- cvfit_FINAL$lambda.min #0.1575089

# Final model with Elastic Net
model_FINAL <- glmnet(x, y, family = "cox", alpha = 0.1, lambda = best_lambda)

       #   Df %Dev Lambda
       #   8  4.77 0.1575

# Model coefficients
c_FINAL <- coef(cvfit_FINAL, s = best_lambda)

coef_list<-as.matrix(c_FINAL)


# Vector with coefficients
coef_vector <- as.numeric(c_FINAL)

names(coef_vector) <- rownames(c_FINAL)

coef_vector<- coef_vector[coef_vector != 0]

# Save the coefficient calculated for each gene
# This object will be used to weight the previously scaled gene expression value of each new patient

save(coef_vector, file = "coef_vector_final.RData")




###--- SCORE CALCULATION

# Select the 8 genes included in the score from train_z_matrix

df_score<- train_z_matrix[,names(coef_vector)]

all(colnames(df_score) == names(coef_vector))  ###TRUE


# Calculate the gene values weighted by the coefficients

df_score_weight <- sweep(df_score, 2, coef_vector, `*`)


# Calculate the linear score for each patient by summing the weighted gene expression values

linear_score <- rowSums(df_score_weight)


# Add the score to the clinical data in order to associate it with PFS 

all(data_score$Sample_ID == Metadata$Sample_ID) ###TRUE

Metadata$Score_8g<-linear_score

# Invert the score (multiply by -1) to align higher values with a protective clinical effect

Metadata$Score_8g_inv<- -1*(Metadata$Score_8g)



###--- To apply the score to a new patients you need
#      Gene expression values of the 8 selected genes ("TNFRSF17", "CD180", "FCER2", "CD19", "CD79B", "PAX5", "MEF2C", "HLA-DOB")
#      "scaling_params_train.RData" to scale normalized gene expression data
#      "coef_vector_final.RData" to weight the previously scaled gene expression value


###--- For the validation set we used:
#      As input dataset, the matrix deposited in Array Express at the accession number E-MTAB-15877
#      Clinical data obtained from the Metadata Template deposited in Array Express at the accession number E-MTAB-15877

