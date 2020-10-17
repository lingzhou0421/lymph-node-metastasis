
setwd("D:\\work\\Data Analysis related\\Zhejiang_Provincial_Cancer_Hospital\\zhou_lin\\Third_甲减与淋巴结分布规律\\Datasets")

HT <- read.csv("HT.csv", header = TRUE)  # import data
head(HT)
tail(HT)
str(HT) # check data structure

varsInNeed <- c("HT", "Age", "PreTSH", "FollowUp", "Sex", "Tstage","Nstage","Mstage","ClinicalStage","NeoadjChemo",
                "Volume","NumberAll", "MaxdiameterAll", "Distance.thyroid", 
                "N.Ib","N.IIa","N.IIb", "N.III","N.IVa", "N.IVb", "N.Va", "N.Vb", "N.Vc", "N.VIIa", 
                "Maxdiameter.Ib", "Maxdiameter.IIa", "Maxdiameter.IIb","Maxdiameter.III","Maxdiameter.IVa",
                "Maxdiameter.IVb",  "Maxdiameter.Va","Maxdiameter.Vb", "Maxdiameter.Vc", "Maxdiameter.VIIa") 
# Variables of interest in this study
HT.Copy <- HT[varsInNeed] # Copy original dataset and only keep those variables of interest

nvars <- varsInNeed[c(2:4,11:13)] # numeric variabels
HT.Copy[nvars] <- lapply(HT.Copy[nvars], as.numeric)
fvars <- varsInNeed[-c(2:4,11:13)] # factor variabels
HT.Copy[fvars] <- lapply(HT.Copy[fvars], factor) 


#### step 2: Data exploration ####

#### 2.1 Demographical characteristics of the dataset (i.e. Table 1 in the study) ####
# 2.1.1 Normality distribution check for numeric variables
var.nonnml <- as.character()
for (i in 1:length(nvars)) {
  if (length(na.omit(HT.Copy[,nvars[i]])) > 3) { 
  nml.test <- shapiro.test(HT.Copy[,which(colnames(HT.Copy)==nvars[i])])
  p <- round(nml.test$p.value,3)
    if(p > 0.1) { print(paste(nvars[i],"=",p,"NormalDistribution")) }
    else {  
      print(paste(nvars[i],"=",p,"Non-NormalDistribution")) 
      var.nonnml[i] <- nvars[i] 
    }
  }
}
var.nonnml <- as.character(na.omit(var.nonnml)) # Results shows: only age is normal distributed

# 2.1.2 chisq/Fisher exact check for categorical variables
chi <- list()
j <- 0
vars.fisher <- as.character()
`%notin%` <- Negate(`%in%`)
for (i in 1:(length(fvars)-1)) {
  chi[[i]] <- chisq.test(x=HT.Copy[,fvars[i+1]],y=HT.Copy$HT,correct = FALSE) # only fvars[2:6] need to be check
  
  if ( length(HT.Copy$HT)>= 40L && ("FALSE" %notin% as.factor(chi[[i]]$expected >= 5L)) ) next  
  else if (length(HT.Copy$HT)>= 40L && ("TRUE" %in% as.factor(chi[[i]]$expected < 5L)) && 
           ("FALSE" %notin% as.factor(chi[[i]]$expected >= 1L)) ) {
    chi[[i]] <- chisq.test(x=HT.Copy[,fvars[i+1]],y=HT.Copy$HT,correct = TRUE) # chi-square correction
    print(paste(fvars[i+1],"___method=Chisq Correction",round(chi[[i]]$p.value,3)))
  }
  else if ( length(HT.Copy$HT)< 40L || ("TRUE" %in% as.factor(chi[[i]]$expected < 1L)) ) {
    chi[[i]] <- fisher.test(x=HT.Copy[,fvars[i+1]],y=HT.Copy$HT) # Fisher
    print(paste(fvars[i+1],"___method=Fisher Exact Probability",round(chi[[i]]$p.value,3)))
    j <- j+1
    vars.fisher[j] <- fvars[i+1]
  }
} 
## 2.1.3 Table 1: The general clinical characteristics of 206 patients with nasopharyngeal carcinoma
library(tableone)
table1 <- CreateTableOne(vars=c(fvars[-1],nvars), strata = c("HT"),  
                         includeNA=TRUE, addOverall=TRUE, data=HT.Copy)
summary(table1)
table1 <- print(table1, nonnormal=var.nonnml,   
                  exact=vars.fisher, showAllLevels = TRUE)
write.csv(table1,file = "D:\\work\\Data Analysis related\\Zhejiang_Provincial_Cancer_Hospital\\zhou_lin\\Third_甲减与淋巴结分布规律\\Analysis\\Results\\table1.csv")

#### 2.2  Correlation Matrix and heatmap(i.e. Figure 1 in the study)####
# try to measure the relationship between maximum diameters and the number of lymphatic nodes 
library(corrplot)
names(HT.Copy)

rltnvars <- c("N.Ib","Maxdiameter.Ib","N.IIa","Maxdiameter.IIa","N.IIb","Maxdiameter.IIb",
  "N.III","Maxdiameter.III","N.IVa","Maxdiameter.IVa","N.IVb","Maxdiameter.IVb",                       
  "N.Va","Maxdiameter.Va", "N.Vb", "Maxdiameter.Vb", "N.Vc","Maxdiameter.Vc", "N.VIIa","Maxdiameter.VIIa")             
corMatrix <- cor(as.data.frame(lapply(HT.Copy[,rltnvars], as.numeric))) 
c <- c("Ib","IIa","IIb","III","IVa","IVb","Va","Vb","Vc","VIIa")
for (i in 1:10) {
  rownames(corMatrix)[2*i-1] <- paste("No. of PLN in",c[i],sep = " ")
  rownames(corMatrix)[2*i] <- paste("Max. diameter in",c[i],sep = " ")
  colnames(corMatrix)[2*i-1] <- paste("No. of PLN in",c[i],sep = " ")
  colnames(corMatrix)[2*i] <- paste("Max. diameter in",c[i],sep = " ")
}
Fig1 <- corrplot.mixed(corr=corMatrix,lower = "number", upper = "circle",tl.pos = "lt",
               tl.col="black", tl.srt = 45)

#### step 3: Model construction (i.e Table 2 in the study)####

## 3.1 univariate analysis
allvars <- c(fvars[-1],nvars[-3]) # exclude "HT" (the outcome variable) and "followup"

# cyclical univariate logistic regression
beta <- numeric()
OR <- as.character()
pvalue <- numeric()
varname <- as.character()
j <- 0
for (i in 1:length(allvars)) {
  fmla <- as.formula(paste("HT~",paste("`",allvars[i],"`",sep = "")))
  fit <- glm(fmla, family = binomial(link = "logit"), data=HT.Copy)
  smry <- summary(fit)
    if ( (allvars[i] %in% nvars) || (allvars[i] %in% fvars &&  (length(levels(HT.Copy[[allvars[i]]])) ==2)) ) {
    # in this if condition, we try to find numeric variables or categorical variables with only two levels
      j <- j+1
      varname[j] <- rownames(smry$coefficients)[2]
      beta[j] <- sprintf("%.2f", smry$coefficients[2,1]) # beta value
      OR[j] <- paste(sprintf("%.2f",exp(smry$coefficients[2,1])),"(", sprintf("%.2f",exp(smry$coefficients[2,1]-smry$coefficients[2,2]*1.96)),"-",
                     sprintf("%.2f",exp(smry$coefficients[2,1]+smry$coefficients[2,2]*1.96)),")",sep = "") # Odds Ratio
      pvalue[j] <- sprintf("%.3f", smry$coefficients[2,4]) # p-value
    }
    else if ( allvars[i] %in% fvars &&  (length(levels(HT.Copy[[allvars[i]]])) > 2) ) {
    # in this if condition, we try to find categorical variables with more than two levels
      for (m in 1:(length(levels(HT.Copy[[allvars[i]]]))-1)) {
        j <- j+1
        varname[j] <- rownames(smry$coefficients)[m+1]
        beta[j] <- sprintf("%.2f", smry$coefficients[m+1,1])
        OR[j] <- paste(sprintf("%.2f",exp(smry$coefficients[m+1,1])),"(", sprintf("%.2f",exp(smry$coefficients[m+1,1]-smry$coefficients[m+1,2]*1.96)),"-",
                       sprintf("%.2f",exp(smry$coefficients[m+1,1]+smry$coefficients[m+1,2]*1.96)),")",sep = "")
        pvalue[j] <- sprintf("%.3f", smry$coefficients[m+1,4])
      }
    } 
}
d1 <- data.frame(varname,beta,OR,pvalue)
colnames(d1) <- c("variables", "beta", "OR(95%CI)", "p value")
write.csv(d1,file = "D:\\work\\Data Analysis related\\Zhejiang_Provincial_Cancer_Hospital\\zhou_lin\\Third_甲减与淋巴结分布规律\\Analysis\\Results\\univariate analysis.csv")

## 3.2 Multivariate analysis with differenct models 
## 3.2.1 Model 1: only clinical indicators included 
fit <- glm(HT~ Sex + Nstage + PreTSH + Volume, family = binomial(link = "logit"), data=HT.Copy)
smry <- summary(fit)  # all clinical indicators first
model.aic <- step(fit, scale = 0, direction = c("both"),trace = 0)
mdlAIC <- summary(model.aic) # selected the best model using AIC criteria
n=nrow(HT.Copy)
model.bic <- step(fit, k=log(n), direction = c("both"),trace=0)
mdlBIC <- summary(model.bic) # selected the best model using BIC criteria
model1 <- glm(mdlBIC$call$formula, family = binomial(link = "logit"), data=HT.Copy)

## 3.2.2 Model 2: clinical+lympathic nodes related factors (diameter)
HT.Copy$Distance.thyroid1 <- as.numeric(HT.Copy$Distance.thyroid)
fit2 <- glm(HT~ PreTSH + Volume + Distance.thyroid1 + 
              Maxdiameter.Ib + Maxdiameter.IIb + Maxdiameter.III + Maxdiameter.IVa + Maxdiameter.IVb + 
              Maxdiameter.Va + Maxdiameter.Vb + Maxdiameter.Vc + Maxdiameter.VIIa + 
              MaxdiameterAll, family = binomial(link = "logit"), data=HT.Copy)
model.aic <- step(fit2, scale = 0, direction = c("both"),trace = 0)
mdlAIC <- summary(model.aic)
model.bic <- step(fit2, k=log(n), direction = c("both"),trace=0)
mdlBIC <- summary(model.bic)

model2 <- glm(mdlBIC$call$formula, family = binomial(link = "logit"), data=HT.Copy)

# 3.2.3 Model 3: clinical+lympathic nodes related factors (numbers of lympatic nodes  model)
fit3 <- glm(HT~ PreTSH + Volume + Distance.thyroid1 +
              N.Ib + N.IIa + N.IIb + N.III + N.IVa + N.IVb + N.Va + 
              N.Vb + N.Vc + N.VIIa + NumberAll, family = binomial(link = "logit"), data=HT.Copy)
model.aic <- step(fit3, scale = 0, direction = c("both"),trace = 0)
mdlAIC <- summary(model.aic)
model.bic <- step(fit3, k=log(n), direction = c("both"),trace=0)
mdlBIC <- summary(model.bic)
model3 <- glm(mdlBIC$call$formula, family = binomial(link = "logit"), data=HT.Copy)
# model3 <- glm(HT ~ PreTSH + Volume + Distance.thyroid1 + N.IIb , family = binomial(link = "logit"), data=HT.Copy)

# summary
model1
model2
model3

k <- length(summary(model1)$coefficients[,1]) # change model names three times and repeat the following circulation three times
beta <- numeric()
OR <- as.character()
pvalue <- numeric()
varname <- as.character()
for (i in 1:k) {
    smry <- summary(model1) # 
    varname[i] <- rownames(smry$coefficients)[i]
    beta[i] <- sprintf("%.2f", smry$coefficients[i,1])
    OR[i] <- paste(sprintf("%.2f",exp(smry$coefficients[i,1])),"(", sprintf("%.2f",exp(smry$coefficients[i,1]-smry$coefficients[i,2]*1.96)),"-",
                   sprintf("%.2f",exp(smry$coefficients[i,1]+smry$coefficients[i,2]*1.96)),")",sep = "")
    pvalue[i] <- sprintf("%.3f", smry$coefficients[i,4])
}
d1 <- data.frame(varname,beta,OR,pvalue)
colnames(d1) <- c("variables", "beta", "OR(95%CI)", "p value")
write.csv(d1,file = "D:\\work\\Data Analysis related\\Zhejiang_Provincial_Cancer_Hospital\\zhou_lin\\Third_甲减与淋巴结分布规律\\Analysis\\Results\\multivariate1.csv")
write.csv(d1,file = "D:\\work\\Data Analysis related\\Zhejiang_Provincial_Cancer_Hospital\\zhou_lin\\Third_甲减与淋巴结分布规律\\Analysis\\Results\\multivariate2.csv")
write.csv(d1,file = "D:\\work\\Data Analysis related\\Zhejiang_Provincial_Cancer_Hospital\\zhou_lin\\Third_甲减与淋巴结分布规律\\Analysis\\Results\\multivariate3.csv")


#### Step 4: Model comparison (i.e. Figure 2 in the study) ####

library(pROC)
library(ggplot2)
## K-fold cross validation
k <- 10

pre0 <- as.numeric()
ori0 <- as.numeric()
pre1 <- as.numeric()
ori1 <- as.numeric()
pre2 <- as.numeric()
ori2 <- as.numeric()

train_pre0 <- as.numeric()
train_pre1 <- as.numeric()
train_pre2 <- as.numeric()
train_ori0 <- as.numeric()
train_ori1 <- as.numeric()
train_ori2 <- as.numeric()

auc0<-as.numeric()
auc1<-as.numeric()
auc2<-as.numeric()

train_auc0 <- as.numeric()
train_auc1 <- as.numeric()
train_auc2 <- as.numeric()

set.seed(123)

for(i in 1:k)
{
  # Train-test splitting
  # 90% of samples -> fitting
  # 10% of samples -> testing
  smp_size <- floor(0.90 * nrow(HT.Copy))
  index <- sample(seq_len(nrow(HT.Copy)),size=smp_size)
  train <- HT.Copy[index, ]
  test <- HT.Copy[-index, ]
  
  # Fitting glm
  model1 <- glm(HT~ Nstage + PreTSH + Volume, family = binomial(link = "logit"), data=train)
  model2 <- glm(HT~ PreTSH + Volume + Distance.thyroid1, family = binomial(link = "logit"), data=train)
  model3 <- glm(HT~ PreTSH + Volume + Distance.thyroid1 + N.IIb, family = binomial(link = "logit"), data=train)
  
  # Predict results in dataset of "test"
  results_prob0 <- predict(model1, newdata=test, type='response')
  results_prob1 <- predict(model2, newdata=test, type='response')
  results_prob2 <- predict(model3, newdata=test, type='response')
  
  pre0 <- append(pre0,as.numeric(results_prob0))
  ori0 <- append(ori0,as.numeric(test$HT))
  
  pre1 <- append(pre1,as.numeric(results_prob1))
  ori1 <- append(ori1,as.numeric(test$HT))
  
  pre2 <- append(pre2,as.numeric(results_prob2))
  ori2 <- append(ori2,as.numeric(test$HT))  
  
  # Predict results in dataset of "train"  
  train_prob0 <- predict(model1, type='response')
  train_prob1 <- predict(model2, type='response')
  train_prob2 <- predict(model3, type='response')
  train_pre0 <- append(train_pre0,as.numeric(train_prob0))
  train_pre1 <- append(train_pre1,as.numeric(train_prob1))
  train_pre2 <- append(train_pre2,as.numeric(train_prob2))
  train_ori0 <- append(train_ori0,as.numeric(train$HT))
  train_ori1 <- append(train_ori1,as.numeric(train$HT)) 
  train_ori2 <- append(train_ori2,as.numeric(train$HT)) 
  
  # AUC
  auc0 <- append(auc0,as.numeric(auc(as.numeric(test$HT),results_prob0)))
  auc1 <- append(auc1,as.numeric(auc(as.numeric(test$HT),results_prob1)))
  auc2 <- append(auc2,as.numeric(auc(as.numeric(test$HT),results_prob2)))
  train_auc0 <- append(train_auc0,as.numeric(auc(as.numeric(train$HT),train_prob0)))
  train_auc1 <- append(train_auc1,as.numeric(auc(as.numeric(train$HT),train_prob1)))
  train_auc2 <- append(train_auc2,as.numeric(auc(as.numeric(train$HT),train_prob2)))
  
  # ROC
  roc0 <- roc(as.numeric(test$HT),results_prob0)
  roc1 <- roc(as.numeric(test$HT),results_prob1)
  roc2 <- roc(as.numeric(test$HT),results_prob2)
}

# AUC
mean(auc0)
mean(auc1)
mean(auc2)
mean(train_auc0)
mean(train_auc1)
mean(train_auc2)

library(plotROC)

# Boxplot of AUC
Boxplot0 <- data.frame(c="Model 1",auc0)
Boxplot1 <- data.frame(c="Model 2",auc1)
Boxplot2 <- data.frame(c="Model 3",auc2)
colnames(Boxplot0) <- c("Model", "AUC")
colnames(Boxplot1) <- c("Model", "AUC")
colnames(Boxplot2) <- c("Model", "AUC")
Boxplot <-rbind(Boxplot0,Boxplot1,Boxplot2) 

f1 <- ggplot(Boxplot, aes(x=Model, y=AUC, color=Model)) + 
  geom_boxplot() +
  scale_color_manual(values=c("#FFCC33","#3366FF", "#FF3366")) +
  theme(legend.position = "bottom") +
  labs(x="Model", y = "AUC") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.75) +
  theme(panel.grid.major=element_line(colour=NA)) +
  theme(panel.grid.minor =element_line(colour=NA)) +
  theme(legend.position="none") # Remove legend

f2 <- ggplot(Boxplot, aes(x=Model, y=AUC, color=Model), width = c(0.25, 2)) + 
  geom_boxplot() + 
  coord_flip() +
  scale_color_manual(values=c("#FFCC33","#3366FF", "#9933FF")) +
  theme(legend.position = "bottom") +
  labs( y = "AUC") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.35) +
  theme_bw() + 
  theme(panel.grid.major=element_line(colour=NA)) +
  theme(panel.grid.minor =element_line(colour=NA)) +
  theme(legend.position="none") + 
  annotate("text", x = 3.45, y = 0.455, label = paste("(a)"))

# plot ROC
# plot ROC for training dataset
length(train_ori0)
ROCplot <-data.frame(c=1:1850, train_ori0)
ROCplot0 <-data.frame(c=1:1850, train_pre0)
ROCplot1 <-data.frame(c=1:1850, train_pre1)
ROCplot2 <-data.frame(c=1:1850, train_pre2)

ROCplot <- merge(ROCplot, ROCplot0, by="c")
ROCplot <- merge(ROCplot, ROCplot1, by="c")
ROCplot <- merge(ROCplot, ROCplot2, by="c")

longtest <- melt_roc(ROCplot, "train_ori0", c("train_pre0", "train_pre1", "train_pre2")) 
head(longtest)
longtest$name <- ifelse(longtest$name=="train_pre0", "Model 1", 
                        ifelse(longtest$name=="train_pre1","Model 2", "Model 3"))
colnames(longtest) <- c("D", "M", "Model")

f3 <- ggplot(longtest, aes(d = D, m = M, color = Model)) + 
  geom_roc(n.cuts = 0) +
  scale_color_manual(values=c("#FFCC33","#3366FF", "#9933FF"), name=NULL,
                     labels=c(paste("Model 1：AUC =",sprintf("%.3f",mean(train_auc0))),
                              paste("Model 2：AUC =",sprintf("%.3f",mean(train_auc1))),
                              paste("Model 3：AUC =",sprintf("%.3f",mean(train_auc2))))) +
  style_roc( xlab = "1 - Specificity", ylab = "Sensitivity")  + 
  annotate("text", x = .10, y = 0.98, label = paste("(b) training set")) +
  theme(legend.background = element_blank()) +
  theme(legend.justification=c(0.5,0.5), legend.position=c(0.75,0.15))


# plot ROC for validation dataset
length(ori0)
ROCplot <-data.frame(c=1:210, ori0)
ROCplot0 <-data.frame(c=1:210, pre0)
ROCplot1 <-data.frame(c=1:210, pre1)
ROCplot2 <-data.frame(c=1:210, pre2)

ROCplot <- merge(ROCplot, ROCplot0, by="c")
ROCplot <- merge(ROCplot, ROCplot1, by="c")
ROCplot <- merge(ROCplot, ROCplot2, by="c")

longtest <- melt_roc(ROCplot, "ori0", c("pre0", "pre1", "pre2")) 
head(longtest)
longtest$name <- ifelse(longtest$name=="pre0", "Model 1",
                 ifelse(longtest$name=="pre1","Model 2", "Model 3"))
colnames(longtest) <- c("D", "M", "Model")

f4 <- ggplot(longtest, aes(d = D, m = M, color = Model)) + 
  geom_roc(n.cuts = 0) +
  scale_color_manual(values=c("#FFCC33","#3366FF", "#9933FF"), name=NULL,
                     labels=c(paste("Model 1：AUC =",sprintf("%.3f",mean(auc0))),
                              paste("Model 2：AUC =",sprintf("%.3f",mean(auc1))),
                              paste("Model 3：AUC =",sprintf("%.3f",mean(auc2))))) +
  style_roc( xlab = "1 - Specificity", ylab = "Sensitivity") + 
  annotate("text", x = .135, y = 0.98, label = paste("(c) validation set")) + 
  theme(legend.background = element_blank()) +
  theme(legend.justification=c(0.5,0.5), legend.position=c(0.75,0.15))

# layout
library(customLayout)
lay1 <- lay_new(mat = matrix(1:2, ncol = 2), widths = c(1,1))    
lay2 <- lay_new(mat = matrix(1, ncol = 1), widths = c(2))
cl_1 <-lay_bind_row(lay2,lay1,  heights = c(6, 12), addmax = TRUE)
lay_grid(list(f2,f3,f4), cl_1)  # Figure 2 in the study


#### Step 5: HT prediction (i.e. figure 3 in the study) ####
require(rms)
finalvars <- c("Nstage","PreTSH","Volume", "Distance.thyroid1")
HT1 <- HT.Copy[,c(1,which(names(HT.Copy) %in% finalvars))]
colnames(HT1) <- c("HT", "Pretreatment TSH concentration", 
               "N.stage", "Thyroid volume", "Distance from thyroid")

dd <- datadist(HT1) # format data for Nomogram
options(datadist = "dd")

model2 <- lrm(HT~`Pretreatment TSH concentration`+ `Thyroid volume` +`Distance from thyroid`, data=HT1)
nomo <- nomogram(model2, fun=plogis,lp=F,
                  fun.at=c(.001,.01,.05,seq(.1,.9,by=.1),.95,.99,.999),
                  vnames=c("names"),
                  funlabel="Risk of Hypothyroidism")
plot(nomo,xfrac=.4, cex.axis=.8, cex.var=0.85,tcl=-0.25, lmgp=.2,
     points.label='Points', total.points.label='Total Points') # figure 3 in the study
