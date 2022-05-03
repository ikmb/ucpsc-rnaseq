library(caret)

prefiltereddata <- rftrain[,c(top_importance$feature[1:25],"UCpositive"),with=FALSE]
prefiltereddata$UCpositive <- as.integer(prefiltereddata$UCpositive)-1

# You can use any threshold you want to deem a correlation too high. Here we use .80
nonColinearData = prefiltereddata[, -findCorrelation(cor(prefiltereddata), cutoff = .8),with=FALSE]


nonColinearData$UCpositive <- factor(nonColinearData$UCpositive)

#rfcvmodel <- rfcv(prefiltereddata[,-"UCpositive"],as.factor(prefiltereddata$UCpositive))
#rfcvmodel$error.cv


# Set RFE control
ctrl = rfeControl(functions = rfFuncs, # "rfFuncs" are built-in to caret
                  method = "repeatedcv", repeats = 10,
                  saveDetails = TRUE)
# By using rfFuncs, caret will use a random forest to evaluate the usefulness of a feature.

# Set a sequence of feature-space sizes to search over:
sizes = seq(sqrt(ncol(nonColinearData))*.5, ncol(nonColinearData), by = 5)
# note, this will fit hundreds of forests (not trees), so it may take a while.

# Use caret's rfe function to fit RF models to these different feature spaces
rfeResults = rfe(x = select(nonColinearData, -UCpositive), y = nonColinearData$UCpositive,
                 sizes = sizes,
                 rfeControl = ctrl)
unique(rfeResults$variables$var)
