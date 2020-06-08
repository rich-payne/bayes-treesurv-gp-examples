################
## BC-EXAMPLE ##
################

The BC-EXAMPLE is a small subset of gene expression profiles originated from the breast cancer data set of Van't Veer (2002) (available at 
http://www.rii.com/publications/2002/vantveer.html). We sampled 50 samples from the total 295 breast cancer patients and applied the kmeans clustering
algorithm to reduce the number of covariates to 300 metagenes. Further, we randomly split the data in training (2/3 of the samples) and test sets (1/3).


*Files:

X.train 
Explanatory variables for training (in sample) data. Matrix with 33 rows corresponding to samples and 500 columns to covariates. 

X.test 
Explanatory variables for test (out of sample) data. Matrix with 17 rows corresponding to samples and 500 columns to covariates. 
 
Y.train 
Response variables for training (in sample) data. Matrix with 33 rows corresponding to samples and 2 columns: "Survival" and "Status". 
"Survival" lists the natural log of time-to-event and "Status" lists the indicator variable for event occurrence, i.e., 1 if event occurred and 0 if censored. 

Y.test
Response variables for test (out of sample) data. Matrix with 17 rows corresponding to samples and 2 columns: "Survival" and "Status". 
"Survival" lists the natural log of time-to-event and "Status" lists the indicator variable for event occurrence, i.e., 1 if event occurred and 0 if censored.
Not needed for running the functions, however might be useful for assessing the fit.


*References
Van't Veer,L.J. et al. (2002) Gene expression profiling predicts clinical outcome of breast cancer. Nature,415,530-536.


*Example
#loading data
X.train<-read.table("X.train.txt")
X.test<-read.table("X.test.txt")
Y.train<-read.table("Y.train.txt")
Y.test<-read.table("Y.test.txt")
