#Emotional overeating in Gemini twin cohort
#Authors: Vaishnavi K. Madhava, Zeynep Nas, and Moritz Herle
#INITIAL DATA PREPERATION

##### residualise the EOE scores 

data$rEOE15m <- residuals (lm(CEBQ_EOE_15mths ~ ageT1_yrs + sexT0, data = data, na.action = "na.exclude"))

data$rEOE5y <- residuals (lm(CEBQ_EOE_5yrs ~ ageT7 + sexT0, data = data, na.action = "na.exclude"))

data$rEOE12y <- residuals (lm(CEBQ_EOE_T9 ~ ageT9 + sexT0, data = data, na.action = "na.exclude"))


#### Shifting the distribution of variables above 0 ("+1") and increasing the SD by multiplying by 0.5
### As the variation is quite small, sometimes increasing the variance can make the model a bit better
data$rEOE15m <- ((data$rEOE15m +1)*0.5)
data$rEOE5y <- ((data$rEOE5y +1)*0.5)
data$rEOE12y <- ((data$rEOE12y +1)*0.5)



#########################################################################################
# Multivariate Cholesky ACE
#####################################



#creating wide dataset
WideTwinData  <- reshape (data, idvar = c("famID", "Zygosity_new"), timevar = "Random", direction = "wide")
describe(WideTwinData)
View(WideTwinData)

#Selecting only the necessary variables
EOEData         <- WideTwinData[,c("Zygosity_new", "rEOE15m.1","rEOE5y.1","rEOE12y.1", "rEOE15m.2", "rEOE5y.2",  "rEOE12y.2")]

describe(EOEData)

#to rename the variable names (each one corresponding to the column names of the dataset)
names(EOEData) <- c("zyg", "rEOE15m1", "rEOE5y1", "rEOE12y1", "rEOE15m2", "rEOE5y2",  "rEOE12y2")

describe(EOEData)

View(EOEData)


# Select variables for analysis
nv        <- 3       # number of traits
ntv       <- nv*2    # number of total variables

nlower	<- ntv*(ntv+1)/2 	# number of free elements in a lower matrix ntv*ntv
ncor	<-(nv*(nv+1)/2)-nv	# number of free elements in a corelation matrix nv*nv

Vars      <- c("EOE15", "EOE5", "EOE12")
selVars   <- c("rEOE15m1", "rEOE5y1", "rEOE12y1","rEOE15m2", "rEOE5y2", "rEOE12y2")

mzData <- subset(EOEData, zyg==1, selVars)
dzData <- subset(EOEData, zyg==2, selVars)  

# To create Labels for Lower Triangular Matrices
aLabs <- paste("a", do.call(c, sapply(seq(1, nv), function(x){ paste(x:nv, x,sep="") })), sep="")
cLabs <- paste("c", do.call(c, sapply(seq(1, nv), function(x){ paste(x:nv, x,sep="") })), sep="")
eLabs <- paste("e", do.call(c, sapply(seq(1, nv), function(x){ paste(x:nv, x,sep="") })), sep="")


mLabs	<- paste("m",1:nv,sep="")

# To create start values 
Stmean <-c(0,0,0,0,0,0)
svSD   <-c(1,1,1,1,1,1)
svMZ <-c(0.5,
         0.5, 0.5, 
         0.5, 0.5, 0.5,
         0.5, 0.5, 0.5, 0.5, 
         0.5, 0.5, 0.5, 0.5, 0.5,
         0.5, 0.5, 0.5, 0.5, 0.5, 0.5)

svDZ <-c(0.2,
         0.2, 0.2, 
         0.2, 0.2, 0.2,
         0.2, 0.2, 0.2, 0.2, 
         0.2, 0.2, 0.2, 0.2, 0.2,
         0.2, 0.2, 0.2, 0.2, 0.2, 0.2)


# ----------------------------------------------------------------------------------------
# 1 Specify Fully Saturated Model (Phenotypic Cholesky Decomposition)
# ----------------------------------------------------------------------------------------

LabsMMz <- c("MZm11", "MZm21", "MZm31", "MZm12", "MZm22","MZm32")
LabsMDz <- c("DZm11", "DZm21", "DZm31", "DZm12", "DZm22","DZm32")

LabsSMz <- c("MZs11", "MZs21", "MZs31", "MZs12", "MZs22", "MZs32")
LabsSDz <- c("DZs11", "DZs21", "DZs31", "DZs12", "DZs22", "DZs32")



# Matrix & Algebra for expected means and covariances 
MZlow		<-mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=T, values=5, name="LowMZ" )
CovMZ		<-mxAlgebra( expression=LowMZ %*% t(LowMZ), name="ExpCovMZ" )
MeanMZ	<-mxMatrix( type="Full", nrow=1, ncol=ntv, free=T, values=Stmean,labels= LabsMMz, name="ExpMeanMZ" )

DZlow		<-mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=T, values=5, name="LowDZ" )
CovDZ		<-mxAlgebra( expression=LowDZ %*% t(LowDZ), name="ExpCovDZ" )
MeanDZ	<-mxMatrix( type="Full", nrow=1, ncol=ntv, free=T, values=Stmean,labels= LabsMDz,  name="ExpMeanDZ" )

# Data objects for Multiple Groups
dataMZ   <- mxData( observed=mzData, type="raw" )
dataDZ   <- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ    <- mxExpectationNormal( covariance="ExpCovMZ", means="ExpMeanMZ", dimnames=selVars)
objDZ    <- mxExpectationNormal( covariance="ExpCovDZ", means="ExpMeanDZ", dimnames=selVars)

fitFunction <- mxFitFunctionML()

# Combine Groups
modelMZ	<- mxModel( MZlow, CovMZ, MeanMZ, dataMZ, objMZ, fitFunction, name="MZ")
modelDZ	<- mxModel( DZlow, CovDZ, MeanDZ, dataDZ, objDZ, fitFunction, name="DZ")
minus2ll	<- mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj		<- mxFitFunctionAlgebra( "m2LL" )
SatModel 	<- mxModel( "Sat", modelMZ, modelDZ, minus2ll, obj)

# -------------------------------------------------------------------------------------------------------------------------------
# 1a)	RUN Saturated Model (Cholesky Decomposition)

SatFit		<- mxRun(SatModel)
(SatSum		<- summary(SatFit))








# -----------------------------------------------------------------------------------------------------------------------------------
# 1a)	Specify Multivariate Saturated Model (Gaussian Decomposition) 
# 	We will use this specification to fit a constrained model
# -------------------------------------------------------------



LabsMMz <- c("MZm11", "MZm21", "MZm31", "MZm12", "MZm22","MZm32")
LabsMDz <- c("DZm11", "DZm21", "DZm31", "DZm12", "DZm22","DZm32")

LabsSMz <- c("MZs11", "MZs21", "MZs31", "MZs12", "MZs22", "MZs32")
LabsSDz <- c("DZs11", "DZs21", "DZs31", "DZs12", "DZs22", "DZs32")

LabsCorMz	<- c("MZCor12",	"MZCor13", "MZCor14", "MZCor15", "MZCor16", 
                          "MZCor23", "MZCor24", "MZCor25", "MZCor26",
                                     "MZCor34", "MZCor35", "MZCor36",
                                                "MZCor45", "MZCor46",
                                                           "MZCor56")

                                    

LabsCorDz	<- c("DZCor12",	"DZCor13", "DZCor14", "DZCor15", "DZCor16", 
                          "DZCor23", "DZCor24", "DZCor25", "DZCor26",
                                     "DZCor34", "DZCor35", "DZCor36",
                                                "DZCor45", "DZCor46",
                                                           "DZCor56")

                                    
# Matrix & Algebra for expected means and covariances 
MeanMZ	<-mxMatrix( "Full", 1, ntv, free=T, values=Stmean, labels= LabsMMz, name="ExpMeanMZ" )
MZsd		<-mxMatrix( "Diag", ntv, ntv, free=T, values=svSD, labels =LabsSMz, name="sdMZ" )
Cormz		<-mxMatrix( "Stand", ntv, ntv, free=T, values=0.5, labels=LabsCorMz, name="MZCor") 
CovMZ		<-mxAlgebra( expression=sdMZ %*% MZCor %*% t(sdMZ), name="ExpCovMZ" )

MeanDZ	<-mxMatrix( "Full", 1, ntv, free=T, values=Stmean, labels=LabsMDz, name="ExpMeanDZ" )
DZsd		<-mxMatrix( "Diag", ntv, ntv, free=T, values=svSD, labels=LabsSDz, name="sdDZ" )
Cordz		<-mxMatrix( "Stand", ntv, ntv, free=T, values=0.5, labels=LabsCorDz,  name="DZCor")
CovDZ		<-mxAlgebra( expression=sdDZ %*% DZCor %*% t(sdDZ), name="ExpCovDZ" ) 

# Data objects for Multiple Groups
dataMZ   <- mxData( observed=mzData, type="raw" )
dataDZ   <- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ    <- mxExpectationNormal( covariance="ExpCovMZ", means="ExpMeanMZ", dimnames=selVars)
objDZ    <- mxExpectationNormal( covariance="ExpCovDZ", means="ExpMeanDZ", dimnames=selVars)

fitFunction <- mxFitFunctionML()

# Combine Groups
modelMZ	<- mxModel( MeanMZ, MZsd, Cormz, CovMZ, dataMZ, objMZ, fitFunction, name="MZ")
modelDZ	<- mxModel( MeanDZ, DZsd, Cordz, CovDZ, dataDZ, objDZ, fitFunction, name="DZ")
minus2ll	<- mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj		<- mxFitFunctionAlgebra( "m2LL" )

#cross twin within trait
Conf1		<- mxCI (c ('MZ.MZCor[4,1]', 'MZ.MZCor[5,2]', 'MZ.MZCor[6,3]') )
Conf2		<- mxCI (c ('DZ.DZCor[4,1]', 'DZ.DZCor[5,2]', 'DZ.DZCor[6,3]') )

# Cross twin cross trait
Conf3		<- mxCI (c ('MZ.MZCor[5,1]', 'MZ.MZCor[6,1]', 'MZ.MZCor[6,2]') )
Conf4		<- mxCI (c ('DZ.DZCor[5,1]', 'DZ.DZCor[6,1]', 'DZ.DZCor[6,2]') )

#rPh
Conf5		<- mxCI (c ('MZ.MZCor[2,1]', 'MZ.MZCor[3,1]', 'MZ.MZCor[3,2]') )
Conf6		<- mxCI (c ('DZ.DZCor[2,1]', 'DZ.DZCor[3,1]', 'DZ.DZCor[3,2]') )

SatGModel 	<- mxModel( "SatG", modelMZ, modelDZ, minus2ll, obj, Conf1, Conf2, Conf3 ,Conf4, Conf5, Conf6)

# -------------------------------------------------------------------------------------------------------------------------------
# 1b)	RUN Saturated Model (Gaussian Decomposition)

SatGFit	<- mxTryHard(SatGModel, intervals=T)
(SatGSum	<- summary(SatGFit))


# --------------------------------------------------------------------------------------------
# 1b)	Specify & Run Suba: Constrained Model
#	Equal Means & Variances across Twin Order 
#	One overall set of Within-person cross-trait correlations
#	Symmetric xtwin-xtrait cor matrices in MZ and DZ group 
# --------------------------------------------------------------------------------------------

# To manipulate  the parameters in matrices we wish to change in the full model, we use the
# 'omxSetParameters' function with the original 'labels' and 'newlabels' to indicate the changes
# i.e. specifying the same label effectively constraints the parameters to be the same 

Sub1aModel	<- mxModel(SatGModel, name="Sub1a" )

# means

Sub1aModel	<- omxSetParameters( Sub1aModel, free=T, values=Stmean, labels=c("MZm11", "MZm21", "MZm31",  "MZm12", "MZm22","MZm32"),
                                newlabels=c("MZm1", "MZm2", "MZm3","MZm1", "MZm2", "MZm3"))

Sub1aModel	<- omxSetParameters( Sub1aModel, free=T, values=Stmean, labels=c("DZm11", "DZm21", "DZm31", "DZm12", "DZm22","DZm32"), 
                                newlabels=c("DZm1", "DZm2", "DZm3","DZm1", "DZm2", "DZm3"))

# SDs

Sub1aModel	<- omxSetParameters( Sub1aModel, free=T, values=svSD, labels=c("MZs11", "MZs21", "MZs31", "MZs12", "MZs22", "MZs32"),
                                newlabels=c("MZs1", "MZs2", "MZs3","MZs1", "MZs2", "MZs3"))

Sub1aModel	<- omxSetParameters( Sub1aModel, free=T, values=svSD, labels=c("DZs11", "DZs21", "DZs31", "DZs12", "DZs22", "DZs32"),
                                newlabels=c("DZs1", "DZs2", "DZs3","DZs1", "DZs2", "DZs3"))


# ------------------------------------------------------------------------------
# 2) RUN Sub1Model
Sub1aFit     <- mxTryHard( Sub1aModel, intervals=F )
(Sub1aSum     <- summary( Sub1aFit ))

mxCompare(SatGFit, Sub1aFit)


# --------------------------------------------------------------------------------------------
# 1b)	Submodel 1b
#	Equal Means & Variances across Twin Order & zyg group
#	One overall set of Within-person cross-trait correlations
#	Symmetric xtwin-xtrait cor matrices in MZ and DZ group 
# --------------------------------------------------------------------------------------------

# To manipulate  the parameters in matrices we wish to change in the full model, we use the
# 'omxSetParameters' function with the original 'labels' and 'newlabels' to indicate the changes
# i.e. specifying the same label effectively constraints the parameters to be the same 

Sub1bModel	<- mxModel(Sub1aModel, name="Sub1b" )

# means

Sub1bModel	<- omxSetParameters(Sub1bModel, free=T, values=Stmean, labels=c("MZm1", "MZm2", "MZm3"),
                              newlabels=c("m1", "m2", "m3"))

Sub1bModel	<- omxSetParameters( Sub1bModel, free=T, values=Stmean, labels=c("DZm1", "DZm2", "DZm3"), 
                               newlabels=c("m1", "m2", "m3"))

# SDs

Sub1bModel	<- omxSetParameters( Sub1bModel, free=T, values=svSD, labels=c("MZs1", "MZs2", "MZs3"),
                               newlabels=c("s11", "s21", "s31"))

Sub1bModel	<- omxSetParameters( Sub1bModel, free=T, values=svSD, labels=c("DZs1", "DZs2", "DZs3"),
                               newlabels=c("s11", "s21", "s31"))

# cor 

Sub1bModel	<- omxSetParameters( Sub1bModel, free=T, values=svMZ, labels= c("MZCor12",	"MZCor13", "MZCor14", "MZCor15", "MZCor16", 
                                                                         "MZCor23", "MZCor24", "MZCor25", "MZCor26",
                                                                         "MZCor34", "MZCor35", "MZCor36",
                                                                         "MZCor45", "MZCor46",
                                                                         "MZCor56"),
                               
                               
                               newlabels=c("rph12","rph13", "MZ15m", "MZxtxt_12","MZxtxt_13", 
                                           "rph23", "MZxtxt_12", "MZ5y", "MZxtxt_23",
                                           "MZxtxt_13", "MZxtxt_23", "MZ12y",
                                           "rph12", "rph13",
                                           "rph23"))
                                          



Sub1bModel	<- omxSetParameters( Sub1bModel, free=T, values=svMZ, labels=c("DZCor12",	"DZCor13", "DZCor14", "DZCor15", "DZCor16", 
                                                                        "DZCor23", "DZCor24", "DZCor25", "DZCor26",
                                                                        "DZCor34", "DZCor35", "DZCor36",
                                                                        "DZCor45", "DZCor46",
                                                                        "DZCor56"),

                               
                               newlabels=c("rph12","rph13", "DZ15m", "DZxtxt_12","DZxtxt_13", 
                                           "rph23", "DZxtxt_12", "DZ5y", "DZxtxt_23",
                                           "DZxtxt_13", "DZxtxt_23", "DZ12y",
                                           "rph12", "rph13",
                                           "rph23"))
                                           


# ------------------------------------------------------------------------------
# 2) RUN Sub1Model
Sub1bFit      <- mxTryHard(Sub1bModel,extraTries = 20, intervals=T )
(Sub1bSum     <- summary( Sub1bFit, verbose=T ))

mxCompare(SatGFit, Sub1aFit)
mxCompare(Sub1aFit, Sub1bFit)
mxCompare(SatGFit, Sub1bFit)





#########################################
#### specify model 	Cholesky ACE model 
#########################################


# Matrix & Algebra for expected means vector
# Matrices declared to store a, c, and e Path Coefficients
pathA	<- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=.4, labels=aLabs, name="a" )
pathC	<- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=.1, labels=cLabs, name="c" )
pathE	<- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=.2, labels=eLabs, name="e" )

# Matrices generated to hold A, C, and E computed Variance Components
covA	<- mxAlgebra( expression=a %*% t(a), name="A" )
covC	<- mxAlgebra( expression=c %*% t(c), name="C" ) 
covE	<- mxAlgebra( expression=e %*% t(e), name="E" )
covP	<- mxAlgebra( expression=A+C+E, name="V" )
StA	<- mxAlgebra( expression=A/V, name="h2" )
StC	<- mxAlgebra( expression=C/V, name="c2" )
StE	<- mxAlgebra( expression=E/V, name="e2" )

# Algebra to compute Correlations
matI	<- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I")
Rph	<- mxAlgebra( expression=solve(sqrt(I*V)) %&% V , name="Phcor")
Rg	<- mxAlgebra( expression=solve(sqrt(I*A)) %&% A , name="Acor")
Rc	<- mxAlgebra( expression=solve(sqrt(I*C)) %&% C , name="Ccor")
Re	<- mxAlgebra( expression=solve(sqrt(I*E)) %&% E , name="Ecor")

stv  <- mxAlgebra( expression = solve(sqrt(I*V)), name= "StV")

stCholA <- mxAlgebra( expression = cbind(StV %*% a), name = "stcholA")
stCholA2 <- mxAlgebra( expression = cbind(stcholA*stcholA), name = "stcholA2")

stCholC <- mxAlgebra( expression = cbind(StV %*% c), name = "stcholC")
stCholC2 <- mxAlgebra( expression = cbind(stcholC*stcholC), name = "stcholC2")

stCholE <- mxAlgebra( expression = cbind(StV %*% e), name = "stcholE")
stCholE2 <- mxAlgebra( expression = cbind(stcholE*stcholE), name = "stcholE2")



# Algebra for expected Mean and Variance/Covariance Matrices in MZ & DZ twins
Mean	<- mxMatrix( type="Full", nrow=1, ncol=ntv, free=T, values=Stmean, labels=c(mLabs,mLabs), name="ExpMean" )

covMZ	<- mxAlgebra( expression= rbind( cbind(A+C+E , A+C),
                                       cbind(A+C , A+C+E)),		name="ExpCovMZ" )
covDZ	<- mxAlgebra( expression= rbind( cbind(A+C+E       , 0.5%x%A+C),
                                       cbind(0.5%x%A+C , A+C+E)),	name="ExpCovDZ" )

# Data objects for Multiple Groups
dataMZ	<- mxData( observed=mzData, type="raw" )
dataDZ	<- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ	<- mxExpectationNormal( covariance="ExpCovMZ", means="ExpMean", dimnames=selVars )
objDZ	<- mxExpectationNormal( covariance="ExpCovDZ", means="ExpMean", dimnames=selVars )

fitFunction <- mxFitFunctionML()

# Combine Groups
pars	<- list( pathA, pathC, pathE, covA, covC, covE, covP, StA, StC, StE, matI, Rph, Rg, Rc, Re, stv, stCholA, stCholA2, stCholC, stCholC2, stCholE, stCholE2) 
modelMZ	<- mxModel( pars, covMZ, Mean, dataMZ, objMZ, fitFunction, name="MZ" )
modelDZ	<- mxModel( pars, covDZ, Mean, dataDZ, objDZ, fitFunction, name="DZ" )
minus2ll<- mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj	<- mxFitFunctionAlgebra( "m2LL" )
conf1	<- mxCI (c ('MZ.h2[1,1]', 'MZ.h2[2,2]', 'MZ.h2[3,3]' ) )		# h2
conf2	<- mxCI (c ('MZ.c2[1,1]', 'MZ.c2[2,2]', 'MZ.c2[3,3]') )		# c2
conf3	<- mxCI (c ('MZ.e2[1,1]', 'MZ.e2[2,2]', 'MZ.e2[3,3]') )		# e2
#conf4	<- mxCI (c ('MZ.Acor[2,1]','MZ.Acor[3,1]','MZ.Acor[3,2]' )) #Rg
#conf5	<- mxCI (c ('MZ.Ccor[2,1]','MZ.Ccor[3,1]','MZ.Ccor[3,2]' )) #Rc			
#conf6	<- mxCI (c ('MZ.Ecor[2,1]','MZ.Ecor[3,1]','MZ.Ecor[3,2]' )) #Re
conf7	<- mxCI (c ('Phcor[2,1]', 'Phcor[3,1]', 'Phcor[3,2]')) #Rph
ciCholA		<- mxCI(c('stcholA2[1,1]', 'stcholA2[2,1]', 'stcholA2[2,2]', 'stcholA2[3,1]', 'stcholA2[3,2]', 'stcholA2[3,3]'))			
ciCholC		<- mxCI(c('stcholC2[1,1]', 'stcholC2[2,1]', 'stcholC2[2,2]', 'stcholC2[3,1]', 'stcholC2[3,2]', 'stcholC2[3,3]'))			
ciCholE		<- mxCI(c('stcholE2[1,1]', 'stcholE2[2,1]', 'stcholE2[2,2]', 'stcholE2[3,1]', 'stcholE2[3,2]', 'stcholE2[3,3]'))			

CholAceModel<- mxModel( "CholACE", pars, modelMZ, modelDZ, minus2ll, obj, conf1, conf2, conf3, conf7, ciCholA, ciCholC,ciCholE )

# -------------------------------------------------
# 	RUN full ACE MODEL
CholAceFit	<- mxRun(CholAceModel, intervals=T)
(CholAceSum	<- summary(CholAceFit))

mxCompare(SatGFit   , Sub1aFit)
mxCompare(Sub1aFit   , Sub1bFit)
mxCompare(SatGFit   , CholAceFit)


