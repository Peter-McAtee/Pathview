
#Explicitly set proxy for PFR
Sys.setenv("http_proxy"="http://proxy.pfr.co.nz:8080")


# Download source for Bioconductor if it is not already present
#source("https://bioconductor.org/biocLite.R")
#biocLite("pathview")

#Load the library
library(pathview)

##Set directory
myMainDir <- "C:/Users/hrapym/Desktop/Pathview_Demo"

setwd(myMainDir)
dir()

#Set and make output directories
myCSV_OUTPUT_dir <- paste0(myMainDir, "/CSV_OUTPUTS/")
myGraph_INPUT_dir <- paste0(myMainDir, "/Graph_INPUT/")

dir.create(myCSV_OUTPUT_dir)
dir.create(myGraph_INPUT_dir)


#### DATA MANIPULATION ####


myDataInput <- read.csv("Dwarfing_Demo_RPKM.csv",header=TRUE)

myDataInput <- myDataInput[1:65535,]

head(myDataInput)
dim(myDataInput)

myKeggID <- myDataInput[,2]
head(myKeggID)
length(myKeggID)

# Mold Gene data into more applicable format by subsetting data values and/or applying mathematical operations

myDataValues <- myDataInput[,3:5] 

head(myDataValues)
dim(myDataValues)

##Calculate ratio expression in log space

M9vM27 <- log(myDataValues[,3]/myDataValues[,1],2)
M9VM793 <- log(myDataValues[,3]/myDataValues[,2],2)


#Attach KEGGID to data
myPlotValues <- cbind(M9vM27, M9VM793)
head(myPlotValues)
dim(myPlotValues)

#Remove duplicates as these cause downstream problems with Pathview
#myPlotData <- subset(myPlotValues, !duplicated(myPlotValues[,1], select=myPlotValues[,1:2]))

myPlotData <- myPlotValues[!duplicated(myKeggID),]
myKeggIDUniq <- myKeggID[!duplicated(myKeggID)]

# Attached KEGG ID's as rownames . This is important for the plot function to work
rownames(myPlotData) <- myKeggIDUniq

head(myPlotData)
dim(myPlotData)

##Clean up ratio by flooring data

myPlotData[myPlotData  =="NaN"] <-0
myPlotData[myPlotData  =="NaN"] <-0
myPlotData[myPlotData =="-Inf"] <- -8.1
myPlotData[myPlotData =="Inf"] <-8.1

myDataMatrix <- as.matrix(myPlotData)
head(myDataMatrix)



### PART 1 ###


##Obtain metabolic pathways of interest
#PNG and XML files are downloaded into working directory 

myMapList <- c("00130","00900","00902","00904","00909","04075","00261","00350","00360","00400",
               "00790","00940","00945","00950","00909","00941","00944","03060","03050","03040",
               "03020","03018","03015","03013")

myMapList <- unique(myMapList)


for(i in 1:length(myMapList)){
  
  download.kegg(pathway.id = paste(myMapList[i]), species = "mdm", kegg.dir=myGraph_INPUT_dir)
}  



### PART 2 ###

## Raster expression values to the graphs
 
# Logic switch to Write Expression table for each pathway. 1 = ON, Zero = OFF
myLogicSwitch = 1

myNodeSum = "sum" # Commonly used values "sum", "mean", "max", "max.abs"

  for(i in 1:length(myMapList)){
    
      pv.OUT <- NULL
    
      pv.OUT<-  pathview(gene.data=myDataMatrix,
                pathway.id=paste(myMapList[i]),
                kegg.dir=myGraph_INPUT_dir,
                species ="mdm",
                gene.idtype ="entrez",
                out.suffix = paste("Dwarfing_EXP_",myNodeSum), 
                kegg.native = TRUE,
                multi.state = TRUE,
                node.sum = myNodeSum,
                res=600,
                bins=c(12,8),
                limit=c(-8,8))
    
    
      if(myLogicSwitch == 1){
        
        myExpresTable <- NULL
        myExpresTable <- pv.OUT$plot.data.gene
        write.csv(myExpresTable, file=paste0(myCSV_OUTPUT_dir ,"Expression_table_for_metabolicPath-",myMapList[i],".csv")) # Writes the input values that are subset into pathway tables
        
      }
    
  }

