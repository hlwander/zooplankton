#########################################
#### Zoop Summary code -Density and Biomass
### Code modified from DCR --> generates summary csv files to use for fig generation

##This runs the code for zooplankton files to summarize and do calculations
#It will upload two zooplankton files from the template
#Summarize data and then output a file that can be merged with other lake data
#########################################

#######################START OF READING IN DATA#######################################
###Function to get the median if odd number in a column, or one entry below median if even number####
median.zoop<-function(vector){
  #Remove nas and then sort
  vector2<-vector[!is.na(vector)]
  vector2<-sort(vector2)
  #If the vector contains an odd number of entries, return the median
  if(length(vector2)%%2==1){return(median(vector2))}else{
    return(vector2[length(vector2)/2]) #If the vector is even #, return the entry immediately below the median
  }
}

###############################End of function

#Read in the csv file that contains all the taxa from the zooplankton key that we use:####
#http://cfb.unh.edu/cfbkey/html/index.html
ZoopTaxonomy<-read.csv("Zooplankton-TaxonomicIdentification.csv", fill = TRUE,as.is=TRUE)
#record which rows are NA
SpeciesNA.rows<-which(is.na(ZoopTaxonomy$Species)) 
#Replace the species with Genus and species, e.g., D. magna
ZoopTaxonomy$Species<-paste0(substring(ZoopTaxonomy$Genus,1,1),". ",ZoopTaxonomy$Species)
#If there are any NA, then fill those back in
ZoopTaxonomy$Species[SpeciesNA.rows]<-NA   

#Read in csv file for size to biomass conversions for rotifers with all taxa ID'ed####
RotiferBiomassConversions<-read.csv("ZooplanktonLengthBiomassConversion-Rotifers.csv", fill = TRUE,as.is=TRUE) 
#Any of the blank species get NA 
RotiferBiomassConversions$Species[RotiferBiomassConversions$Species==""]<-NA
#Replace the species with Genus and species, e.g., K. earlinae
RotiferBiomassConversions$Species[!is.na(RotiferBiomassConversions$Species)]<-paste0(substring(RotiferBiomassConversions$Genus[!is.na(RotiferBiomassConversions$Species)],1,1),". ",RotiferBiomassConversions$Species[!is.na(RotiferBiomassConversions$Species)])

#Read in the csv file with the occular lens conversions####
ZoopOccularLensConversions<-read.csv("Zooplankton-OccularLensConversions.csv", fill = TRUE,as.is=TRUE)

#Read in csv file for size to biomass conversions for crustaceans with all taxa ID'ed####
CrustaceanBiomassConversions<-read.csv("ZooplanktonLengthBiomassConversion-Crustaceans.csv", fill = TRUE,as.is=TRUE) 
#Any of the blank species get NA 
CrustaceanBiomassConversions$Species[CrustaceanBiomassConversions$Species==""]<-NA
#Replace the species with Genus and species, e.g., D. magna
CrustaceanBiomassConversions$Species[!is.na(CrustaceanBiomassConversions$Species)]<-paste0(substring(CrustaceanBiomassConversions$Genus[!is.na(CrustaceanBiomassConversions$Species)],1,1),". ",CrustaceanBiomassConversions$Species[!is.na(CrustaceanBiomassConversions$Species)])
#Move all Cladacora ids to suborder
CrustaceanBiomassConversions$Suborder[CrustaceanBiomassConversions$Order=="Cladocera"]<-"Cladocera"
CrustaceanBiomassConversions$Order[CrustaceanBiomassConversions$Order=="Cladocera"]<-""
#Make sure all Daphniidae are subOrder cladocera
CrustaceanBiomassConversions$Suborder[CrustaceanBiomassConversions$Family=="Daphniidae"]<-"Cladocera"

#For loop through all the years(or specify year here)
year<-2019

#Names of the two files for a particular year
dataFiles<-paste("FCR",year,"_ZooplanktonCounting_SizeID_DataEntry.csv",sep="")

#Create SizeID data frames for that year
SizeID<-read.csv(dataFiles, fill = TRUE,as.is=TRUE)

#If any SizeID are missing, add NA
SizeID$LowestTaxonomicLevelOfID[SizeID$LowestTaxonomicLevelOfID==""]<-NA
SizeID$TaxaID[SizeID$TaxaID==""]<-NA
SizeID$MarksInOcularMicrometer_No.[SizeID$MarksInOcularMicrometer_No.==""]<-NA
SizeID$ObjectiveMagnification[SizeID$ObjectiveMagnification==""]<-NA

#convert from character to numeric
SizeID$MarksInOcularMicrometer_No.<- as.numeric(SizeID$MarksInOcularMicrometer_No.)

#Calculate the size based on the marks in the occular micrometer and conversion
#Checks for non-NA and non-numeric strings
#Prints out the row to examine
for(i in 1:length(SizeID$MarksInOcularMicrometer_No.)){
  if(!is.na(SizeID$MarksInOcularMicrometer_No.[i])&is.na(as.numeric(SizeID$MarksInOcularMicrometer_No.[i]))){print(paste("Row",i))}
}

#Calculates the size of each zooplankton
SizeID$Size_mm<-as.numeric(SizeID$MarksInOcularMicrometer_No.)*(ZoopOccularLensConversions[match(SizeID$ObjectiveMagnification,ZoopOccularLensConversions$Objective),2]/ZoopOccularLensConversions[match(SizeID$ObjectiveMagnification,ZoopOccularLensConversions$Objective),3])

#Create a new column for each level of taxonomic richness
SizeID[,names(ZoopTaxonomy)]<-NA

#Check if order name has a first letter capital
substr(SizeID$LowestTaxonomicLevelOfID, 1, 1) <- toupper(substr(SizeID$LowestTaxonomicLevelOfID, 1, 1))

#Figure out which column the taxa ID is
for(j in 1:length(SizeID$TaxaID)){
  #Pick out only values with Taxa Level and ID
  if(!is.na(SizeID$LowestTaxonomicLevelOfID[j])&!is.na(SizeID$TaxaID[j])){
    #Get the column that the lowest level of ID
    TaxaColumn<-match(SizeID$LowestTaxonomicLevelOfID[j],names(ZoopTaxonomy))
    #Get the row for the correct taxa (it will output the first match)
    TaxaRow<-match(SizeID$TaxaID[j],ZoopTaxonomy[,TaxaColumn])
    #Run through each of the levels up to the lowest level and store it in new columns in the SizeID data frame
    #if nauplius, then record lowest taxanomic level as Copepoda (insert NA in order-species columns)
    for(k in 1:TaxaColumn){
      SizeID[j,names(ZoopTaxonomy)[k]]<-ZoopTaxonomy[TaxaRow,k]
    }
  }
  for(rot.i in 5:10){
    if(SizeID$TaxaID[j]=="Cyclopoida" | SizeID$TaxaID[j]=="Calanoida" & is.na(SizeID$Nauplius[j]=="nauplius")){
      SizeID[j,names(ZoopTaxonomy)[rot.i]]<- NA}
  }
}

###################################################
####Biomass conversion in SizeID data frame####
###################################################

#For crustaceans, prepare the conversion equation
CrustaceanBiomassConversions$a<-exp(CrustaceanBiomassConversions$ln.a..intercept.)
#Replace all blanks with NA
CrustaceanBiomassConversions$Subclass[CrustaceanBiomassConversions$Subclass==""]<-NA
CrustaceanBiomassConversions$Order[CrustaceanBiomassConversions$Order==""]<-NA
CrustaceanBiomassConversions$Suborder[CrustaceanBiomassConversions$Suborder==""]<-NA
CrustaceanBiomassConversions$Family[CrustaceanBiomassConversions$Family==""]<-NA
CrustaceanBiomassConversions$Subfamily[CrustaceanBiomassConversions$Subfamily==""]<-NA
CrustaceanBiomassConversions$Genus[CrustaceanBiomassConversions$Genus==""]<-NA
CrustaceanBiomassConversions$Species[CrustaceanBiomassConversions$Species==""]<-NA

#Loop through each individual
for(k in 1:length(SizeID$TaxaID)){
  #If to find if the individual is a rotifer
  if(!is.na(SizeID[k,"Phylum"])&SizeID[k,"Phylum"]=="Rotifera"){
    #loop through each taxa starting with species and ending with Phylum
    for(taxa in length(names(ZoopTaxonomy)):1){
      #If there is an id there and that matches an id in the conversion dataframe, 
      #then do the conversion using the median value 
      if(!is.na(SizeID[k,names(ZoopTaxonomy)[taxa]])&SizeID[k,names(ZoopTaxonomy)[taxa]]  %in% as.matrix(RotiferBiomassConversions[names(ZoopTaxonomy)[taxa]])){
        #Find the median of all values of that taxa
        conversion.parameter.a<-median(RotiferBiomassConversions$Parameter_vol.a3[RotiferBiomassConversions[,names(ZoopTaxonomy)[taxa]]==SizeID[k,names(ZoopTaxonomy)[taxa]]],na.rm=T)
        #Convert to biovolume using the Volume (mm3) = parameter * length^3
        Volume_mm3<-conversion.parameter.a * (SizeID[k,"Size_mm"])^3
        #Convert from volume to fresh weight using a specific gravity of 1 kg/m3 (1000 ug/mm3)
        Freshweight_ug<-Volume_mm3*1000
        #Convert to dry weight assuming ratio of 0.1 (dry:wet) - dry weight = wet * 0.1
        SizeID[k,"Biomass_ug"]<-Freshweight_ug*0.1
        #break the for loop
        break
        #End of if loop to find which id matches
      }
      #End of looping through taxonomy
    }
    
    
    #End of Rotifer if else
  }else if(!is.na(SizeID[k,"Subphylum"])&SizeID[k,"Subphylum"]=="Crustacea"){
    #If it is a copepod nauplius
    if(SizeID$Nauplius[k]=="nauplius"&SizeID$Subclass[k]=="Copepoda"){  
      #Find the conversion for Copepod nauplius
      conversion.parameter.b<-CrustaceanBiomassConversions$b[is.na(CrustaceanBiomassConversions$Species)]
      #choose the a that matches that b
      conversion.parameter.a<-CrustaceanBiomassConversions$a[is.na(CrustaceanBiomassConversions$Species)]
      #Convert to dry mass using the equation w=e^(a+b(ln(L)) where L is length in mm
      SizeID[k,"Biomass_ug"]<-conversion.parameter.a*(SizeID[k,"Size_mm"]^conversion.parameter.b)
    }else{
      #Else for non-nauplii  
      #loop through each taxa starting with species and ending with Phylum
      for(taxa in length(names(ZoopTaxonomy)):1){
        #If there is an id there and that matches an id in the conversion dataframe, 
        #then do the conversion using the median value 
        if(!is.na(SizeID[k,names(ZoopTaxonomy)[taxa]]) & SizeID[k,names(ZoopTaxonomy)[taxa]]%in% as.matrix(CrustaceanBiomassConversions[names(ZoopTaxonomy)[taxa]])){
          #Remove NAs from the data frame
          temp.CrustaceanBiomassConversions<-CrustaceanBiomassConversions[!is.na(CrustaceanBiomassConversions[,names(ZoopTaxonomy)[taxa]]),]
          #Find the median of all values of that taxa for b
          conversion.parameter.b<-median.zoop(temp.CrustaceanBiomassConversions$b[temp.CrustaceanBiomassConversions[,names(ZoopTaxonomy)[taxa]]==SizeID[k,names(ZoopTaxonomy)[taxa]]])
          #choose the a that matches that b, use median.zoop in case there are several of that taxa with the same slope
          conversion.parameter.a<-median.zoop(temp.CrustaceanBiomassConversions$a[temp.CrustaceanBiomassConversions$b==conversion.parameter.b&temp.CrustaceanBiomassConversions[,names(ZoopTaxonomy)[taxa]]==SizeID[k,names(ZoopTaxonomy)[taxa]]])
          #Convert to dry mass using the equation w=e^(a+b(ln(L)) where L is length in mm
          SizeID[k,"Biomass_ug"]<-conversion.parameter.a*(SizeID[k,"Size_mm"]^conversion.parameter.b)
          
          #break the for loop
          break
          #End of if loop to find which id matches
        }
        #End of looping through taxonomy
      }
      
      #End of else for non-nauplii
    }
    
    
    #End of crustacea if else  
  }else{
    #All others get NA biomass
    SizeID$Biomass_ug[k]<-NA} 
  #End of loop through the individuals
}


####################################################
####THIS SECTION IS USED TO CONFIRM THE BIOMASS CONVERSION IS WORKING####

hist(SizeID$Biomass_ug[SizeID$Subphylum=="Crustacea"])
max(SizeID$Biomass_ug[SizeID$Phylum=="Rotifera"],na.rm=T)
head(SizeID[order(SizeID$Biomass_ug,decreasing=TRUE),c("LowestTaxonomicLevelOfID","TaxaID","Size_mm","Biomass_ug")],15)
temp<-SizeID[order(SizeID$Biomass_ug,decreasing=TRUE),]
head(temp[temp$Phylum=="Rotifera"&!is.na(temp$Phylum),c("LowestTaxonomicLevelOfID","TaxaID","Size_mm","Biomass_ug")],15)
head(temp[temp$Subphylum=="Crustacea"&!is.na(temp$Phylum),c("LowestTaxonomicLevelOfID","TaxaID","Size_mm","Biomass_ug")],15)
median(temp[temp$Subphylum=="Crustacea"&!is.na(temp$Phylum),c("Biomass_ug")],na.rm=T)
length(temp[temp$Subphylum=="Crustacea"&!is.na(temp$Phylum),c("Biomass_ug")])
length(temp[temp$Phylum=="Rotifera"&!is.na(temp$Phylum),c("Biomass_ug")])
hist(SizeID$Biomass_ug[SizeID$Phylum=="Rotifera"])
SizeID[!is.na(SizeID$Biomass_ug>100),c("LowestTaxonomicLevelOfID","TaxaID","Biomass_ug")]
names(ZoopTaxonomy)

#Check total lengths to see how many in each taxa of rotifer vs. crustacean vs. neither
length(SizeID$Biomass_ug)
sum(SizeID$Phylum=="Rotifera"&!is.na(SizeID$Phylum))
sum(SizeID$Subphylum=="Crustacea"&!is.na(SizeID$Subphylum))
sum(!is.na(SizeID$Biomass_ug))
##########################################################

#Create a vector of the sample ID
uniqueSampleID<-unique(SizeID$sample_ID)

#Initialize a data frame with the correct columns
ZoopBiomass<-data.frame(matrix(ncol = length(names(SizeID)), nrow = length(uniqueSampleID)))
names(ZoopBiomass)<-names(SizeID)
#Store the sample_no
ZoopBiomass$sample_ID<-uniqueSampleID[order(uniqueSampleID)]

#Pull the unique entries for some of the string columns
ZoopBiomass$Project<-aggregate(Project~sample_ID,data=SizeID,FUN=unique,na.action=na.omit)[,2]
ZoopBiomass$site_no<-aggregate(site_no~sample_ID,data=SizeID,FUN=unique,na.action=na.omit)[,2]
ZoopBiomass$collect_date<-as.Date(aggregate(collect_date~sample_ID,data=SizeID,FUN=unique,na.action=na.omit)[,2],"%d-%b-%y")
ZoopBiomass$Hour<-aggregate(Hour~sample_ID,data=SizeID,FUN=unique,na.action=na.omit)[,2]

ZoopBiomass$INT<-ifelse(sum(!is.na(SizeID$INT))==0,NA,aggregate(INT~sample_ID,data=SizeID,FUN=unique,na.action=na.omit)[,2])

#Keep important colums
keep<-c("Project","site_no","collect_date","Hour","sample_ID")
ZoopBiomass<-ZoopBiomass[,keep]

#Create new dataframe for biomass
ZoopFinal<-ZoopBiomass

#Overall size stats for each sample (moved biomass to after taxa calcs)
#Average size by site_ID
ZoopFinal<- merge(ZoopFinal,setNames(aggregate(SizeID$Size_mm~SizeID$sample_ID,FUN=mean,na.rm=TRUE),c("sample_ID","OverallMeanSize_mm")),by=c("sample_ID"),all.x=T)
#SD of size by site_ID
ZoopFinal<-merge(ZoopFinal,setNames(aggregate(Size_mm~sample_ID,data=SizeID,FUN=sd,na.action=na.omit),c("sample_ID","OverallSESize_mm")),by=c("sample_ID"),all.x=T)
#count of size by site_ID
ZoopFinal<-merge(ZoopFinal,setNames(aggregate(TaxaID~sample_ID,data=SizeID,FUN=length,na.action=na.pass),c("sample_ID","OverallCount_n")),by=c("sample_ID"),all.x=T)

#total biomass is calculated after taxa biomass (see below)
ZoopFinal$TotalBiomass_ug <- NA
ZoopFinal$BiomassConcentration_ugpL <- NA
###original calc was biom_unadj * count/zoops_measured, but this means that smaller zoop biomass is largely overestimated because I generally measure a lot more of the large zoops (i.e., total biomass/zoops measured is skewed towards the larger zoops). 

#Create columns broken down by the following taxa
taxaOfInterest<-c("Copepoda","Cladocera","Rotifera")
CorrespondingTaxaLevel<-c("Subclass","Suborder","Phylum")

#For loop that runs through all the taxa of interest and generates output (including nauplius!)
for(taxa.index in 1:length(taxaOfInterest)){
  taxa<-taxaOfInterest[taxa.index]
  temp<-SizeID[SizeID[,CorrespondingTaxaLevel[taxa.index]]==taxa &!is.na(SizeID[,CorrespondingTaxaLevel[taxa.index]]),]
 
  #If there are no individuals in that taxa, then export as NA
  if(nrow(temp)==0){
    ZoopFinal[,c(paste(taxa,"MeanSize_mm",sep=""),paste(taxa,"SESize_mm",sep=""))]<-NA
    ZoopFinal[,c(paste(taxa,"Count_n",sep=""),paste(taxa,"_PercentOfTotal",sep=""),paste(taxa,"_totalbiomass_ug",sep=""),paste(taxa,"_PercentOftotalbiomassConcentration",sep=""))]<-0
  }else{
    #Average size by sample_ID
    ZoopFinal<-merge(ZoopFinal,setNames(aggregate(Size_mm~sample_ID,data=temp,FUN=mean,na.action=na.omit),c("sample_ID",paste(taxa,"MeanSize_mm",sep=""))),by=c("sample_ID"),all.x=T)
    #SD of size by sample_ID
    ZoopFinal<-merge(ZoopFinal,setNames(aggregate(Size_mm~sample_ID,data=temp,FUN=sd,na.action=na.omit),c("sample_ID",paste(taxa,"SESize_mm",sep=""))),by=c("sample_ID"),all.x=T)
    #biomass calculation
    ZoopFinal<-merge(ZoopFinal,setNames(aggregate(Biomass_ug~sample_ID,data=temp,FUN=sum,na.action=na.omit),c("sample_ID",paste(taxa,"_total_biomass_ug",sep=""))),by=c("sample_ID"),all.x=T)
  }   
  #End of loop through taxa of interest  
}


#just need to sum these three taxa to calculate total biomass
TaxaforTotalCalcs<- c("Copepoda","Cladocera","Rotifera")
#Calculate raw biomass 
ZoopFinal$TotalBiomass_ug<- rowSums(ZoopFinal[,paste0(TaxaforTotalCalcs,"_total_biomass_ug")], na.rm=TRUE)
#Calculate biomass concentration as biomass relative to total volume counted 
ZoopFinal$BiomassConcentration_ugpL<-ZoopFinal$TotalBiomass_ug /ZoopFinal$Volume_L

#Replace NA with 0's for the relevant cells (count, density, biomass, biomass concentration)
ZoopFinal[c(paste(taxaOfInterest,"_total_biomass_ug",sep=""))][is.na(ZoopFinal[c(paste(taxaOfInterest,"_total_biomass_ug",sep=""))])] <- 0

#Export the final Zoop table
#write.csv(ZoopFinal,paste("FCR_ZooplanktonSummary",year,".csv",sep=""),row.names = FALSE)


ZoopFinal$BiomassConcentration_ugpL
ZoopFinal$Cladocera_BiomassConcentration_ugpL
ZoopFinal$Calanoida_BiomassConcentration_ugpL
ZoopFinal$Cladocera_PercentOftotalbiomassConcentration
ZoopFinal$Calanoida_PercentOftotalbiomassConcentration

#End of for loop going through the years
#}

