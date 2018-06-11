# Script to pull Hemlock sites and multi-taxa site from the Neotoma Paleoecology Database
# And eliminate sites that don't meet minimum resolution and chronology requirements
# Copyright 2017 by M. Allison Stegner

########################################################################################
###############################    FUNCTIONS     #######################################
########################################################################################

# has.chron______________________________________________________
# function to select ids for datasets that have a chonolorgy in Neotoma db
# tax_dec_dl is a Neotoma download object

has.chron<-function(tax_dec_dl){
	#chron.table<-c()
	#chron.type<-c()
	ind<-c()
	for (i in 1:length(tax_dec_dl)) {
		sitei<-tax_dec_dl[[i]]
		sitei.chron<-tax_dec_dl[[i]]$chronologies #all the chronologies for site i
		i.chrons<-names(sitei.chron)
		n.pollen.samples<-nrow(tax_dec_dl[[i]]$counts) #number of samples for site 
		if (i.chrons=="No chronology"){ #toss sites with no chronologies 
			ind[i]<-NA 
		} else {
			ind[i]<-sitei$dataset$dataset.meta$dataset.id
		}
	}
return(has.chron=ind)
}


#ts.min.length_____________________________________________________
# function to select ids for datasets that have number of data points >= minimum samples
# tax_dec_dl is a Neotoma download object
# min.samples is an integer: sites with number of samples less than min.samples are excluded

ts.min.length<-function(tax_dec_dl,min.samples){
	ind<-c()
	for (i in 1:length(tax_dec_dl)) {
		sitei<-tax_dec_dl[[i]]
		n.pollen.samples<-nrow(sitei$counts) #number of samples for site i
		if (n.pollen.samples<min.samples){ 
			ind[i]<-NA
		} else {
			ind[i]<-sitei$dataset$dataset.meta$dataset.id
		}
	}
		return(adequate.n=ind)
}
	
	
#select.high.res________________________________________________________
# function to select ids for datasets where number of years represented per pollen sample is less than max.grain
# tax_dec_dl is a Neotoma download object
# max.grain is the maximum allowable number of years represented per pollen sample

select.high.res<-function(tax_dec_dl,max.grain){ 
	chron.table<-c()
	chron.type<-c()
	for (i in 1:length(tax_dec_dl)) {
		sitei.chron<-tax_dec_dl[[i]]$chronologies #all the chronologies for site i
		i.chrons<-names(sitei.chron)
		n.pollen.samples<-nrow(tax_dec_dl[[i]]$counts) #number of samples for site i
			
		#develop a table summarizing various aspects including average resolution
		chron.data<-c()
		for (j in 1:length(i.chrons)){
			chronj<-sitei.chron[[i.chrons[j]]]
			if (nrow(chronj)<2){	next }
			i.names<-tax_dec_dl[[i]]$dataset$site.data$site.name
			n.dates<-nrow(chronj)
			chronj.dur<-max(chronj$age)-min(chronj$age)
			dates.per.time<-chronj.dur/n.dates
			yrs.per.pollen<-chronj.dur/n.pollen.samples
			chron.data1<-c(i.names,unique(chronj$dataset.id),names(sitei.chron)[j],n.dates,unique(chronj$age.type)[1],dates.per.time,yrs.per.pollen)
			chron.data<-rbind(chron.data,chron.data1)
			}
		chron.table<-rbind(chron.table,chron.data)
	}
	rownames(chron.table)<-c(1:nrow(chron.table))
	chron.table<-as.data.frame(chron.table)
	colnames(chron.table)<-c("dataset.name","dataset.id","chronology.type","n.dates","age.type","n.year.per.date","yrs.per.pollen")

	#choose sites where duration/n.samples is less than max.grain yrs
	resolution.include.exclude<-ifelse(as.numeric(as.vector(chron.table$yrs.per.pollen))>=max.grain,"exclude","include") 
	dated_tax<-cbind(chron.table,resolution.include.exclude)
	dated_tax<-dated_tax[dated_tax$resolution.include.exclude %in% "include",]
	dataset.list<-unique(as.character(dated_tax$dataset.id)) #list of dataset ids
	
	high.res.sites<-unique(as.numeric(as.character(dated_tax$dataset.id)))
	return(high.res.sites)
}


#min.chron.control________________________________________________
# function to select ids for datasets where number of chron controls is adequate
# tax_dec_dl is a Neotoma download object
# min.chron.grain is the minimum allowable number of chronolgy controls for the length of the record

min.chron.control<-function(tax_dec_dl,max.chron.grain){
	well.dated.sites<-c()
	tsuga.chrons<-c()
	for (i in 1:length(tax_dec_dl)){
		print(i)
		sitei<-tax_dec_dl[[i]]
		
		sitei.id<-sitei$dataset$dataset.meta$dataset.id
		controls<-get_chroncontrol(sitei) #this is also slow, but unavoidable (?)
		print(controls)
		
		#get duration of the time series
		duration<-max(sitei$sample.meta$age,na.rm=T)-min(sitei$sample.meta$age,na.rm=T)
		
		bad.chrons<-c("Tsuga decline","Biostratigraphic, pollen","Sediment stratigraphic","Guess","Interpolated")
		controls.temp<-cbind(as.character(controls[[1]]$control.type),controls[[1]]$age)
		
		if (length(which(as.numeric(controls.temp[,2])>200))==1) {
			controls.temp2<-controls.temp[which(as.numeric((controls.temp[,2]))>200)]
			xx<-(controls.temp2 %in% bad.chrons)
		} else if (length(which(as.numeric(controls.temp[,2])>200))==0){
			xx<-FALSE
		} else {
			controls.temp2<-controls.temp[which(as.numeric((controls.temp[,2]))>200),]
			xx<-(controls.temp2[,1] %in% bad.chrons)
		}
		
		if (length(which(xx==TRUE))==0) {
			nchrons<-length(controls[[1]]$control.type)
		} else {
			nchrons<-length(controls[[1]]$control.type)-length(which(xx==TRUE))
		} 
		
		if (nchrons==1) { #if there is only 1 chron control, toss
			next
		} else if (nchrons==0) {
			next
		} else if (duration/nchrons>max.chron.grain) { #toss sites with a max grain of more than 1 chron control per XXX years (2000 is reasonable)
			next
		} else { #if N chron controls is adequate, send ids into a vector
			well.dated.sites<-c(well.dated.sites,sitei$dataset$dataset.meta$dataset.id)			
		}
	}
	return(well.dated.sites)
}
	
	
	
#min_pol_pct________________________________________________________
# function to select ids for datasets where a single taxon pollen reaches a minimum %
# tax_dec_dl is a Neotoma download object
# eco.group is a vector of Neotoma database ecological group codes.
# taxon is a species or group name. * indicates partial matching
# min.pct is the pollen % cut off; if the site never acheives min.pct of taxon pollen, it is excluded
# pct.zeros is the maximum allowable number of pollen samples where taxon is unsampled
# if eco.sort="pulished" function chooses a published taxon list to use for calculating pollen %, else eco.group is used

min_pol_pct<-function(tax_dec_dl,eco.group,taxon,min.pct,pct.zeros,eco.sort){
	dataset.list<-c()
	for (i in 1:length(tax_dec_dl)){
		print(i)
		sitei<-tax_dec_dl[[i]]
		dataset.list[i]<-sitei$dataset$dataset.meta$dataset.id
	}
	high.pct<-c()
	max.pct<-c()
	Hpct.ids<-c()
	Hpct.names<-c()
	for (k in 1:length(dataset.list)){
		sitei<-tax_dec_dl[[k]]
		sitei.counts<-sitei$counts
		all_taxa <- do.call(rbind.data.frame, lapply(tax_dec_dl, function(x)x$taxon.list[,1:6]))
		all_taxa <- all_taxa[!duplicated(all_taxa),]
	if (eco.sort=="published"){
		counts.subset<-compile_taxa(sitei,list.name="WS64")
		sitei.counts<-counts.subset$counts[,-grep("Lycopodium", colnames(sitei.counts))]
	} else if (eco.sort=="ecological.group") {
		good_cols<-c(which(colnames(sitei.counts) %in% all_taxa[all_taxa$ecological.group %in% eco.group,1]))	
		sitei.counts<-sitei.counts[,good_cols]
	} else {
		print("choose eco.sort = 'published' or 'ecological group'")
	}

		tax_pct<-sitei.counts[,1:ncol(sitei.counts)]/rowSums(sitei.counts[,1:ncol(sitei.counts)],na.rm=TRUE)
		sitei.taxon<-tax_pct[,grep(taxon, colnames(tax_pct))]
		max.pct[k]<-max(sitei.taxon[])
	
		if (max(sitei.taxon,na.rm=TRUE)<min.pct){
			next
		} else if ((sum(sitei.taxon==0)/length(sitei.taxon))>pct.zeros) {
			next
		} else {
			high.pct<-c(high.pct,k)
		}
	}
	Hpct.ids<-dataset.list[high.pct]
}


#trim_to________________________________________________________
# function to trim dataset to a minimum age and determine if there are enough remaining datapoints
# tax_dec_dl is a Neotoma download object
# min.age is the age cut off: datapoints younger than min.age will by trimmed
# cutoff is the minimum allowable number of remaining points 

trim_to<-function(tax_dec_dl,min.age,cutoff){
	ind<-c()
	for (i in 1:length(tax_dec_dl)) {
		sitei<-tax_dec_dl[[i]]
		age<-sitei$chronologies[[1]]$age
		length(age[which(age>min.age)])
		if (length(age[which(age>min.age)])<cutoff){ 
			ind[i]<-NA
		} else {
			ind[i]<-sitei$dataset$dataset.meta$dataset.id
		}
	}
		return(adequate.n=ind)
}

# pull.multi.spp________________________________________________________
pull.multi.spp<-function(hem_dec_dl,j,TAXON){
	sitei<-hem_dec_dl[[j]]$counts
	sort(colnames(sitei))
	if (length(grep("Lycopodium", colnames(sitei)))==0){
		sitei<-sitei
	} else {
		sitei<-sitei[,-grep("Lycopodium", colnames(sitei))]
	}
	
	all_taxa <- do.call(rbind.data.frame, lapply(hem_dec_dl, function(x)x$taxon.list[,1:6]))
	all_taxa <- all_taxa[!duplicated(all_taxa),]
	good_cols<-c(which(colnames(sitei) %in% all_taxa[all_taxa$ecological.group %in% c("TRSH","UPHE"),1]))	
	sitei<-sitei[,good_cols]
	#sitei<-compile_taxa(sitei,list.name="WS64")
	sitei.pct<-sitei[,1:ncol(sitei)]/rowSums(sitei[,1:ncol(sitei)],na.rm=TRUE)

	#taxa.list<-TAXON
	taxoni<-sitei.pct[,grep(TAXON,colnames(sitei.pct),ignore.case=TRUE)]
	if (length(taxoni)==0) taxoni<-rep(NA,nrow(sitei.pct)) 
	if (is.matrix(taxoni)) taxoni<-rowSums(taxoni)
	
	taxon.matx<-taxoni
	
	AM<-hem_dec_dl[[j]]$sample.meta$age
	
	taxon.matx2<-cbind(AM,taxon.matx)
	colnames(taxon.matx2)<-c("age",TAXON)

	write.csv(taxon.matx2,paste(TAXON,"_",names(hem_dec_dl[j]),".csv",sep=""),row.names=F)
}

#map_dl________________________________________________________
# function to map location of sites
# tax_dec_dl is a Neotoma download object
# X is a vector with 2 elements: min and max longitude
# Y is a vector with 2 elements: min and max latitude
# add: should points be added to an existing map?
# color: color to use for points
# label.sites: should site names be added as text to the map

map_dl<-function(tax_dec_dl,X,Y,add,color,label.sites){
	if (add==FALSE){
		map("world",xlim=X,ylim=Y)
	} 
	lat<-c()
	long<-c()
	dataset.id<-c()
	site.name<-c()
	site.id<-c()
	for (i in 1:length(tax_dec_dl)){
		long[i]<-tax_dec_dl[[i]]$dataset$site.data$long
		lat[i]<-tax_dec_dl[[i]]$dataset$site.data$lat
		points(long[i],lat[i],pch=16,col=color)
		dataset.id[i]<-tax_dec_dl[[i]]$dataset$dataset.meta$dataset.id
		site.name[i]<-tax_dec_dl[[i]]$dataset$site.data$site.name
		site.id[i]<-tax_dec_dl[[i]]$dataset$site.data$site.id
	}
	if (label.sites==TRUE){
		text(long,lat,site.id,cex=0.5,pos=4,offset=0.2)
	}
	
	out<-cbind(site.name,dataset.id,site.id,lat,long)
	return(out)
}

# a simpler version on map_dl that does not produce a map
map_dl_2<-function(tax_dec_dl){
	# if (add==FALSE){
	# 	map("world",xlim=X,ylim=Y)
	# } 
	lat<-c()
	long<-c()
	dataset.id<-c()
	site.name<-c()
	site.id<-c()
	for (i in 1:length(tax_dec_dl)){
		long[i]<-tax_dec_dl[[i]]$dataset$site.data$long
		lat[i]<-tax_dec_dl[[i]]$dataset$site.data$lat
		# points(long[i],lat[i],pch=16,col=color)
		dataset.id[i]<-tax_dec_dl[[i]]$dataset$dataset.meta$dataset.id
		site.name[i]<-tax_dec_dl[[i]]$dataset$site.data$site.name
		site.id[i]<-tax_dec_dl[[i]]$dataset$site.data$site.id
	}
	# if (label.sites==TRUE){
	# 	text(long,lat,site.id,cex=0.5,pos=4,offset=0.2)
	# }
	
	out<-cbind(site.name,dataset.id,site.id,lat,long)
	return(out)
}

#END FUNCTIONS________________________________________


library(neotoma)
#library(dplyr)
library(maps)

########################################################################################
############################  multi-taxa analyses  #####################################
########################################################################################

source("~/map_dl_7May2018.R")
source("~/Lottery_Model_functions_7May2018.R")

ALLTAXON <- c("Alnus*", "Fagus*", "Juglans*", "Picea*", "Platanus*",  "Tilia*", "Tsuga*", "Ulmus*", "Pinus banksiana", "Pinus banksiana-type", "Pinus resinosa-type", "Pinus subg. Pinus", "Pinus strobus", "Pinus subg. Strobus", "Pinus subg. Strobus undiff.", "Pinus det. P. strobus")

alltaxon <- c("Alnus", "Fagus", "Juglans", "Picea", "Platanus",  "Tilia", "Tsuga", "Ulmus", "Pinus banksiana", "Pinus banksiana-type", "Pinus resinosa-type", "Pinus subg. Pinus", "Pinus strobus", "Pinus subg. Strobus", "Pinus subg. Strobus undiff.", "Pinus det. P. strobus")


setwd("~/final_data/") ##choose where to save the new files
# contains all names from tn (above)  that contains "Pinus" 

LL <- length(ALLTAXON)
for(i in 1:LL){
	print(i)
	# i <- 1
	TAXON <- as.character(ALLTAXON[i])
	print(TAXON)
	tax_dec<-get_dataset(taxonname = TAXON, datasettype = "pollen", loc = c(-105, 25, -60, 55))
	if(length(tax_dec) == 0) next

	# download datasets from Neotoma database	
	# this next line may take up to 30 min to run	
	tax_dec_dl<-get_download(tax_dec, verbose = FALSE)

	#limit to sites with chronologies
	pol_chron<-has.chron(tax_dec_dl) #expect warnigns here. Not an issue
	pol_chron<-pol_chron[complete.cases(pol_chron)]
	pol_dl_sub1<-tax_dec_dl[as.character(pol_chron)]
	if(length(pol_dl_sub1) == 0) next

	#limit to sites with at least n samples
	pol_n<-ts.min.length(pol_dl_sub1,20)
	pol_n<-pol_n[complete.cases(pol_n)]
	pol_dl_sub2<-pol_dl_sub1[as.character(pol_n)]
	if(length(pol_dl_sub2) == 0) next

	#limit to sites with minimum level of average resolution
	pol_high<-select.high.res(pol_dl_sub2,200)
	pol_dl_sub3<-pol_dl_sub2[as.character(pol_high)]
	if(length(pol_dl_sub3) == 0) next

	#limit to sites with adequate chron controls
	pol_chroncont<-min.chron.control(pol_dl_sub3,2000)
	pol_dl_sub4<-pol_dl_sub3[as.character(pol_chroncont)]
	if(length(pol_dl_sub4) == 0) next

	# limit to sites that reach minimum % pollen
	# and for which Tsuga pollen is sampled in at least 50% of samples
	eco<-c("TRSH", "UPHE")
	min.pct<-0.0
	pct.zeros<-0.5

	taxon <- alltaxon[i]
	min.tax<-min_pol_pct(pol_dl_sub4,eco,taxon,min.pct,pct.zeros,"ecological.group")	
	pol_dl_sub5<-pol_dl_sub4[as.character(min.tax)]
	if(length(pol_dl_sub5) == 0) next


	# limit to sites with at least 20 datapoints older than 1000 years BP
	ids<-trim_to(pol_dl_sub5,1000,20)
	id.list<-ids[complete.cases(ids)]
	pol_dl_sub6<-pol_dl_sub5[as.character(id.list)]
	length(pol_dl_sub6)  
	if(length(pol_dl_sub6) == 0) next

	# generate csv files
	tax_dec_dl<-pol_dl_sub6

	for (j in 1:length(tax_dec_dl)){
		print(j)
		pull.multi.spp(tax_dec_dl,j, taxon)
	}
	nameout <- paste(taxon, "_site_info.csv", sep = "")
	out <- map_dl_2(tax_dec_dl)
	write.csv(out, nameout)	
}

###########################      FUNCTIONS ends     ####################################
########################################################################################

########################################################################################
########################### TSUGA specific analyses ####################################
########################################################################################

# query Neotoma database
# requires internet connection
hem_dec<-get_dataset(taxonname = "Tsuga*", 
	datasettype = "pollen",
	loc = c(-105, 25, -60, 55))

# download datasets from Neotoma database	
# this next line may take up to 30 min to run	
hem_dec_dl<-get_download(hem_dec)

#limit to sites with chronologies
pol_chron<-has.chron(hem_dec_dl) #expect warnigns here. Not an issue
pol_chron<-pol_chron[complete.cases(pol_chron)]
pol_dl_sub1<-hem_dec_dl[as.character(pol_chron)]

#limit to sites with at least n samples
pol_n<-ts.min.length(pol_dl_sub1,20)
pol_n<-pol_n[complete.cases(pol_n)]
pol_dl_sub2<-pol_dl_sub1[as.character(pol_n)]

#limit to sites with minimum level of average resolution
pol_high<-select.high.res(pol_dl_sub2,200)
pol_dl_sub3<-pol_dl_sub2[as.character(pol_high)]
length(pol_dl_sub3)

#limit to sites with adequate chron controls
pol_chroncont<-min.chron.control(pol_dl_sub3,2000)
pol_dl_sub4<-pol_dl_sub3[as.character(pol_chroncont)]
length(pol_dl_sub4)

# limit to sites that reach minimum % Tsuga pollen
# and for which Tusga pollen is sampled in at least 50% of samples
eco<-c("TRSH", "UPHE")
taxon<-"Tsuga*"
min.pct<-0.1
pct.zeros<-0.5

min.tsuga<-min_pol_pct(pol_dl_sub4,eco,"Tsuga*",0.1,0.5,"ecological.group")	
pol_dl_sub5<-pol_dl_sub4[as.character(min.tsuga)]
length(pol_dl_sub5)

# limit to sites with at least 20 datapoints older than 1000 years BP
ids<-trim_to(pol_dl_sub5,1000,20)
id.list<-ids[complete.cases(ids)]
pol_dl_sub6<-pol_dl_sub5[as.character(id.list)]

# map sites
X<-c(-180,-50)
Y<-c(10,90)
mapX<-map_dl(pol_dl_sub6,X,Y,add=FALSE,color="blue",label.sites=FALSE)


###################################################################################################################
############################                        BCP ANALYSES            #######################################
###################################################################################################################


### installl bcp package 
# install.packages(bcp) 

# loads the package to run 'bcp'
library('bcp')


########## MAKE SURE TO CHANGE THE DIRECTORY ###############

#############################################################
################ HEMLOCK DATA FINAL BCP FIT #################
#############################################################

#------------------------------------------------------------------------#
## This code is used for figures that showed hemlock empirical data set ##
#------------------------------------------------------------------------#
setwd("/Users/gola/Box Sync/projects/ACES_paleo_lottery/data/hemlock_data")

# sites ID extracted from neotoma
ID <- c(238, 494, 503, 516, 532, 795, 871, 971, 973, 986, 1105, 1136, 1137,1615, 1698, 1820, 1981, 2023, 2292, 2332, 3058, 3131, 3454, 13047, 14410, 15032,15209, 15350, 15660, 15682, 15732, 15764, 15866, 15892, 15904, 15909, 15931, 16189, 16231, 16271, 17387, 17404, 18110, 19842,20397, 20512, 20627)

for (id in ID){
	# import data
	name <- paste('tsuga_',id,'.csv', sep = "")
	dat <- read.csv(name, header = T)
	
	r <- bcp( y = as.numeric(dat[,2]), w0 = 0.2, p0 = 0., burnin = 1000, mcmc = 10000) # runs the bcp analysis

	nr <- dim(dat)[1]

	res1 <- cbind(dat[,1],r$posterior.prob)
	res2 <- cbind(dat[,1],r$posterior.mean)

	# export the data where _pp_ and _pm_ differentiate posterior probability and posterior mean
	name1 <- paste('tsuga_pp_',id,'.csv', sep = "")
	name2 <- paste('tsuga_pm_',id,'.csv', sep = "")
	write.csv(res1,name1, row.names = FALSE)
	write.csv(res2,name2, row.names = FALSE)
}


#----------------------------------------------------------------------------#
## This code is used for multi-taxa analyses ##
#----------------------------------------------------------------------------#
setwd("/Users/gola/Box Sync/projects/ACES_paleo_lottery/data/clean_allsp_data/")
library("bcp")

idSP <- read.csv("/Users/gola/Box Sync/projects/ACES_paleo_lottery/data/all_clean_ID.csv", header = F)

for(i in 1:nrow(idSP)){
	sp <- idSP[i,1]
	id <- idSP[i,2]
	name <- paste(sp,'_',id,'.csv', sep = "")
	print(name)
	dat <- read.csv(name, header = F)
	res1 <- dat
	res2 <- dat
	r <- bcp( y = as.numeric(dat[,2]) , w0 = 0.2, p0 = 0., burnin = 1000, mcmc = 10000) # runs the bcp analysis
	res1[,2] <- r$posterior.prob
	res2[,2] <- r$posterior.mean
	
	# export the data where _pp_ and _pm_ differentiate posterior probability and posterior mean
	name1 <- paste(sp,'_pp_',id,'.csv', sep = "")
	name2 <- paste(sp,'_pm_',id,'.csv', sep = "")
	write.table(res1,name1, sep = ",", row.names = FALSE, col.names = FALSE)
	write.table(res2,name2, sep = ",", row.names = FALSE, col.names = FALSE)
}


######################################################
################ SIMULATED  EXAMPLE ##################
######################################################

#--------------------------------------------------------------------#
##                  This code is used for figure 1B                 ##
#--------------------------------------------------------------------#


# use bcp to fit simulated data that resembles one hemlock time-series
dat <- read.csv("/Users/gola/Box Sync/projects/ACES_paleo_lottery/data/simulated_data/fig1/sim_s_0.2_seed_3.csv", header = F)

r <- bcp(dat$V1, w0 = .2, p0 = 0., burnin = 1000, mcmc = 10000)

res1 <- r$posterior.prob
res2 <- r$posterior.mean

name1 <- paste('/Users/gola/Box Sync/projects/ACES_paleo_lottery/data/simulated_data/fig1/sim_pp_s_0.2_seed_3.csv', sep = "")
name2 <- paste('/Users/gola/Box Sync/projects/ACES_aleo_lottery/data/simulated_data/fig1/sim_pm_s_0.2_seed_3.csv', sep = "")
write.csv(res1,name1,  row.names = FALSE)
write.csv(res2,name2,  row.names = FALSE)

# for the Fagus site used in fig 1C, use the fitted data from the mult-taxa analyes

######################################################
############## FIT SIMULATED DATA ####################
#####    WITH DIFFERENT DENSITY DEPENDENCE      ######
######################################################

#--------------------------------------------------------------------#
##                  This code is used for figure 2E                 ##
#--------------------------------------------------------------------#

setwd('/Users/gola/Box Sync/projects/ACES_paleo_lottery/data/simulated_data/fig2')

SS <- seq(0,10,1) # SS is the vector ID of different strength of the density-dependences (s).

for(ss in SS){
	# sed is the ID of the seed used in the simulation
	for (sed in 1:100) {
		name <- paste('sim_s_',ss,'_seed_',sed,'.csv', sep = "")
		dat <- read.csv(name, header = F)
		
		# nr is the number of sites exported from the simulated data (should be 32)
		# nc is the number of points in the time series
		nr <- dim(dat)[1]
		nc <- dim(dat)[2]
		
		# initialize matrix that contains the fitted bcp for each posterior probability and posterior mean
		res1 <- matrix(0., nrow = nr, ncol =  nc)
		res2 <- matrix(0., nrow = nr, ncol =  nc)
		# fits each time-series for each site
		for(i in 1:nr){
			r <- bcp(as.numeric(dat[i,]), w0 = .2, p0 = 0., burnin = 1000, mcmc = 10000)
			res1[i,] <- r$posterior.prob
			res2[i,] <- r$posterior.mean

		}
		name1 <- paste('sim_pp_s_',ss,'_seed_',sed,'.csv', sep = "")
		name2 <- paste('sim_pm_s_',ss,'_seed_',sed,'.csv', sep = "")
		write.csv(res1,name1,  row.names = FALSE)
		write.csv(res2,name2,  row.names = FALSE)
	}
}


#--------------------------------------------------------------------#
##                  This code is used for figure 2F                 ##
#--------------------------------------------------------------------#
setwd('/Users/gola/Box Sync/projects/ACES_paleo_lottery/data/simulated_data/fig2_same_I')

SS <- seq(0,10,1) # SS is the vector ID of different strength of the density-dependences (s).

for(ss in SS){
	# sed is the ID of the seed used in the simulation
	for (sed in 1:100) {
		name <- paste('sim_s_',ss,'_seed_',sed,'.csv', sep = "")
		dat <- read.csv(name, header = F)
		
		# nr is the number of sites exported from the simulated data (should be 32)
		# nc is the number of points in the time series
		nr <- dim(dat)[1]
		nc <- dim(dat)[2]
		
		# initialize matrix that contains the fitted bcp for each posterior probability and posterior mean
		res1 <- matrix(0., nrow = nr, ncol =  nc)
		res2 <- matrix(0., nrow = nr, ncol =  nc)
		# fits each time-series for each site
		for(i in 1:nr){
			r <- bcp(as.numeric(dat[i,]), w0 = .2, p0 = 0., burnin = 1000, mcmc = 10000)
			res1[i,] <- r$posterior.prob
			res2[i,] <- r$posterior.mean

		}
		name1 <- paste('sim_pp_s_',ss,'_seed_',sed,'.csv', sep = "")
		name2 <- paste('sim_pm_s_',ss,'_seed_',sed,'.csv', sep = "")
		write.csv(res1,name1,  row.names = FALSE)
		write.csv(res2,name2,  row.names = FALSE)
	}
}

######################################################
############## FIT SIMULATED DATA ####################
########  WITH SPATIAL AUTOCORRELATION    ############
######################################################

#--------------------------------------------------------------------#
##      This code is used for figures  4BD, 5BDF and 6A             ##
#--------------------------------------------------------------------#

setwd('/Users/gola/Box Sync/projects/ACES_paleo_lottery/data/simulated_data/fig5a')

RR <- seq(0,8, by = 1) # RR is the vector ID for different parameter values for spatial autocorrelation (rho).
lr <- length(RR)

for(k in 1:lr){
	# sed is the ID of the seed used in the simulation
	for(sed in 1:100){
	rr <- RR[k]
	name <- paste('sim_r_',rr,'_seed_',sed,'.csv', sep = "")
	dat <- read.csv(name, header = F)
	
	# nr is the number of sites exported from the simulated data (should be 32)
	# nc is the number of points in the time series
	nr <- dim(dat)[1]
	nc <- dim(dat)[2]
		
	# initialize matrix that contains the fitted bcp for each posterior probability and posterior mean			
	res1 <- matrix(0., nrow = nr, ncol =  nc)
	res2 <- matrix(0., nrow = nr, ncol =  nc)

	# fits each time-series for each site
	for(i in 1:nr){
		r <- bcp(as.numeric(dat[i,]), w0 = .2, p0 = 0., burnin = 1000, mcmc = 10000)
		res1[i,] <- r$posterior.prob
		res2[i,] <- r$posterior.mean
	}
	name1 <- paste('sim_pp_r_',rr,'_seed_',sed,'.csv', sep = "")
	name2 <- paste('sim_pm_r_',rr,'_seed_',sed,'.csv', sep = "")
	write.csv(res1,name1, row.names = FALSE)
	write.csv(res2,name2, row.names = FALSE)
	}
}


######################################################
############## FIT SIMULATED DATA ####################
########    WITH DIFFERENT DISPERSAL      ############
######################################################

#--------------------------------------------------------------------#
##                 This code is used for figure 6B                  ##
#--------------------------------------------------------------------#

setwd('/Users/gola/Box Sync/projects/ACES_paleo_lottery/data/simulated_data/fig5b')

FF <- seq(0,10, by = 1) # FF is the vector ID for different parameter values for dispersal (f).
lf <- length(FF)

for(k in 1:lf){
	# sed is the ID of the seed used in the simulation
	for(sed in seq(1,100, by = 1)){
		ff <- FF[k]
		name <- paste('sim_f_',ff,'_seed_',sed,'.csv', sep = "")
		dat <- read.csv(name, header = F)
		
		# nr is the number of sites exported from the simulated data (should be 32)
		# nc is the number of points in the time series
		nr <- dim(dat)[1]
		nc <- dim(dat)[2]
				
		# initialize matrix that contains the fitted bcp for each posterior probability and posterior mean		
		res1 <- matrix(0., nrow = nr, ncol =  nc)
		res2 <- matrix(0., nrow = nr, ncol =  nc)
		# fits each time-series for each site
		for(i in 1:nr){
			r <- bcp(as.numeric(dat[i,]), w0 = .2, p0 = 0., burnin = 1000, mcmc = 10000)
			res1[i,] <- r$posterior.prob
			res2[i,] <- r$posterior.mean
		}
		name1 <- paste('sim_pp_f_',ff,'_seed_',sed,'.csv', sep = "")
		name2 <- paste('sim_pm_f_',ff,'_seed_',sed,'.csv', sep = "")
		write.csv(res1,name1, row.names = FALSE)
		write.csv(res2,name2, row.names = FALSE)
	}
}


