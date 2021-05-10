
# loads the necessary tables and saves to disk
library(data.table)
library(RODBC)
library(Hmisc)
library(eeptools)


options(stringsAsFactors = FALSE)

myconn <- odbcDriverConnect('dsn=mehqv')

# list all tables!
sqlQuery(myconn, "select table_name from [mis_oe13].[information_schema].tables")

#grep('ethni', sqlQuery(myconn, "select table_name from [dw_live].[information_schema].tables"),value=T,ignore.case=TRUE)

#grep('episode', sqlQuery(myconn, "select table_name from [dw_live].[information_schema].tables"),value=T,ignore.case=TRUE)

myconn <- odbcDriverConnect('dsn=mehqv')

genes <- read('genes.csv')
rownames(genes) <- as.character(genes$gene_id)

effects <- read('effects.csv')
rownames(effects) <- as.character(effects$id)

inheritance <- data.frame(rbind(
			c("1","Autosomal Dominant"),
			c("2","Autosomal Dominant with Reduced Penetrance"),
			c("3","Autosomal Recessive"),
			c("4","Mitochondrial"),
			c("5","Simplex"),
			c("6","Unknown/other"),
			c("7","X-linked Recessive")))
colnames(inheritance) <- c('id','name')
rownames(inheritance) <- as.character(inheritance$id)

dimpat <- sqlQuery(myconn, "select * from [DW_LIVE].[patient].[dimpatient]")
write.csv(dimpat, file='dimpat.csv', row.names=FALSE)


#gene <- sqlQuery(myconn, "select * from [MIS_OE13].[ext13].[gene]")
#write.csv(gene,file='gene.csv',row.names=FALSE)


disorder <- sqlQuery(myconn, "select * from [MIS_OE13].[ext13].[disorder]")
rownames(disorder) <- as.character(disorder$id)
write.csv(disorder,file='disorder.csv',row.names=FALSE)
 
event_type <- sqlQuery(myconn, "select * from [mis_oe13].ext13.event_type")
rownames(event_type) <- as.character(event_type$id)
write.csv(event_type,file='event_type.csv',row.names=FALSE)

service <- sqlQuery(myconn,'select * from [mis_oe13].[ext13].service')
write.csv(service,file='service.csv',row.names=FALSE)


pas <- sqlQuery(myconn, "select * from [DW_LIVE].[outpatient].[PAS_appointment]")
#pas <- as.data.frame(fread('~/People/PraveenPatel/PAS.csv'))
write.csv(pas,file='PAS.csv',row.names=FALSE)

print(head(dnaextraction <- sqlQuery(myconn, "select * from [MIS_OE13].[ext13].[et_ophindnaextraction_dnaextraction]")))
rownames(dnaextraction) <- dnaextraction$id

write.csv(methods,file='ophingeneticresults_test_method.csv',row.names=FALSE)


print(head(methods <- sqlQuery(myconn, "select * from [MIS_OE13].[ext13].[ophingeneticresults_test_method]")))
rownames(methods) <- methods$id
write.csv(methods,file='ophingeneticresults_test_method.csv',row.names=FALSE)

print(head(genetics_patient <- sqlQuery(myconn, "select * from [MIS_OE13].[ext13].[genetics_patient]")))
write.csv(genetics_patient,file='genetics_patient.csv', row.names=FALSE)

print(head(genetics_patient_diagnosis <- sqlQuery(myconn, "select * from [MIS_OE13].[ext13].[genetics_patient_diagnosis]")))
write.csv(genetics_patient_diagnosis, file='genetics_patient_diagnosis.csv', row.names=FALSE)

print(head(genetics_patient_pedigree <- sqlQuery(myconn, "select * from [MIS_OE13].[ext13].[genetics_patient_pedigree]")))
write.csv(genetics_patient_pedigree, file='genetics_patient_pedigree.csv', row.names=FALSE)

print(head(genetics_patient_relationship <- sqlQuery(myconn, "select * from [MIS_OE13].[ext13].[genetics_patient_relationship]")))
write.csv(genetics_patient_relationship, file='genetics_patient_relationship.csv', row.names=FALSE)

print(head(genetics_relationship <- sqlQuery(myconn, "select * from [MIS_OE13].[ext13].[genetics_relationship]")))
rownames(genetics_relationship) <- as.character(genetics_relationship$id)
write.csv(genetics_relationship, file='genetics_relationship.csv', row.names=FALSE)

print(head(genetics_study <- sqlQuery(myconn, "select * from [MIS_OE13].[ext13].[genetics_study]")))
write.csv(genetics_study, file='genetics_study.csv', row.names=FALSE)

print(head(genetics_study_subject <- sqlQuery(myconn, "select * from [MIS_OE13].[ext13].[genetics_study_subject]")))
write.csv(genetics_study_subject,file='genetics_study_subject.csv',row.names=FALSE )

print(head(genetics_study_proposer <- sqlQuery(myconn, "select * from [MIS_OE13].[ext13].[genetics_study_proposer]")))
write.csv(genetics_study_proposer ,file='genetics_study_proposer.csv', row.names=FALSE)

print(head(et_ophindnasample_sample <- sqlQuery(myconn, "select * from mis_oe13.ext13.et_ophindnasample_sample")))
write.csv(et_ophindnasample_sample,file='et_ophindnasample_sample.csv',row.names=FALSE)

print(head(event <- sqlQuery(myconn, "select * from mis_oe13.ext13.event")))
write.csv(event,file='event.csv',row.names=FALSE)

print(head(patient <- sqlQuery(myconn, "select *  from mis_oe13.ext13.patient")))
rownames(patient) <- as.character(patient$id)
write.csv(patient,file='patient.csv',row.names=FALSE)

#sqlQuery(myconn, "select * from mis_oe13.ext13.genetics_patient")
#sqlQuery(myconn, "select * from mis_oe13.ext13.genetics_patient_pedigree")

print(head(pedigree <- sqlQuery(myconn, "select * from mis_oe13.ext13.pedigree")))
pedigree$gene_names <- genes[as.character(pedigree$gene_id),'gene_names']
pedigree$inheritance <- inheritance[as.character(pedigree$inheritance_id),'name']
colnames(pedigree) <- paste('pedigree_',colnames(pedigree),sep='')
write.csv(pedigree, file='pedigree.csv', row.names=FALSE)

print(head(ophindnasample_sample_type <- sqlQuery(myconn, "select * from MIS_OE13.ext13.ophindnasample_sample_type")))
write.csv(ophindnasample_sample_type,file='ophindnasample_sample_type.csv',row.names=FALSE)

print(head(pedigree_status <- sqlQuery(myconn, "select * from MIS_OE13.ext13.pedigree_status")))

print(head(contact <- sqlQuery(myconn,"select * from mis_oe13.ext13.contact")))

print(head(episode <- sqlQuery(myconn,"select * from mis_oe13.ext13.episode")))

print(head(service <- sqlQuery(myconn,"select * from mis_oe13.ext13.service")))

print(head(subspeciality <- sqlQuery(myconn,"select * from mis_oe13.ext13.subspecialty")))

print(head(service_subspecialty_assignment <- sqlQuery(myconn,"select * from mis_oe13.ext13.service_subspecialty_assignment")))

print(head(period <- sqlQuery(myconn,"select * from mis_oe13.ext13.period")))

print(head(event_group <- sqlQuery(myconn,"select * from mis_oe13.ext13.event_group")))

print(head(event_dna_extraction <- sqlQuery(myconn,"select * from [MIS_OE13].[ext13].[et_ophindnaextraction_dnatests]" )))

print(head(event_dna_extraction <- sqlQuery(myconn,"select * from [MIS_OE13].[ext13].[et_ophindnasample_sample_genetics_studies]" )))

# Genetic studies
#SPEED
study.speed <- na.omit(read('Studies/SPEED.csv')[,c('MEH','Barcode')])
study.speed$hos_num <- as.numeric(study.speed$MEH)
study.speed$study_id <- as.character(study.speed$Barcode)
study.speed <- do.call('rbind', by( study.speed, study.speed$hos_num, function(x) {data.frame(hos_num=unique(x$hos_num),study_name='SPEED',study_id=paste(x$study_id,collapse=';'))}))
#NIHR
study.bioresource <- na.omit(read('Studies/Bioresource_Tracker.csv')[,c('MEH no','Bardcode')])
study.bioresource$hos_num <- as.numeric(study.bioresource[,'MEH no'])
study.bioresource$study_id <- as.character(study.bioresource[,'Bardcode'])
study.bioresource <- do.call('rbind', by( study.bioresource, study.bioresource$hos_num, function(x) {data.frame(hos_num=unique(x$hos_num),study_name='NIHR',study_id=paste(x$study_id,collapse=';'))}))
#GEL
study.gel <- na.omit(read('Studies/GEL_results.csv')[,c('hosp no','GEL ID')])
study.gel$hos_num <- as.numeric(study.gel[,'hosp no'])
study.gel$study_id <- as.character(study.gel[,'GEL ID'])
study.gel <- do.call('rbind', by( study.gel, study.gel$hos_num, function(x) {data.frame(hos_num=unique(x$hos_num),study_name='GEL',study_id=paste(x$study_id,collapse=';'))}))
#UKIRDC
study.ukirdc <- na.omit(read('Studies/UKIRDC.csv')[,c('meh N','Study ID')])
study.ukirdc$hos_num <- as.numeric(study.ukirdc[,'meh N'])
study.ukirdc$study_id <- as.character(study.ukirdc[,'Study ID'])
study.ukirdc <- na.omit(study.ukirdc)
study.ukirdc <- do.call('rbind', by( study.ukirdc, study.ukirdc$hos_num, function(x) {data.frame(hos_num=unique(x$hos_num),study_name='UKIRDC',study_id=paste(x$study_id,collapse=';'))}))
studies <- rbind(study.speed,study.bioresource,study.gel,study.ukirdc)

