library(data.table)
library(RODBC)
library(Hmisc)
library(eeptools)

# loads the genetic_results main dataset

options(stringsAsFactors = FALSE)

myconn <- odbcDriverConnect('dsn=mehqv')

genetic_results <- sqlQuery(myconn, " select *
from [MIS_OE13].[ext13].[genetics_patient] gp
inner join [mis_oe13].[ext13].[patient] p on p.id=gp.patient_id  
inner join [DW_LIVE].[patient].[dimpatient] dp on p.hos_num=dp.PatKey
left join [MIS_OE13].[ext13].[genetics_patient_pedigree] gpp on gpp.patient_id=gp.id
left join [MIS_OE13].[ext13].[genetics_patient_diagnosis] gpd on gpd.patient_id=gp.id
left join [MIS_OE13].[ext13].[genetics_patient_relationship] gpr on gp.id=gpr.patient_id
left join [MIS_OE13].[ext13].[episode] ep on ep.patient_id=gp.patient_id 
left join [MIS_OE13].[ext13].[event] ev on ep.id=ev.episode_id
left join [MIS_OE13].[ext13].[et_ophingeneticresults_test] res on res.event_id=ev.id
left join [MIS_OE13].[ext13].[et_ophindnasample_sample] sample on sample.event_id=ev.id
left join [MIS_OE13].[ext13].[et_ophindnaextraction_dnaextraction] dna_extraction on dna_extraction.event_id=ev.id
left join [MIS_OE13].[ext13].[event_type] ev_type on ev_type.id=ev.event_type_id
where p.hos_num!='' and p.hos_num is not null
")
#left join [DW_LIVE].[outpatient].[PAS_appointment] pas on pas.PatKey=p.hos_num")

genetic_results <- genetic_results[ which(genetic_results$deleted==0), ]
genetic_results$effect_type <- effects[as.character(genetic_results$effect_id),'effect']
genetic_results$event_year <- format(as.Date(genetic_results$event_date),'%Y')
genetic_results$gene_names <- genes[as.character(genetic_results$gene_id),'gene_names']
genetic_results$status <- as.character(pedigree_status[as.character(genetic_results$status_id),'name'])
genetic_results$method <- methods[as.character(genetic_results$meth),'name']
genetic_results$disorder_term <- as.character(disorder[as.character(genetic_results$disorder_id),'term'])
genetic_results$relationship <- genetics_relationship[as.character(genetic_results$relationship_id),'relationship']
genetic_results$patient_id <- genetic_results$id.1
genetic_results$event_id <- genetic_results$id.6
genetic_results$event_type <- as.character(genetic_results$class_name)

# Genetic Results
genetic_results$result_comment <- genetic_results$comments.1
genetic_results$result_year <- year(genetic_results$result_date)
genetic_results$patient_age_at_result <- as.numeric(genetic_results$event_year)-as.numeric(year(genetic_results$dob))
genetic_results$dna_extraction_year <- year(genetic_results$extracted_date)
#head( genetic_results[ which( genetic_results$class_name == 'OphInDnaextraction' ) , 'extracted_date' ] )
# time  till result

# Merge with pedigree
genetic_results <- merge(genetic_results,pedigree,by.x='pedigree_id',by.y='pedigree_id',all.x=TRUE)
# fix gene names
genetic_results$pedigree_gene <- unlist(lapply(strsplit(gsub(';',',',genetic_results$pedigree_gene_names),','),`[[`,1))
genetic_results$gene <- unlist(lapply(strsplit(gsub(';',',',genetic_results$gene_names),','),`[[`,1))
genetic_results$gene <- gsub('Multiple screen \\(negative\\)','',genetic_results$gene)
# clean amino acids, HGVSp
genetic_results$amino_acid_change <- gsub('splice','',genetic_results$amino_acid_change)
genetic_results$amino_acid_change <- gsub('\\(hom\\)','',genetic_results$amino_acid_change)
genetic_results$amino_acid_change <- gsub('RNA','', gsub('N/A','',gsub(' ','',gsub('\\)','',gsub('\\(','',genetic_results$amino_acid_change)))))
genetic_results$amino_acid_change[which(is.na(genetic_results$amino_acid_change))] <- ''
genetic_results$hgvsp <- paste(genetic_results$gene, genetic_results$amino_acid_change, sep=':')
# clean base change, cDNA, HGVSc
genetic_results$base_change <- gsub('\t','',genetic_results$base_change)
genetic_results$base_change <- gsub('\\(hom\\)','',genetic_results$base_change)
genetic_results$base_change <- gsub('RNA','', gsub('N/A','',gsub(' ','',gsub('\\)','',gsub('\\(','',genetic_results$base_change)))))
genetic_results$base_change[which(is.na(genetic_results$base_change))] <- ''
genetic_results$hgvsc <- paste(genetic_results$gene_transcript, genetic_results$base_change, sep=':')
genetic_results[ which(genetic_results$hgvsc=='NA:'), 'hgvsc' ] <- NA
genetic_results[ which(genetic_results$hgvsp=='NA:'), 'hgvsp' ] <- NA
# years till result
years_till_result <- na.omit(as.data.frame(rbindlist( by(genetic_results, genetic_results$hos_num, function(x) {
return(data.frame(hos_num=unique(x$hos_num), pedigree_created_year=max(year(x$pedigree_created_date)), result_year=max(x$result_year,na.rm=TRUE), dna_extracted_year=min(x$dna_extraction_year,na.rm=TRUE))) }))))
years_till_result <- years_till_result[-which(years_till_result$dna_extracted_year==Inf),]
years_till_result$max_result <- apply( years_till_result[,c('pedigree_created_year','result_year')], 1,  max, na.rm=TRUE)
years_till_result$years_till_result <- years_till_result$max_result-years_till_result$dna_extracted_year
years_till_result[ which(years_till_result$years_till_result<0) ,'years_till_result'] <- NA


genetic_results <- merge(genetic_results,years_till_result,all.x=TRUE)

# DNA sample extraction
#genetic_results[ which(genetic_results$event_type == 'OphInDnaextraction'),]

# Genetic studies
genetic_results <- merge( genetic_results, studies, all.x=TRUE)

# merge with PAS for appointments
genetic_results$PatKey <- as.numeric(genetic_results$hos_num)
print(length(unique(genetic_results$PatKey)))
pas$PatKey <- as.numeric(pas$PatKey)
genetic_results <- merge( genetic_results, pas, by=c('PatKey'), all.x=TRUE )
print(length(unique(genetic_results$PatKey)))

# patient age
genetic_results$patient_age <- ifelse(is.na(year(genetic_results$date_of_death)), year(Sys.Date())-year(genetic_results$dob), year(genetic_results$date_of_death)-year(genetic_results$dob))
genetic_results$patient_age_at_appointment <- year(genetic_results$AppointmentDate)-year(genetic_results$dob)

# try to calculate age at first appointment
genetic_results$event_year<-as.numeric(genetic_results$event_year)
patient_age_at_first_appointment <- as.data.frame(rbindlist( by(genetic_results, genetic_results$hos_num, function(x) {
return(data.frame(
hos_num=unique(x$hos_num),
patient_age_at_first_appointment=min(min(year(x$AppointmentDate),na.rm=TRUE)-min(year(x$dob)), min(min(x$event_year,na.rm=T)-min(year(x$dob))),min(x$patient_age),na.rm=TRUE)))})))
genetic_results <- merge(genetic_results,patient_age_at_first_appointment,all.x=TRUE)


# remove columns with all NAs
genetic_results <- genetic_results[, -which(colSums(is.na(genetic_results))==nrow(genetic_results)) ]
# remove columns with all same id
#genetic_results <- genetic_results[, -which( apply(genetic_results,2,function(x){length(unique(x))}) == 1)]

# upon's Michel's recommendation
# Leber's Congential Amaurosis / Severe Early-Onset Severe Retinal Dystrophy
genetic_results$disorder_term <- gsub("Leber's amaurosis","LCA / SEORD", genetic_results$disorder_term)
# 
genetic_results$disorder_term <- gsub( "Rod monochromatism","Achromatopsia", genetic_results$disorder_term )
#
genetic_results$disorder_term <- gsub( "Adult vitelliform macular dystrophy","Vitelliform macular dystrophy", genetic_results$disorder_term )


write.csv(genetic_results,file='genetics_results.csv',row.names=FALSE)


#d <- read('genetics_results.csv')

# exon numbers
#d2 <- read('hgvs_variants.csv')

#d3 <- merge(d,d2,all.x=TRUE,by='hgvs')

#write.csv(d3,file='genetics_results2.csv',quote=TRUE,row.names=FALSE)

# time till result
#unique(genetic_results[,c('hos_num','patient_age_at_result','gene','pedigree_gene')])
#table(genetic_results$method)
#table(genetic_results$assay)


