
options(stringsAsFactors = FALSE)

gr <- merge(genetic_results,all.mutations,all.x=TRUE)
gr$hgvsp <- gr$hgvsp.new
gr$hgvsc <- gr$hgvsc.new
gr$hgvsc[which(gr$hgvsc=='')]<-NA
gr$hgvsp[which(gr$hgvsp=='')]<-NA
 
MOI <- list( 'Autosomal Dominant'='AD', 'Autosomal Dominant with Reduced Penetrance'='AD-RP', 'Autosomal Recessive'='AR', 'Mitochondrial'='MT', 'X-linked Recessive'='XL')

produce_gene_summary <- function (d.mut) {
cat('num of solved families:', solved.families.count <- (nrow(unique(d.mut[which(d.mut$status=='Affected' & !is.na(d.mut$pedigree_gene)),c('pedigree_id','pedigree_gene')]))),'\n' )
cat('num of tested families:', all.families.count <- ( length(unique(d.mut[which(d.mut$event_type=='OphInGeneticresults'),'pedigree_id'])) ),'\n')
cat('solve rate:', solved.families.count / all.families.count, '\n')
	d.mut.gene_count <- as.data.frame(do.call('rbind',as.list(by(d.mut,d.mut$pedigree_gene,function(x) {
			   cat('Pedigree gene:',unique(x$pedigree_gene),'\n')
			   cat('Other genes screened:', other.genes.screened <- paste(unique(na.omit(x$gene)),collapse=', '),'\n')
                           print(nrow(x.affected <- x[which(x$status=='Affected'),]))
                           x.affected.hom <- x.affected[which(x.affected$homo==1),]
                           x.affected.nonhom <- unique(x.affected[which(x.affected$homo==0),c('hos_num','hgvsc')])
                           x.gender <- unique(x[,c('hos_num','gender')])
                           t.gender <- table(x.gender$gender)
                           gender_count<-paste(paste(names(t.gender),as.numeric(t.gender),sep=': '),collapse=', ')
                           x.affected_gender <- unique(x[which(x$status=='Affected'),c('hos_num','gender')])
                           t.affected_gender <- table(x.affected_gender$gender)
                           affected_gender_count<-paste(paste(names(t.affected_gender),as.numeric(t.affected_gender),sep=': '),collapse=', ')
                           x.ethnicity_group <- unique(x[,c('hos_num','EthnicityGroup')])
                           t.ethnicity_group <- sort(table(x.ethnicity_group$EthnicityGroup),decreasing=TRUE)
                           ethnicity_group_count <- paste(paste(names(t.ethnicity_group),as.numeric(t.ethnicity_group),sep=': '),collapse=', ')
                           print(most_frequent_ethnicity_group <- ifelse(length(names(t.ethnicity_group)),as.character(names(t.ethnicity_group)[[1]]),''))
                           x.study <- unique(x[,c('hos_num','study_name')])
                           t.study <- sort(table(x.study$study_name),decreasing=TRUE)
                           study_count <- paste(paste(names(t.study),as.numeric(t.study),sep=': '),collapse=', ')
                           x.ethnicity_name <- unique(x[,c('hos_num','EthnicityName')])
                           t.ethnicity_name <- sort(table(x.ethnicity_name$EthnicityName),decreasing=TRUE)
                           ethnicity_name_count <- paste(paste(names(t.ethnicity_name),as.numeric(t.ethnicity_name),sep=': '),collapse=', ')
                           print(most_frequent_ethnicity_name <- ifelse(length(names(t.ethnicity_name)),as.character(names(t.ethnicity_name)[[1]]),''))
                           x.disorder_term <- unique(x[,c('hos_num','disorder_term')])
                           t.disorder_term <- sort(table(x.disorder_term$disorder_term),decreasing=TRUE)
                           disorder_count <- paste(paste(names(t.disorder_term),as.numeric(t.disorder_term),sep=': '),collapse=', ')
                           print(most_frequent_disorder <- ifelse(length(names(t.disorder_term)),as.character(names(t.disorder_term)[[1]]),''))
			   x$pedigree_inheritance <- as.character(x$pedigree_inheritance)
			   x$pedigree_inheritance[which(x$pedigree_inheritance=='Unknown/other')] <- NA
			   x.inheritance  <- unique(x[,c('hos_num','pedigree_inheritance')])
                           t.inheritance <- sort(table(x.inheritance$pedigree_inheritance),decreasing=TRUE)
                           inheritance_count <- as.character(paste(paste(names(t.inheritance),as.numeric(t.inheritance),sep=': '),collapse=', '))
			   i <- which(MOI[names(t.inheritance)] %in% c('AD','AD-RP','AR','XL'))
                           print(inheritance_pct <- as.character(paste(paste(MOI[names(t.inheritance)[i]], round(100*as.numeric(t.inheritance)[i]/sum(as.numeric(t.inheritance)[i])) ,sep=': '),collapse=', ')))
                           print(most_frequent_inheritance <- ifelse(length(names(t.inheritance)),as.character(names(t.inheritance)[[1]]),''))
                           x.method <- unique(x[,c('hos_num','method')])
                           t.method <- sort(table(x.method$method),decreasing=TRUE)
                           method_count <- as.character(paste(paste(names(t.method),as.numeric(t.method),sep=': '),collapse=', '))
                           x.event_year <- unique(x[,c('hos_num','event_year')])
                           t.event_year <- table(x.event_year$event_year)
			   event_year_count <- as.character(paste(paste(names(t.event_year),as.numeric(t.event_year),sep=': '),collapse=', '))
                           x.result_year <- unique(x[,c('hos_num','result_year')])
                           t.result_year <- table(x.result_year$result_year)
			   result_year_count <- as.character(paste(paste(names(t.result_year),as.numeric(t.result_year),sep=': '),collapse=', '))
                           hospital_numbers <- paste(unique(x.affected$hos_num), collapse=';')
                           hom_patients <- unique(x.affected.hom$hos_num)
                           comp_het_patients <- setdiff(names(which(table(x.affected.nonhom$hos_num)>1)),hom_patients)
                           het_patients <- setdiff(names(which(table(x.affected.nonhom$hos_num)==1)),hom_patients)
                           min_age_hosnum <- x.affected$hos_num[which.min(x.affected$patient_age_at_first_appointment)]
                           min_age_hgvsp <- x.affected$hgvsp[which.min(x.affected$patient_age_at_first_appointment)]
                           min_age_hgvsc <- x.affected$hgvsp[which.min(x.affected$patient_age_at_first_appointment)]
                           x.hgvsp <- unique(x.affected[,c('hos_num','hgvsp')])
                           t.hgvsp <- sort(table(x.hgvsp$hgvsp),decreasing=TRUE)
                           hgvsp_counts <- paste(paste(names(t.hgvsp),as.numeric(t.hgvsp),sep=': '),collapse=', ')
                           top.hgvsp_counts <- paste(paste(head(names(t.hgvsp),5),head(as.numeric(t.hgvsp),5),sep=': '),collapse=', ')
			   x.hgvsc <- unique(x.affected[,c('hos_num','hgvsc')])
                           t.hgvsc <- sort(table(x.hgvsc$hgvsc),decreasing=TRUE)
                           hgvsc_counts <- paste(paste(names(t.hgvsc),as.numeric(t.hgvsc),sep=': '),collapse=', ')
                           top.hgvsc_counts <- paste(paste(head(names(t.hgvsc),5),head(as.numeric(t.hgvsc),5),sep=': '),collapse=', ')
                           data.frame(
			     other_genes_screened=other.genes.screened,
                             hgvsp_count=length(unique(x.affected$hgvsp)),
                             hgvsc_count=length(unique(x.affected$hgvsc)),
			     top.hgvsp_counts=top.hgvsp_counts,
			     hgvsp_counts=hgvsp_counts,
			     top.hgvsc_counts=top.hgvsc_counts,
			     hgvsc_counts=hgvsc_counts,
                             hom_count=length(hom_patients),
                             compound_het_count=length(comp_het_patients),
                             het_count=length(het_patients),
                             affected_count=length(unique(x.affected$hos_num)),
                             pedigree_count=length(unique(x.affected$pedigree_id)),
                             patient_count=length(unique(x$hos_num)),
			     inheritance_count=inheritance_count,
			     inheritance_pct=inheritance_pct,
		             most_frequent_inheritance=most_frequent_inheritance,
                             disorder_count=disorder_count,
			     most_frequent_disorder=most_frequent_disorder, 
                             ethnicity_name_count=ethnicity_name_count,
			     most_frequent_ethnicity_name=most_frequent_ethnicity_name,
                             ethnicity_group_count=ethnicity_group_count,
			     most_frequent_ethnicity_group=most_frequent_ethnicity_group,
                             gender_count=gender_count,
			     affected_gender_count=affected_gender_count,
                             method_count=method_count,
			     study_count=study_count,
			     min_years_till_result=min(x.affected$years_till_result,na.rm=TRUE),
	                     median_years_till_result=median(x.affected$years_till_result,na.rm=TRUE),
			     max_years_till_result=max(x.affected$years_till_result,na.rm=TRUE),
                             #min_age=min(x.affected$patient_age_at_first_appointment,na.rm=TRUE),
                             min_age=min(x.affected$patient_age,na.rm=TRUE),
                             min_age_hosnum=ifelse(length(min_age_hosnum)==0, '', min_age_hosnum),
                             min_age_hgvsp=ifelse(length(min_age_hgvsp)==0, '', min_age_hgvsp),
                             min_age_hgvsc=ifelse(length(min_age_hgvsc)==0, '', min_age_hgvsc),
                             median_age=median(x.affected$patient_age,na.rm=T),
                             mean_age=mean(x.affected$patient_age,na.rm=T),
                             #max_age=max(x.affected$patient_age_at_first_appointment,na.rm=T),
                             max_age=max(x.affected$patient_age,na.rm=T),
			     min_age_at_result=min(x.affected$patient_age_at_result,na.rm=T),
                             max_age_at_result=max(x.affected$patient_age_at_result,na.rm=T),
			     event_year_count=event_year_count,
                             min_event_year=min(x$event_year),
                             max_event_year=max(x$event_year),
                             result_year_count=result_year_count,
                             hospital_numbers=hospital_numbers,
                             hom_patients=paste(hom_patients,collapse=';'),
                             het_patients=paste(het_patients,collapse=';'),
                             comp_het_patients=paste(comp_het_patients,collapse=';')
                            )
                     }))))
	d.mut.gene_count$gene <- rownames(d.mut.gene_count)
	genetic_results_summary <- d.mut.gene_count
	# add transcript lengths and transcript counts
	biomart.genes <- read('~/People/OmarMahroo/biomart_genes.txt')
	gene.info <- unique(data.frame( gene=as.character(biomart.genes[,'Gene name']), gene_start=as.numeric(biomart.genes[,'Gene start (bp)']), gene_end=as.numeric(biomart.genes[,'Gene end (bp)']), transcript.length=biomart.genes[,'Transcript length (including UTRs and CDS)'], transcript.count=biomart.genes[,'Transcript count']))
	gene.info <- gene.info[which(gene.info$gene %in% genetic_results_summary$gene),]
	#pick the min and max transcript length
	gene.info <- do.call('rbind', by(gene.info,(gene.info$gene),function(x) {
		return (data.frame( gene=as.character(unique(x$gene)), gene.length=max(x$gene_end-x$gene_start), max.transcript.length=max(x$transcript.length), min.transcript.length=min(x$transcript.length), number.of.transcripts=max(x$transcript.count)))
	}))
	genetic_results_summary <- merge(genetic_results_summary,gene.info,all.x=TRUE)
	# pLI from exac
	pli <- unique(read('fordist_cleaned_nonpsych_z_pli_rec_null_data.txt')[,c('gene','pLI')])
	pli <- pli[which(pli$gene %in% genetic_results_summary$gene),]
	pli <- do.call('rbind',by(pli,pli$gene,function(x){data.frame(gene=unique(x$gene),pLI=max(x$pLI))}))
	genetic_results_summary <- merge(genetic_results_summary,pli,all.x=TRUE)
	genetic_results_summary <- genetic_results_summary[order(genetic_results_summary$pedigree_count,decreasing=TRUE),]
	rownames(genetic_results_summary) <- genetic_results_summary$gene
	genetic_results_summary$pedigree_count <- paste(genetic_results_summary$pedigree_count , ' (', round( 100* genetic_results_summary$pedigree_count / sum(genetic_results_summary$pedigree_count),2), '%)', sep='' )
	genetic_results_summary$pct_affected <- round(100*genetic_results_summary$affected_count/genetic_results_summary$patient_count,2)
	genetic_results_summary$affected_age_distribution <- paste( genetic_results_summary$min_age, genetic_results_summary$median_age, genetic_results_summary$max_age, sep='; ' )
	genetic_results_summary$affected_gender_count <- paste( genetic_results_summary$affected_count , ' (', genetic_results_summary$affected_gender_count, ')', sep='' )
	return(genetic_results_summary)
}


print('All:')
genetic_results_summary <- produce_gene_summary(gr)
write.csv(genetic_results_summary,file='gene_summary.csv',row.names=FALSE)

print('Under 18:')
genetic_results_under_18 <- gr[which(gr$patient_age<18),]
genetic_results_summary_under_18 <- produce_gene_summary(genetic_results_under_18)
write.csv(genetic_results_summary_under_18,file='gene_summary_under_18.csv',row.names=FALSE)

print('Since 2017:')
genetic_results_since_2017 <- gr[which( year(gr$AppointmentDate)>=2017),]
genetic_results_summary_since_2017 <- produce_gene_summary( genetic_results_since_2017 )
write.csv(genetic_results_summary_since_2017,file='gene_summary_since_2017.csv',row.names=FALSE)

#dim(unique(gr[which(!is.na(gr$dna_extraction_year)),c('pedigree_id','dna_extraction_year')]))
#dim(unique(gr[which(!is.na(gr$dna_extraction_year)&!is.na(gr$pedigree_gene)),c('pedigree_id','pedigree_gene')]))

#unique(gr[which(!is.na(gr$dna_extraction_year)),c('hos_num','dna_extraction_year','pedigree_gene')])

#X <- gr[which(!is.na(gr$dna_extraction_year)),]

#length(unique(gr[which(gr$event_type=='OphInGeneticresults'),'pedigree_id']))

#X$ConsultantFirmName %in% c('Moore Tony', 'Webster Andrew', 'Michaelides Michel', 'Holder Graham', 'mahroo omar')

#sum( sort( table(unique(X[which(!is.na(X$pedigree_gene)),c('hos_num','ConsultantFirmName')])$ConsultantFirmName)))




