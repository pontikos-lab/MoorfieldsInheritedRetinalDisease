
## numbers and figures for the paper

#genetic_results <- read('genetics_results.csv')

sprintf("A total of %s patients (%s families) out of %s (%s families) had genetic screening data available.\n",
length(unique(genetic_results[which(!is.na(genetic_results$method)&genetic_results$status=='Affected'),]$PatKey)),
length(unique(genetic_results[which(!is.na(genetic_results$method)&genetic_results$status=='Affected'),]$pedigree_id)),
length(unique(genetic_results[which(genetic_results$status=='Affected'),]$PatKey)),
length(unique(genetic_results[which(genetic_results$status=='Affected'),]$pedigree_id))
)

cat ('Number of distinct mutations:', length(unique(gr[which( gr$effect_type=='Possible mutation'), 'hgvsp'])), '\n')

cat('num of solved individuals (<18):', solved.individuals.count <- (nrow(unique(gr[which(gr$patient_age<18&gr$status=='Affected' & !is.na(gr$pedigree_gene)),c('hos_num','pedigree_gene')]))),'\n' )
cat('num of solved families (<18):', solved.families.count <- (nrow(unique(gr[which(gr$patient_age<18&gr$status=='Affected' & !is.na(gr$pedigree_gene)),c('pedigree_id','pedigree_gene')]))),'\n' )
cat('num of tested families (<18):', all.families.count <- ( length(unique(gr[which(gr$patient_age<18&gr$event_type=='OphInGeneticresults'),'pedigree_id'])) ),'\n')


cat ('Number of distinct mutations in under 18:', length(unique(gr[which(gr$patient_age<18 & gr$effect_type=='Possible mutation'), 'hgvsp'])), '\n') 

cat('num of individuals see in clinic since 2017:', length(unique(gr[which(year(gr$AppointmentDate)>=2017 &gr$status=='Affected'&gr$event_type=='OphInGeneticresults'),c('hos_num')])),'\n' )
cat ('Number of distinct mutations since 2017:', length(unique(gr[which( year(gr$AppointmentDate)>=2017
 & gr$effect_type=='Possible mutation'), 'hgvsp'])), '\n') 

affected.per.family <- sort(table(unique(genetic_results[which(genetic_results$status=='Affected'),c('hos_num','status','pedigree_id')])[,c('pedigree_id')]),decreasing=T)
cat('Number of families with more than one affected individual:',length(which(affected.per.family>1)),'\n')

cat("Number of pedigrees per individual, should not be more than 2", head(sort(table(unique(genetic_results[,c('hos_num','pedigree_id')])$hos_num),decreasing=T)), "\n")

cat('Number of patients:',length(unique(genetic_results$hos_num)),'\n')
cat('Number of families:',length(unique(genetic_results$pedigree_id)),'\n')
cat('Number of genes:',length(unique(genetic_results$gene_names)),'\n')

cat('Number of affected patients:',length(affected.patients <- unique(genetic_results[which(genetic_results$status=='Affected'),'hos_num'])), '\n')
cat('Number of affected females/males/unknown:',table(unique(genetic_results[which(genetic_results$status=='Affected'),c('hos_num','gender')])$gender), '\n')
cat('Number of affected families:',length(affected.families <- unique(genetic_results[which(genetic_results$status=='Affected'),'pedigree_id'])), '\n')


cat('Number of affected patients with genetic screening available:',

table(unique(genetic_results[,c('hos_num','method')])$method)

, '\n')

cat('Number of affected solved patients:', length(solved.patients <- unique(genetic_results[which(genetic_results$status=='Affected' & !is.na(genetic_results$pedigree_gene_names)),'hos_num'])), '\n')
cat('Number of affected solved families:', length(solved.families <- unique(genetic_results[which(genetic_results$status=='Affected' & !is.na(genetic_results$pedigree_gene_names)),'pedigree_id'])), '\n')

solved.families.gene_names <- unique(genetic_results[which(genetic_results$status=='Affected' & !is.na(genetic_results$pedigree_gene_names)),c('pedigree_id','pedigree_gene_names')])

cat('Number of genes in affected solved families:', length(unique(solved.families.gene_names$pedigree_gene_names)), '\n')


cat('Solve rate:',round(100*(length(solved.patients)/length(affected.patients))),'\n')

# some patients have 2 pedigrees (one may be solved while the other is unsolved, so they should not count as unsolved
cat('Number of affected unsolved patients:', length(unsolved.patients <- setdiff(unique(genetic_results[which(genetic_results$status=='Affected' & is.na(genetic_results$gene_names)),'hos_num']),solved.patients)), '\n')
cat('Number of affected unsolved families:', length(unsolved.families <- setdiff(unique(genetic_results[which(genetic_results$status=='Affected' & is.na(genetic_results$gene_names)),'pedigree_id']),solved.families)), '\n')


cat('Ten most common genes in solved families\n')
print( head( round(100*prop.table(sort( table(unique(genetic_results[genetic_results$pedigree_id %in% solved.families,c('pedigree_id','pedigree_gene')])$pedigree_gene), decreasing=TRUE))),10) )


cat('Number of unique disorders in affected families',length(disorder_terms <- na.omit(unique(unique(genetic_results[genetic_results$pedigree_id %in% affected.families,c('pedigree_id','disorder_term')])$disorder_term))),'\n')


cat('Ten most common disorders in affected families\n')
print( head( round(100*prop.table(sort( table(unique(genetic_results[genetic_results$pedigree_id %in% affected.families,c('pedigree_id','disorder_term')])$disorder_term), decreasing=TRUE)))) )


cat('Number of unique disorders in solved families',length(disorder_terms <- na.omit(unique(unique(genetic_results[genetic_results$pedigree_id %in% solved.families,c('pedigree_id','disorder_term')])$disorder_term))),'\n')


cat('Ten most common disorders in solved families\n')
print( head( round(100*prop.table(sort( table(unique(genetic_results[genetic_results$pedigree_id %in% solved.families,c('pedigree_id','disorder_term')])$disorder_term), decreasing=TRUE)))) )



genetic_results_summary <- unique(genetic_results[,c('hos_num','pedigree_id','status','pedigree_gene','gene','disorder_term')])
abca4.disorders <- unique(genetic_results_summary[which(genetic_results_summary$pedigree_gene=='ABCA4'),c('hos_num','disorder_term')])
cat('ABCA4 cases with retinitis pigmentosa\n')
print(abca4.disorders[grep('pigment',abca4.disorders$disorder_term,ignore.case=TRUE),])

cat('Number of cases screened by Manchester NGS panel',length(unique(genetic_results[ which(genetic_results$method=='Manchester NGS Retinal Panel'&genetic_results$status=='Affected'), 'hos_num'])),'\n')

cat('Most common mutations in solved families\n')

X <- unique(all.mutations[ which(all.mutations$status=='Affected'&all.mutations$effect_type=='Possible mutation' & !is.na(all.mutations$pedigree_gene)), c('pedigree_id','hgvsp','hgvsc')])
head(sort(table(paste(X$hgvsp,X$hgvsc)),decreasing=TRUE),30)

cat('Ten most common disorders in unsolved families\n')
print( head( round(100*prop.table(sort( table(unique(genetic_results[genetic_results$pedigree_id %in% unsolved.families,c('pedigree_id','disorder_term')])$disorder_term), decreasing=TRUE))),10) )

cat('Ten most common modes of inheritance in unsolved families\n')
print( head( round(100*prop.table(sort( table(unique(genetic_results[genetic_results$pedigree_id %in% unsolved.families,c('pedigree_id','pedigree_inheritance')])$pedigree_inheritance), decreasing=TRUE)))) )


cat('Most common mutations in unsolved families\n')

X <- unique(all.mutations[ which(all.mutations$status=='Affected'&all.mutations$effect_type=='Possible mutation' & is.na(all.mutations$pedigree_gene)), c('pedigree_id','gene','hgvsp','hgvsc')])
write.table(head(sort(table(paste(X$gene, X$hgvsp,X$hgvsc)),decreasing=TRUE),30),file='',quote=FALSE, row.names=FALSE, sep='\t')


X <- unique(all.mutations[ which(all.mutations$status=='Affected'&all.mutations$effect_type=='Possible mutation' & is.na(all.mutations$pedigree_gene)), c('pedigree_id','gene','hgvsp')])
write.table(head(sort(table(paste(X$gene, X$hgvs)),decreasing=TRUE),100),file='',quote=FALSE, row.names=FALSE, sep='\t')


cat('Possible mutations in unsolved families\n')

genetic_results[ which(genetic_results$pedigree_id %in% unsolved.families$pedigree_id), 'hgvsp']

genetic_results$status=='Affected'


range(year(genetic_results$DateRegistered))

cat('Comorbities?\n')
print(sort(table(genetic_results$ServiceName)))


table(genetic_results$EthnicityName)
