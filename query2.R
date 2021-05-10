
library(data.table)
library(RODBC)
library(Hmisc)
library(eeptools)

myconn <- odbcDriverConnect('dsn=mehqv')

genes <- read('genes.csv')
rownames(genes) <- as.character(genes$gene_id)

effects <- read('effects.csv')
rownames(effects) <- as.character(effects$id)

inheritance <- list("1"="Autosomal Dominant","2"="Autosomal Dominant with Reduced Penetrance","3"="Autosomal Recessive","4"="Mitochondrial","5"="Simplex","6"="Unknown/other","7"="X-linked Recessive")

pas <- as.data.frame(fread('PAS.csv'))
#pas <- sqlQuery(myconn, 'select * from 892039' )

print(head(methods <- sqlQuery(myconn, "select * from [MIS_OE13].[ext13].[ophingeneticresults_test_method]")))
rownames(methods) <- methods$id
print(head(genetics_patient <- sqlQuery(myconn, "select * from [MIS_OE13].[ext13].[genetics_patient]")))
print(head(genetics_patient_diagnosis <- sqlQuery(myconn, "select * from [MIS_OE13].[ext13].[genetics_patient_diagnosis]")))
print(head(genetics_patient_pedigree <- sqlQuery(myconn, "select * from [MIS_OE13].[ext13].[genetics_patient_pedigree]")))
print(head(genetics_patient_relationship <- sqlQuery(myconn, "select * from [MIS_OE13].[ext13].[genetics_patient_relationship]")))
print(head(genetics_relationship <- sqlQuery(myconn, "select * from [MIS_OE13].[ext13].[genetics_relationship]")))
print(head(genetics_study <- sqlQuery(myconn, "select * from [MIS_OE13].[ext13].[genetics_study]")))
print(head(genetics_study_subject <- sqlQuery(myconn, "select * from [MIS_OE13].[ext13].[genetics_study_subject]")))
print(head(genetics_study_proposer <- sqlQuery(myconn, "select * from [MIS_OE13].[ext13].[genetics_study_proposer]")))
print(head(et_ophindnasample_sample <- sqlQuery(myconn, "select * from mis_oe13.ext13.et_ophindnasample_sample")))
#print(head(event <- sqlQuery(myconn, "select * from mis_oe13.ext13.event")))
#write.csv(event,file='event.csv',row.names=FALSE)
print(head(patient <- sqlQuery(myconn, "select *  from mis_oe13.ext13.patient")))
rownames(patient) <- as.character(patient$id)
#sqlQuery(myconn, "select * from mis_oe13.ext13.genetics_patient")
#sqlQuery(myconn, "select * from mis_oe13.ext13.genetics_patient_pedigree")
print(head(pedigree <- sqlQuery(myconn, "select * from mis_oe13.ext13.pedigree")))
pedigree$gene_names <- genes[as.character(pedigree$gene_id),'gene_names']
print(head(ophindnasample_sample_type <- sqlQuery(myconn, "select * from MIS_OE13.ext13.ophindnasample_sample_type")))
print(head(pedigree_status <- sqlQuery(myconn, "select * from MIS_OE13.ext13.pedigree_status")))
print(head(contact <- sqlQuery(myconn,"select * from mis_oe13.ext13.contact")))
print(head(episode <- sqlQuery(myconn,"select * from mis_oe13.ext13.episode")))


pedigree$inheritance <- inheritance[as.character(pedigree$inheritance_id)]

d <- sqlQuery(myconn, "select EV.created_date as event_created_date, S.id as SampleId, GP.id as Subject_id,PPD.pedigree_id as FamilyId ,
          PS.name as PedigreeStatus, P.hos_num as HospitalNo, P.gender, P.dob, P.ethnic_group_id, C.last_name, C.first_name,
          CONVERT(date,S.created_date) as DateTaken, ST.name as SampleType, S.volume as Volume, S.comments as comment, PD.*
          FROM mis_oe13.ext13.et_ophindnasample_sample S
          join mis_oe13.ext13.event EV on EV.id= S.event_id and EV.deleted=0
          join mis_oe13.ext13.episode EP on EP.id= EV.episode_id and EP.deleted=0
          join mis_oe13.ext13.patient P on P.id= EP.patient_id
          join mis_oe13.ext13.contact C on C.id= P.contact_id
          left join mis_oe13.ext13.genetics_patient GP on GP.patient_id= P.id
          left join mis_oe13.ext13.genetics_patient_pedigree PPD on PPD.patient_id=GP.id
          left join mis_oe13.ext13.pedigree PD on PD.id=PPD.pedigree_id
          left join MIS_OE13.ext13.ophindnasample_sample_type ST on ST.id=S.type_id
          left join MIS_OE13.ext13.pedigree_status PS on PS.id=PPD.status_id order by S.id")

print(gene.count <- sort(table(pedigree$gene_names)))
gene.count <- as.data.frame(t(t(gene.count)))[,-2]
colnames(gene.count) <- c('gene_name','pedigree_count')
#pedigree
#genetics_patient_pedigree
i <- intersect( genetics_patient$patient_id, rownames(patient) )
p <- patient[i,]

rownames(genes) <- as.character(genes$gene_id)
d$gene_names <- genes[as.character(d$gene_id),'gene_names']

aa_change_type <- sqlQuery(myconn,'select * from [MIS_OE13].[ext13].[pedigree_amino_acid_change_type]')
base_change_type <- sqlQuery(myconn,'select * from [MIS_OE13].[ext13].[pedigree_base_change_type]')

rownames(aa_change_type) <- aa_change_type$id
rownames(base_change_type) <- base_change_type$id

d$aa_change <- (aa_change_type[d[,'amino_acid_change_id'],'change'])
d$base_change <- (base_change_type[d[,'base_change_id'],'change'])

d$PedigreeStatus <- as.character(d$PedigreeStatus)
d <- d[ which(!is.na(d$gene_names) & d$PedigreeStatus=='Affected' & !is.na(d$HospitalNo)), ]
print(dim(d <- unique(d)))


d$age <- floor(age_calc(as.Date(d$dob),unit='years'))

sort(table(paste(d$gene_names, d$base_change)))
sort(table(paste(d$gene_names,d$amino_acid_change)))

'RP9, PIM1K; His137Leu'

X <- d[ which(d$base_change==''&d$amino_acid_change==''), ]

rowSums(table(format(as.Date(d$created_date),'%Y'),d$gene_names))

event_type <- sqlQuery(myconn, "select * from mis_oe13.ext13.event_type")
rownames(event_type) <- event_type$id

event$name <- event_type[ as.character(event$event_type_id), 'name']

genetic.results <- event[which(event$name=='Genetic Results' & event$deleted==0),]

genetic.results$year <- format(as.Date(genetic.results$event_date),'%Y')

table(format(as.Date(genetic.results$event_date),'%Y'))


episode$hos_num <- patient[ as.character(episode$patient_id), 'hos_num']


episode$disorder_term <- disorder[as.character(episode$disorder_id),'term']

disorder.count <- sqlQuery(myconn, "select term, count(distinct patient_id) as c from [MIS_OE13].[ext13].[disorder] disorder join [MIS_OE13].[ext13].[episode] episode on  disorder.id=episode.disorder_id group by term order by c desc")

disorders.patients <- sqlQuery(myconn, "select hos_num, term from [MIS_OE13].[ext13].[disorder] disorder join [MIS_OE13].[ext13].[episode] episode on disorder.id=episode.disorder_id join [MIS_OE13].[ext13].[patient] patient on patient.id=episode.patient_id")

write.csv(disorders.patients,file='disorder_patient.csv',row.names=FALSE)

print(head(geneticresults_test <- sqlQuery(myconn, "select * from [MIS_OE13].[ext13].[et_ophingeneticresults_test]")))

genetic_results <- sqlQuery(myconn, "select * from
                            [MIS_OE13].[ext13].[genetics_patient] gp,
                            [mis_oe13].[ext13].[patient] p,
		      [MIS_OE13].[ext13].[genetics_patient_pedigree] gpp,
		[MIS_OE13].[ext13].[et_ophingeneticresults_test] res,
                            [MIS_OE13].[ext13].[event] ev,
                            [MIS_OE13].[ext13].[episode] ep
		where
		p.id=gp.patient_id and
		res.event_id=ev.id and
                ep.id=ev.episode_id and
		ep.patient_id=gp.patient_id and
		gpp.patient_id=gp.id ;" )


genetic_results <- sqlQuery(myconn, " select * from
[MIS_OE13].[ext13].[genetics_patient] gp
inner join  [mis_oe13].[ext13].[patient] p on p.id=gp.patient_id
inner join [MIS_OE13].[ext13].[genetics_patient_pedigree] gpp on gpp.patient_id=gp.id
left join [MIS_OE13].[ext13].[genetics_patient_diagnosis] gpd on gpd.patient_id=gp.id
left join [MIS_OE13].[ext13].[genetics_patient_relationship] gpr on gp.id=gpr.patient_id
left join [MIS_OE13].[ext13].[episode] ep on ep.patient_id=gp.patient_id
left join [MIS_OE13].[ext13].[event] ev on ep.id=ev.episode_id
left join [MIS_OE13].[ext13].[et_ophingeneticresults_test] res on res.event_id=ev.id
left join [MIS_OE13].[ext13].et_ophindnasample_sample dna_sample on dna_sample.event_id=ev.id")


Y2 <- genetic_results[which(genetic_results$status=='Affected'&!is.na(genetic_results$hos_num)&genetic_results$effect_type=='Possible mutation' & !is.na(genetic_results$amino_acid_change) & genetic_results$amino_acid_change!='' & genetic_results$amino_acid_change!='N/A'),]
Y2<-data.frame(ped=Y2$pedigree_id,mehno=Y2$hos_num,mut=paste(Y2$gene_names,Y2$amino_acid_change,Y2$hom,sep=':'))
head(sort(table(unique(Y2[,c('mehno','mut')])$mut),decreasing = T))
head(sort(table(unique(Y2[,c('ped','mut')])$mut),decreasing = T))

X <- as.list(by(Y2$mut,Y2$mehno,paste,collapse=' '))
head(t(t(sort(table(unlist(X)),decreasing=T))))
X <- unique(genetic_results[,c('hos_num','gene_names')])
gene.count2 <- sort(table(X$gene_names))
gene.count2 <- as.data.frame(t(t(gene.count2)))[,-2]
colnames(gene.count2) <- c('gene_name','person_count')

Y <- unique(genetic_results[,c('hos_num','gene_names','effect_type','homo','base_change','amino_acid_change')])
Y2 <- Y[which(Y$effect_type=='Possible mutation' & !is.na(Y$amino_acid_change) & Y$amino_acid_change!='' & Y$amino_acid_change!='N/A'),]

X <- genetic_results[which(!is.na(genetic_results$hos_num)),]
length(unique(X$hos_num))
length(unique(X$pedigree_id))
table(unique(X[,c('hos_num','status')])$status)
X <- X[which(X$status=='Affected' & !is.na(X$amino_acid_change) & X$amino_acid_change!='' & X$amino_acid_change!='N/A' & X$amino_acid_change!='normal' & X$effect_type=='Possible mutation'),]
length(unique(paste(X$gene_names,X$amino_acid_change,sep=':')))
length(unique(X$gene_names))
head(t(t(sort(table(paste(X$gene_names,X$amino_acid_change,sep=':')),decreasing = T))))
X$PatKey<-X$hos_num
pas <- pas[which(pas$ConsultantFirmName %in% c("Webster Andrew","Michaelides Michel")),]
X <- merge(X,pas,by='PatKey')
X <- X[which(X$IsLastAppointment==1),]
t(t(table(substr(X$AppointmentDate,1,4))))
t(t(table(substr(X$result_date,1,4))))
method.year <- t(t(table(as.character(X$method),substr(X$result_date,1,4))))
write.csv(method.year,file='')
patient.dob <- t(t(table(substr(X$dob,1,4))))
write.csv(patient.dob,file='')

hist(2018-as.numeric(substr(X$dob,1,4)),main='age distribution of IRD patients',xlab='age')

age.at.death <- as.numeric(substr(X$date_of_death,1,4))-as.numeric(substr(X$dob,1,4))
hist(age.at.death,main='age at death distribution of IRD patients',xlab='age')


write.csv(table(X$gene_names,age.at.death),file='')

hist(as.numeric(substr(X$AppointmentDate,1,4)),main='date of last appointment',xlab='date')

X.2017 <- X[which(substr(X$AppointmentDate,1,4)=='2017'),]
write.csv(t(t(sort((table(X.2017$gene_names)),decreasing=T))),file='')

write.csv(t(t(sort(table(unique(X[,c('gene_names','amino_acid_change')])$gene_names),decreasing = T))),file='')

table(X$gender)

dim(Y2)

head(as.data.frame(t(t(sort(table(paste(Y2$gene_names,Y2$amino_acid_change,sep=':')),decreasing=T))))->xx)
write.csv(xx,file='common_mutations.csv')


effect_type.gene_names <- (t(table(Y$effect_type,Y$gene_names)))
write.csv(effect_type.gene_names,file='effect_type_by_gene.csv')

gene.counts <- merge(gene.count,gene.count2,by='gene_name')
write.csv(gene.counts,file='gene_counts2.csv',row.names=FALSE)


library(googlesheets)
k <- '1YFpQssZQykuUZVa0LK8B8vKjdo_flccxHsFsSOedV8U'
gsheet <- gs_key(k)
gs_ws_new (gsheet, input=gene.counts, trim=TRUE, verbose=TRUE, ws_title='gene_counts')
#gs_edit_cells(gsheet, input=gene.counts, trim=TRUE, verbose=TRUE)
gsheet %>% gs_read()

gs_ws_new (gsheet, input=disorder.count, trim=TRUE, verbose=TRUE, ws_title='disorder_counts')
#gs_edit_cells(gsheet, input=gene.counts, trim=TRUE, verbose=TRUE)
gsheet %>% gs_read()

gs_ws_new (gsheet, input=method.year, trim=TRUE, verbose=TRUE, ws_title='method_year')


#write.csv( , file='', row.names=FALSE)

library(plotly)

p <- plot_ly(
  x = disorder.count$term,
  y = disorder.count$c,
  name = "Disorder Count",
  type = "bar"
)

p %>% offline()

# Create a shareable link to your chart
# Set up API credentials: https://plot.ly/r/getting-started
chart_link = api_create(p, filename="bar-basic")
chart_link

library(ggplot2)
p<-ggplot(data=disorder.count, aes(x=disorder.count$term, y=disorder.count$c)) +
  geom_bar(stat="identity")
p

# Horizontal bar plot
p + coord_flip()
