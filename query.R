

library(data.table)
library(RODBC)
library(Hmisc)
library(eeptools)

myconn <- odbcDriverConnect('dsn=mehqv')

myconn <- odbcDriverConnect('dsn=DW_DW')

xx <- sqlQuery(myconn, "SELECT TABLE_NAME FROM INFORMATION_SCHEMA.TABLES")

#sqlQuery(myconn, "SELECT TABLE_NAME, COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS")


#dimpat <- sqlQuery(myconn, "select top 10 * from [DW\\DW].[patient].[dimpatient]")


d2 <- sqlQuery(myconn, "select
		EV.created_date as event_created_date,
		S.id as SampleId,
		GP.id as Subject_id,
		PPD.pedigree_id as FamilyId , 
		PS.name as PedigreeStatus,
		P.hos_num as HospitalNo,
		P.gender,
		P.dob,
		P.ethnic_group_id,
		C.last_name,
		C.first_name,
                CONVERT(date,S.created_date) as DateTaken,
		ST.name as SampleType,
		S.volume as Volume,
		S.comments as comment,
		PD.*,
		DP.*
        FROM
	mis_oe13.ext13.et_ophindnasample_sample S
	join mis_oe13.ext13.event EV on EV.id= S.event_id and EV.deleted=0
        join mis_oe13.ext13.episode EP on EP.id= EV.episode_id and EP.deleted=0
        join mis_oe13.ext13.patient P on P.id= EP.patient_id
        join mis_oe13.ext13.contact C on C.id= P.contact_id
        left join mis_oe13.ext13.genetics_patient GP on GP.patient_id=P.id
        left join mis_oe13.ext13.genetics_patient_pedigree PPD on PPD.patient_id=GP.id
	left join mis_oe13.ext13.pedigree PD on PD.id=PPD.pedigree_id
        left join MIS_OE13.ext13.ophindnasample_sample_type ST on ST.id=S.type_id
        left join MIS_OE13.ext13.pedigree_status PS on PS.id=PPD.status_id
	left join [DW_LIVE].[patient].[dimpatient] DP on DP.patkey=P.hos_num
	order by S.id") 


d3 <- merge( d2, pas, by=c('PatKey'), all.x=TRUE )


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

which(genetic.results$year=='2015')


episode$hos_num <- patient[ as.character(episode$patient_id), 'hos_num']


disorder <- sqlQuery(myconn, "select * from [MIS_OE13].[ext13].[disorder]")
rownames(disorder) <- disorder$id


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
		gpp.patient_id=gp.id " )


cols <- c("PatKey", "id", "patient_id", "patient_id.2", "patient_id.4", 
"comments", "last_modified_user_id", 'Study',
"last_modified_date", "created_user_id", "created_date", "gender_id",
"CreateDT", "SeqId", "id.1", "pas_key",
"dob", "gender", "hos_num", "nhs_num", "last_modified_user_id.1",
"last_modified_date.1", "created_user_id.1", "created_date.1",
"gp_id", "date_of_death", "practice_id", "ethnic_group_id", "contact_id",
"no_allergies_date", "no_family_history_date", "no_risks_date",
"deleted.1", "nhs_num_status_id",
"is_deceased","dod"
"ExtractDT.x",
"pedigree_id", "status_id",
'term',"disorder_id"
"related_to_id", "relationship_id", "last_modified_user_id.4",
"firm_id", "start_date", "episode_status_id", "legacy", "deleted.2",
"ExtractDT.1", "SeqId.5", "id.6", "episode_id",
"event_type_id", "last_modified_user_id.6", "last_modified_date.6",
"created_date.6", "event_date", "deleted.3", "ExtractDT.2", "SeqId.6",
"id.7", "event_id", "gene_id", "method_id",
"comments.1",
"exon", "base_change", "amino_acid_change", "assay", "effect_id", "homo",
"result", "result_date", "last_modified_user_id.7", 
"created_user_id.7", "created_date.7", "base_change_id", "amino_acid_change_id",
"genomic_coordinate", "genome_version", "gene_transcript",
"SeqId.7", "effect_type", "event_year", "gene_names", "status",
"method", "AgeOnAppointmentKeyN", "AppointmentDate", "AppointmentTime",
"BookingDateKeyN", "BookingTime", "FirstBookingDateKeyN", "FirstBookingTime", "CancellationDate",
"CancelledBy", "AttendanceTime", "DepartureTime", "SeenTime",
"NewFUP", "PatientCategoryKeyN", "PatientCategoryCode", "PatientCategoryName",
"ClinicKeyN", "ClinicCode", "ClinicName", "SessionElementKeyN",
"SessionElementCode", "SessionElementName", "ClinicLocationKeyN",
"ClinicLocationCode", "ClinicLocationName", "SiteKeyN", "SiteCode",
"SiteName", "ServiceKeyN", "ServiceCode", "ServiceName", "ManagerStaffKeyN",
"ManagerStaffCode", "ManagerStaffName", "SeeingStaffKeyN", "SeeingStaffCode",
"SeeingStaffName", "ConsultantFirmKeyN", "ConsultantFirmName",
"ChargeableKeyN", "ChargeableName", "IsChargeable", "CommissionerKeyN",
"CommissionerCode", "CommissionerName", "CCGofResidenceKeyN",
"CCGofResidenceCode", "CCGofResidenceName", "TreatmentFunctionCode",
"AppointmentStatusCode", "SuspensionReasonCode", "SuspensionReasonName",
"SuspensionReasonKeyN", "AppointmentReasonCode", "AppointmentReasonName",
"AppointmentReasonKeyN", "SuspensionStartDate", "DirectorateKeyN",
"DirectorateName", "CancelledByName", "NewFUPName", "AppointmentStatusName",
"IsAttendance", "IsDNA", "IsWaiting", "IsCancelled", "IsFuture",
"IsNoOutcome", "PAS_appointmentKey", "JourneyTimeMI", "PAS_EP_NO",
"PAS_UBRN", "IsConsultantLed", "ReferralNumber", "EpisodeNumber",
"AppointmentNumber", "RefSeqId", "EpSeqId", "EpisodeDischargeDateKeyN",
"IsDischarged", "IsValid", "IsExcluded", "RuleID", "IsLastAttendance",
"IsLastAppointment", "IntendedFollowUpDate", "headerGrp", "AppSeqId")

colnames(genetic_results) <- gsub('\\.y$','',colnames(genetic_results))
genetic_results <- genetic_results[,intersect(cols,colnames(genetic_results))]


length(unique(genetic_results$hos_num))
length(unique(genetic_results$pedigree_id))

table(unique(genetic_results[,c('hos_num','status')])$status)

#genetic_results <- genetic_results[-which(is.na(genetic_results$hos_num)),]

genetic_results[which(genetic_results$pedigree_id==224),'hos_num']


head(table(unique(genetic_results[,c('hos_num','status','pedigree_id')])[,c('pedigree_id','status')])

# variants in unsolved patients

unsolved.mutations <- unique(genetic_results[ which(genetic_results$status=='Affected' & is.na(genetic_results$gene_names)),c('hos_num','status','effect_type','amino_acid_change','homo')])

cat('Total number of mutations in unsolved patients:',length(unique(unsolved.mutations$amino_acid_change)),'\n')

write.csv(unique(unsolved.mutations[c('amino_acid_change','homo')]),file='',row.names=F)

(sort(table(unsolved.mutations$amino_acid_change),decreasing=TRUE))


edd <- read('EDD_appointments_PAS_numbers_1991_to_2019_AC_03_Jul_2019.csv')[,c('PAS No','FILEDUNDER','DATESEEN')]
colnames(edd) <- c('PatKey','edd.filed_under','edd.dateseen')

genetic_results.edd <-merge(genetic_results,edd,by='PatKey',all.x=TRUE)

edd.genes <- unique(genetic_results.edd[ which(!is.na(genetic_results.edd$edd.filed_under)&!is.na(genetic_results.edd$gene_names)), c('hos_num','gene_names') ])

write.csv( t(t(sort(table(edd.genes$gene_names),decreasing=TRUE))), file='')



t( t( table(unique(cbind(genetic_results2$PatKey,year(genetic_results2$AppointmentDate)))[,2]) ) )

range(genetic_results$AppointmentDate,na.rm=TRUE)

genetic_results.affected.possible_mutation <- genetic_results[which(genetic_results$status=='Affected'&!is.na(genetic_results$hos_num)&genetic_results$effect_type=='Possible mutation' & !is.na(genetic_results$amino_acid_change) & genetic_results$amino_acid_change!='' & genetic_results$amino_acid_change!='N/A'),]

length(unique(genetic_results.affected.possible_mutation$hos_num))


Y2<-data.frame(ped=genetic_results.affected.possible_mutation$pedigree_id,mehno=genetic_results.affected.possible_mutation$hos_num,mut=paste(genetic_results.affected.possible_mutation$gene_names,genetic_results.affected.possible_mutation$amino_acid_change,genetic_results.affected.possible_mutation$hom,sep=':'))
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
#X <- X[which(X$status=='Affected' & !is.na(X$amino_acid_change) & X$amino_acid_change!='' & X$amino_acid_change!='N/A' & X$amino_acid_change!='normal' & X$effect_type=='Possible mutation'),]
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
k <- ''
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

prop.table(sort(table(unique(genetic_results[,c('PatKey','gene_names')])$gene_names)))


genetic_results$AppointmentDateYear <- year(genetic_results$AppointmentDate)

XX <- table( unique(genetic_results[,c('PatKey','gene_names','AppointmentDateYear')])[,c('gene_names','AppointmentDateYear')] )[,as.character(2008:2018)]


XX[order(rowSums(XX),decreasing=TRUE),]

length( unique( genetic_results[ which(genetic_results$AppointmentDateYear==2018), 'PatKey' ]) )

xxx <- t( t(sort(table(unique(genetic_results[ which(genetic_results$AppointmentDateYear==2018), c('PatKey','gene_names') ])$gene_names),decreasing=TRUE)))
rownames(xxx) <- gsub(';.*','',rownames(xxx))
rownames(xxx) <- gsub(',.*','',rownames(xxx))
write.csv(xxx , file='' )





