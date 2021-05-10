
library(data.table)
library(RODBC)
library(Hmisc)
library(eeptools)

myconn <- odbcDriverConnect('dsn=mehqv')

genetic_results <- sqlQuery(myconn, " select *
from [MIS_OE13].[ext13].[genetics_patient] gp
inner join [DW_LIVE].[patient].[dimpatient] dp on hos_num=dp.PatKey
inner join  [mis_oe13].[ext13].[patient] p on p.id=gp.patient_id  
inner join [MIS_OE13].[ext13].[genetics_patient_pedigree] gpp on gpp.patient_id=gp.id
left join [MIS_OE13].[ext13].[genetics_patient_diagnosis] gpd on gpd.patient_id=gp.id
left join [MIS_OE13].[ext13].[genetics_patient_relationship] gpr on gp.id=gpr.patient_id
left join [MIS_OE13].[ext13].[episode] ep on ep.patient_id=gp.patient_id 
left join [MIS_OE13].[ext13].[event] ev on ep.id=ev.episode_id
left join [MIS_OE13].[ext13].[et_ophingeneticresults_test] res on res.event_id=ev.id
left join [MIS_OE13].[ext13].[et_ophindnasample_sample] sample on sample.event_id=ev.id
where hos_num!='' and hos_num is not null
")

print(head(unsolved.patients.comments <- unique(genetic_results[ which(genetic_results$PatKey %in% unsolved.patients), c('PatKey','gene_id','effect_type','base_change','amino_acid_change','comments.1')])))

write.csv(unsolved.patients.comments, file='unsolved_patients_comments.csv',row.names=FALSE)

