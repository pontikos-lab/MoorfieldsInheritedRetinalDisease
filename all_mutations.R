library(stringr)

if (!file.exists('all_mutations.csv')) {
cat('Creating file all_mutations.csv...\n')
X <- unique(genetic_results[,c('hos_num','pedigree_id','status','study_name','study_id','method','assay','effect_type','gene','hgvsp','hgvsc','homo','pedigree_gene','disorder_term')])
mutations <- rbindlist(by(X,X$hos_num,function(x) {
data.frame(hos_num=unique(x$hos_num),
status=unique(x$status),
pedigree_id=paste(unique(na.omit(x$pedigree_id)),collapse=';'),
study_name=paste(unique(na.omit(x$study_name)),collapse=';'),
study_id=paste(unique(na.omit(x$study_id)),collapse=';'),
method=paste(unique(na.omit(x$method)),collapse=';'),
assay=paste(unique(na.omit(x$assay)),collapse=';'),
effect_type=paste(unique(na.omit(x$effect_type)),collapse=';'),
gene=paste(unique(na.omit(x$gene)),collapse=';'),
hgvsp=paste(unique(na.omit(x$hgvsp)),collapse=';'),
hgvsc=paste(unique(na.omit(x$hgvsc)),collapse=';'),
homo=paste(unique(na.omit(x$homo)),collapse=';'),
pedigree_gene=paste(unique(na.omit(x$pedigree_gene)),collapse=';'),
disorder_term=paste(unique(na.omit(x$disorder_term)),collapse=';')
) }))
dim(all.mutations<-X[ -which(is.na(X$hgvsc)&is.na(X$hgvsp)),])
write.csv(all.mutations,file='all_mutations.csv',row.names=FALSE)
} else {
cat('File all_mutations.csv exists...\n')
print(dim(all.mutations <- read('all_mutations.csv')))
}

all.mutations$hgvsp.old <- all.mutations$hgvsp
all.mutations$hgvsc.old <- all.mutations$hgvsc

all.mutations$hgvsp <- gsub(';','',all.mutations$hgvsp)
all.mutations$hgvsp <- gsub('^\\s+;','',gsub('NA:','',all.mutations$hgvsp))
all.mutations$hgvsc <- gsub('^\\s+;','',gsub('NA:','',all.mutations$hgvsc))
all.mutations$hgvsp[which(all.mutations$hgvsp==';')] <- ''
all.mutations$hgvsc[which(all.mutations$hgvsc==';')] <- ''
all.mutations$hgvsp <- gsub('.*:','', all.mutations$hgvsp)
all.mutations$hgvsc <- gsub('.*:','', all.mutations$hgvsc)
all.mutations$hgvsp[which(is.na(all.mutations$hgvsp))] <- ''
all.mutations$hgvsc[which(is.na(all.mutations$hgvsc))] <- ''
all.mutations$hgvsp <- gsub('&#8727','\\*',all.mutations$hgvsp)
all.mutations$hgvsc <- gsub('&#8727','\\*',all.mutations$hgvsc)
all.mutations$hgvsp <- gsub('^P\\.','p\\.', all.mutations$hgvsp)
all.mutations$hgvsc <- gsub('^P\\.','p\\.', all.mutations$hgvsc)
all.mutations$hgvsc <- gsub('\\[','',all.mutations$hgvsc)
all.mutations$hgvsc <- gsub('\\]','',all.mutations$hgvsc)
all.mutations$hgvsp <- gsub('\\[','',all.mutations$hgvsp)
all.mutations$hgvsp <- gsub('\\]','',all.mutations$hgvsp)
all.mutations$hgvsc <- gsub('^(\\d+[ACGT]+>[ACGT]+)','c\\.\\1',all.mutations$hgvsc)
all.mutations$hgvsc <- gsub('^c(\\d+[ACGT]+>[ACGT]+)','c\\.\\1',all.mutations$hgvsc)
all.mutations$hgvsp <- gsub('Het','',all.mutations$hgvsp)


all.mutations$hgvsp <- gsub('p?.?SPLICE','',all.mutations$hgvsp)
all.mutations$hgvsp <- gsub('p?.?acceptor','',all.mutations$hgvsp)

amino.acids <- rbind(
data.frame(full.name='Alanine',code='Ala',code2='A'),
data.frame(full.name='Arginine',code='Arg',code2='R'),
data.frame(full.name='Asparagine',code='Asn',code2='N'),
data.frame(full.name="Aspartic acid",code='Asp',code2='D'),
data.frame(full.name="Asparagine or aspartic acid",code='Asx',code2='B'),
data.frame(full.name='Cysteine', code='Cys', code2='C'),
data.frame(full.name='Glutamic acid',code='Glu',code2='E'),
data.frame(full.name='Glutamine',code='Gln',code2='Q'),
data.frame(full.name='Glutamine or glutamic acid', code='Glx',code2='Z'),
data.frame(full.name='Glycine', code='Gly', code2='G'),
data.frame(full.name='Histidine', code='His', code2='H'),
data.frame(full.name='Isoleucine', code='Ile', code2='I'),
data.frame(full.name='Leucine', code='Leu',code2='L'),
data.frame(full.name='Lysine', code='Lys',code2='K'),
data.frame(full.name='Methionine', code='Met', code2='M'),
data.frame(full.name='Phenylalanine',code='Phe',code2='F'),
data.frame(full.name='Proline', code='Pro',code2='P'),
data.frame(full.name='Serine', code='Ser', code2='S'),
data.frame(full.name='Threonine', code='Thr', code2='T'),
data.frame(full.name='Tryptophan', code='Trp', code2='W'),
data.frame(full.name='Tyrosine', code='Tyr', code2='Y'),
data.frame(full.name='Valine', code='Val', code2='V'),
data.frame(full.name='Stop', code='Ter', code2='X')
)
rownames(amino.acids) <- amino.acids$code2

 #G863A
#H423P
#all.mutations$hgvsp <- gsub('P1401P', 'p.Pro1401=', all.mutations$hgvsp)
#E1033Q
#A52S
#D90H
#L782H
#P1948L
#D2095D
#K192E
#P1380L
#p.C1490Yp.L1850P

#grep('p.291P>LP',all.mutations$hgvsp)

#p.Leu156Pro/Arg
#p.R2107H/P
#p.Lys1255ArgfsX8PCR/diges
#p.Leu156Pro/Arg
#p.Leu156Pro/Arg
#p.del/FS
#p.W1408R/R1640W
#p.I1120L/V
#p.G991R/X
#p.V256V/
#p.P1380LHet/p.D1532NHet
#p.Leu156Pro/Arg
#p.Leu156Pro/Arg
#p.n/a
#p.Leu156Pro/Arg

all.mutations$hgvsp <- gsub('X','Ter', all.mutations$hgvsp)

all.mutations$hgvsp <- gsub('p.thr383Ilefs*13','p.Thr383Ilefs*13', all.mutations$hgvsp, fixed=TRUE)

all.mutations$hgvsp <- gsub('p.Glu767SerfsX21','p.Glu767SerfsTer21',all.mutations$hgvsp)

all.mutations$hgvsp <- gsub('p.Arg677X','p.Arg677Ter',all.mutations$hgvsp)

# EFEMP1
all.mutations$hgvsp <- gsub('R345W','p.Arg345Trp', all.mutations$hgvsp)
all.mutations[ grep('R345W',all.mutations$hgvsc), 'hgvsp'] <- 'p.Arg345Trp'
all.mutations[ grep('c.1033C>T',all.mutations$hgvsp), 'hgvsp'] <- 'p.Arg345Trp'
all.mutations[ grep('p.Arg345Trp',all.mutations$hgvsp), 'hgvsc'] <- 'NM_001039348.3:c.1033C>T'

all.mutations$hgvsp <- gsub('^p?.?G1961E','p.Gly1961Glu',all.mutations$hgvsp)
all.mutations$hgvsp <- gsub('p.E700Ter','p.Glu700Ter',all.mutations$hgvsp)
all.mutations[which('p.903'==all.mutations$hgvsp&all.mutations$gene=='OPA1'),'hgvsp']<-'p.Val903Glyfs'
# p.Thr383Ilefs*13: 11
# p.Thr383IlefsTer13: 10
# p.Thr383llefsTer13: 7 
all.mutations$hgvsp <- gsub('p.Thr383llefsTer13','p.Thr383IlefsTer13',all.mutations$hgvsp)
all.mutations$hgvsp <- gsub('p?.?A269fsTer1','Ala269fsTer1',all.mutations$hgvsp)
all.mutations$hgvsp <- gsub('\\*','Ter',all.mutations$hgvsp)


all.mutations[grep('p.Leu2027Phe',all.mutations$hgvsp),'hgvsc'] <- ':c.6079C>T'

# ABCA4
all.mutations[grep('c.5461-10T>C',all.mutations$hgvsc),'hgvsc'] <- 'NM_000350.2:c.5461-10T>C'
all.mutations[intersect(grep('ABCA4',all.mutations$gene),grep('IVS-10T>C',all.mutations$hgvsc)),'hgvsc'] <- 'NM_000350.2:c.5461-10T>C'
all.mutations[which(all.mutations$hgvsc=='NM_000350.2:c.5461-10T>C'),'hgvsp'] <- ''
all.mutations[grep('c.5714+5G>A',all.mutations$hgvsc),'hgvsc']<-'NM_000350.3:c.5714+5G>A'
all.mutations[grep('c.5714+5G>A',all.mutations$hgvsc),'hgvsp']<-''
all.mutations[grep('c.634C>T',all.mutations$hgvsc),'hgvsc'] <- 'NM_000350.3:c.634C>T'
all.mutations[grep('c.859-9T>C',all.mutations$hgvsc),'hgvsc'] <- 'NM_000350.3:c.859-9T>C'
all.mutations[grep('c.859-9T>C',all.mutations$hgvsc),'hgvsp'] <- ''


# ABCC6
all.mutations[grep('c.3634-3C>A',all.mutations$hgvsc),'hgvsc'] <- 'NM_001171.5:c.3634-3C>A'
all.mutations[grep('c.3634-3C>A',all.mutations$hgvsc),'hgvsp'] <- ''


all.mutations$hgvsp <- gsub('(^p?\\.?[ARNDBCEQZGHILKMFPSTWYV]\\d+[ARNDBCEQZGHILKMFPSTWYVX])(p.[ARNDBCEQZGHILKMFPSTWYV]\\d+[ARNDBCEQZGHILKMFPSTWYVX])$','\\1;\\2', all.mutations$hgvsp)
xx <- str_match(all.mutations$hgvsp,'^p?\\.?([ARNDBCEQZGHILKMFPSTWYV])(\\d+)([ARNDBCEQZGHILKMFPSTWYVX])$')
hgvsp2 <- as.character(paste('p.',amino.acids[xx[,2],][,2],xx[,3],amino.acids[xx[,4],][,2],sep=''))
all.mutations$hgvsp[grep('NA',hgvsp2,invert=TRUE)] <- grep('NA',hgvsp2,invert=TRUE,value=TRUE)
all.mutations$hgvsp <- paste('p.',all.mutations$hgvsp,sep='')
all.mutations$hgvsp <- gsub('^p.p.','p.',all.mutations$hgvsp)
all.mutations$hgvsp <- gsub('^p.p.','p.',all.mutations$hgvsp)
all.mutations$hgvsp <- gsub('^p.$','',all.mutations$hgvsp)
all.mutations$hgvsp <- gsub('\t','', all.mutations$hgvsp)
all.mutations$hgvsp <- gsub( 'p.none', '', all.mutations$hgvsp)
all.mutations$hgvsp <- gsub( 'p.normal', '', all.mutations$hgvsp)


#EYS
all.mutations[grep('c.5928-2A>G',all.mutations$hgvsc),'hgvsc']<- 'NM_001142800.1:c.5928-2A>G'
all.mutations[grep('c.5928-2A>G',all.mutations$hgvsc),'hgvsp']<- ''


clean.mut <- function(all.mutations, split_char) {
	new.mutations <- data.frame()
	lines_to_split <- unique(c(grep( split_char, all.mutations$hgvsp),grep(split_char,all.mutations$hgvsc)))
	cat('Number of lines to split:', length(lines_to_split), '\n')
	for (i in lines_to_split) {
	x.mut <- all.mutations[i,]
	x.mut$hgvsp <- gsub(sprintf('^\\s+%s',split_char),'',gsub('NA:','',gsub(sprintf('%s:',x.mut$gene),'',x.mut$hgvsp)))
	# hgvsp
	if (grepl(split_char ,x.mut[['hgvsp']]) && nchar(x.mut[['hgvsp']])>1) {
		hgvsp <- unlist(strsplit(x.mut[['hgvsp']],split_char))
		hgvsp.1 <- hgvsp[[1]]
		if (length(hgvsp)>1) {
			hgvsp.2 <- hgvsp[[2]]
		}
		all.mutations[i,'hgvsp'] <- hgvsp.1
		x.mut$hgvsp <- hgvsp.1
	}
	# hgvsc
	if (grepl(split_char ,x.mut[['hgvsc']]) && nchar(x.mut[['hgvsc']])>1) {
		hgvsc <- unlist(strsplit(x.mut[['hgvsc']],split_char))
		hgvsc.1 <- hgvsc[[1]]
		if (length(hgvsc)>1) {
			hgvsc.2 <- hgvsc[[2]]
		}
		all.mutations[i,'hgvsc'] <- hgvsc.1
		x.mut$hgvsc <- hgvsc.1
	}
	new.mutations <- rbind(new.mutations,x.mut)
	}
	all.mutations <- rbind( all.mutations, new.mutations)
	return(all.mutations)
}
cat('Splitting ";" \n')
print(dim(all.mutations <- clean.mut(all.mutations, ';')))
cat('Splitting "," \n')
print(dim(all.mutations <- clean.mut(all.mutations, ',')))
#cat('Splitting "/" \n')
#print(dim(all.mutations <- clean.mut(all.mutations, '/')))

#sort(unique(paste(all.mutations$gene,all.mutations$hgvsp,sep=':')))
#head( sort(table( paste(all.mutations$gene, all.mutations$hgvsp , sep=':')), decreasing=TRUE) , 100 )

#all.mutations[i,]

#grep('ABCA4.*1961.*',all.mutations$hgvsp)

#ABCA4:p.Gly1961Glu
#ABCA4:G1961E
#ABCA4:p.G1961E

all.mutations$hgvsp.new <- all.mutations$hgvsp
all.mutations$hgvsc.new <- all.mutations$hgvsc

all.mutations$hgvsp <- all.mutations$hgvsp.old
all.mutations$hgvsc <- all.mutations$hgvsc.old

all.mutations <- all.mutations[-which(all.mutations$hgvsp.new==''&all.mutations$hgvsc.new==''),]

print(dim(all.mutations))

write.csv(all.mutations,file='all_mutations_clean.csv',row.names=FALSE)

write.csv(unique(paste(all.mutations$gene,all.mutations$hgvsp.new,sep=':')), file='all_mutations_clean_hgvsp.csv', row.names=FALSE, quote=FALSE)


