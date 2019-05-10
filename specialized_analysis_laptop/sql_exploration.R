
library(SRAdb)
#sqlfile <- file.path(system.file('extdata', package='SRAdb'), 'SRAmetadb_demo.sqlite')
sqlfile <- getSRAdbFile()

sra_con <- dbConnect(SQLite(),sqlfile)

conversion <- sraConvert( c('SRR1182504','SRR1182500'), sra_con = sra_con )
