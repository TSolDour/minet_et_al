ESs = read.table('editing_sites_Shoshan2021.txt', header = TRUE, sep = '\t', quote = '')

library(dplyr)
library(vroom)
library(stringr)

AG = c('A', 'G')
CT = c('C', 'T')

for(i in list.files(pattern = '_pile.txt')){
	
	sample_id = gsub('_pile.txt', '', i)
	mpile = vroom(file = i, delim = '\t', col_names = FALSE, quote = '')
	colnames(mpile) = c('ref_seq', 'position', 'ref_base', paste0('cov_', sample_id), 'read_bases', paste0('read_quals_', sample_id), paste0('align_quals_', sample_id))
	
	d = left_join(ESs, mpile, by = c('Trinity_id' = 'ref_seq', 'Editing_location' = 'position'))
	d = d[which(!is.na(d[, paste0('cov_', sample_id)])), ] # NAs exist when not a single read aligns to any part of a transcript. Thus it's not in the mpileup file at all.
	
	d$read_bases = gsub(',', '.', d$read_bases, fixed = TRUE) # make all matches the same regardless of to which strand the reads aligned
	d$read_bases = str_replace_all(d$read_bases, '[acgt]', toupper) # make all bases uppercase regardless of to which strand the reads aligned
	d$read_bases = str_replace_all(d$read_bases, fixed('.'), d$ref_base) # convert matches to the nucleotide letter
	
	
	d[, paste0('EL_', sample_id)] = apply(d, 1, function(j){
		
		tmp = table(strsplit(j['read_bases'], '')$read_bases) # make a counts table of each nucleotide
		
		if(j[['ref_base']] %in% AG){
			if('G' %in% names(tmp)) EL = tmp['G'] / sum(tmp[AG], na.rm = TRUE) else EL = 0 # if transcript CDS is on the top strand, calculate EL with G/(A+G)
		}
		if(j[['ref_base']] %in% CT){
			if('C' %in% names(tmp)) EL = tmp['C'] / sum(tmp[CT], na.rm = TRUE) else EL = 0 # if transcript CDS is on the bottom strand, calculate EL with C/(T+C)
		}
		return(as.numeric(EL))
	})
	
	d = d[, c('Trinity_id', 'Editing_location', paste0('EL_', sample_id), paste0('cov_', sample_id))] # return just the necessary columns
	
	ESs = left_join(ESs, d, by = c('Trinity_id', 'Editing_location'))
}

write.csv(ESs, 'Soff_head_ESs.csv')
