
import pandas
import pyreadr

home_dir = '/home/stephen/Documents/classes/bme/237/final_project'
exp_gdc = pyreadr.read_r(home_dir+ "/data/exp_gdc.rds")
exp_kal = pyreadr.read_r(home_dir+ "/data/exp_kal.rds")
meta_data = pyreadr.read_r(home_dir+ "/data/meta_data.rds")

exp_gdc = exp_gdc[None]
exp_kal = exp_kal[None]
meta_data = meta_data[None]

# convert exp_gdc to gene names
gdc_recode = pyreadr.read_r(home_dir+ "/data/gdc_recode.rds")
gdc_recode = gdc_recode[None]
gdc_recode = dict(set((zip(gdc_recode['ensembl_gene_id'], gdc_recode['hgnc_symbol']))))
gdc_recode_rows = list(map(gdc_recode.get, [x.split('.')[0] for x in exp_gdc.index]))
exp_gdc.index = gdc_recode_rows
exp_gdc = exp_gdc.loc[exp_gdc.index.notnull(),:]

# convert exp_kal to gene names
kal_recode = pyreadr.read_r(home_dir+ "/data/kal_recode.rds")
kal_recode = kal_recode[None]
kal_recode = dict(set((zip(kal_recode['ensembl_transcript_id'], kal_recode['hgnc_symbol']))))
kal_recode_rows = list(map(kal_recode.get, [x.split('.')[0] for x in exp_kal.index]))
exp_kal.index = kal_recode_rows
exp_kal = exp_kal.loc[exp_kal.index.notnull(),:]

# save to tsv
exp_gdc.to_csv(home_dir+ "/data/exp_gdc.csv", index=True)
exp_kal.to_csv(home_dir+ "/data/exp_kal.csv", index=True)






