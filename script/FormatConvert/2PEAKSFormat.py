
#usage: python 2PEAKSex.py result-location name type
# python 2PEAKSex.py ./ M1 DIA-NN
# python 2PEAKSex.py ./ Sm0 pFind
# python 2PEAKSex.py ./ Z MaxQuant

import pandas as pd
import os
import warnings
import argparse

warnings.filterwarnings('ignore')

pattern=r"\b[a-zA-Z]"
pd.set_option('display.unicode.east_asian_width',True)

def doPF(input_path, name):
	pd.set_option('display.unicode.east_asian_width',True)
	dff=pd.read_excel(input_path,header=1,usecols=['ID','AC','No.Peptide','Description']) 
	dff=dff.dropna()
	delete_string='REV_'
	dff=dff[~dff['AC'].str.contains(delete_string)]
	dff.columns=['Protein Group','Accession','#Peptides','Description']
	dff['Gene']=''
	dff['GN']='0'
	order=['Protein Group','Accession','Gene','#Peptides','Description','GN']
	dff=dff[order]
	dff['Description'] = dff['Description'].str.replace('>', '')
	dff['Description'] = dff['Description'].str.replace('CK820\s', 'CK820_',regex=True)
	dff['GN']=dff['Accession'].str.split('_').str[0] 
	dff['GN']=dff['GN'].str.split('-').str[0] 
	dff['GN']=dff['GN'].astype(str)
	dff.loc[dff['GN'].str.contains('sp,',na=True),['GN']]=dff['Description'].str.extract(r'GN=(\w+)\sPE', expand=False)
	dff['GN']==dff['GN'].fillna('XXXX')
	dff.loc[dff['GN'].str.contains('XXXX',na=True),['GN']]=dff['Description'].str.extract(r'GN=(\w+-\w+)\sPE', expand=False)
	dff['Gene']=dff['GN']
	dff.to_excel(name+".proteins-ori.xlsx",index=False)

	df=pd.read_excel(input_path,header=2,usecols=range(0,19))  
	df=df.dropna(subset=['Sequence'])  
	delete_strings=['SubSet','SameSet']
	df=df[~df['Unnamed: 1'].str.contains('|'.join(delete_strings),na=False)]  
	df['Unnamed: 0']=df['Unnamed: 0'].fillna(method='ffill')
	df['Unnamed: 1']=df['Unnamed: 1'].fillna(method='ffill')
	delete_string='REV_'
	df=df[~df['Unnamed: 1'].str.contains(delete_string)]
	delete_string1='decoy'
	df=df[~df['Target/Decoy'].str.contains(delete_string1,na=False)]
	delete_string2='CON_'
	df=df[~df['Unnamed: 1'].str.contains(delete_string2)]
	df=df[['Unnamed: 0','Unnamed: 1','Sequence','Calc.MH+','Mass_Shift(Exp.-Calc.)','Final_Score','Modification','Positions','File_Name','Charge','Spec_Num']]
	df.columns=['Protein Group','Accession','Peptide','Mass','ppm','-10LgP','PTM','Start','Scan','z','#Spec']
	df['Length']=''
	df['Source File']=''
	df['GN']='0'
	df['Gene']=''
	df['Area']='1'
	df['Intensity']='10000000'
	order=['Protein Group','Accession','Gene','Peptide','Mass','Length','ppm','z','Area','Intensity','Source File','Scan','#Spec','Start','PTM','GN']
	df=df[order].reset_index(drop=True)
	df=df.dropna(subset=['#Spec'])
	df['#Spec']=df['#Spec'].astype(int)
	df['Source File']=df['Scan'].str.split('.').str[0]
	df['Length']=df['Peptide'].str.len()
	df2=pd.read_excel(name+".proteins-ori.xlsx",header=0)
	df['GN']=df['Accession'].map(df2.set_index('Accession')['GN'])
	df['GN']=df['GN'].str.replace('_(.+)', '')
	df=df[~df['Accession'].str.contains('CON_')]
	df['Gene']=df['GN']
	df.to_csv(name+".protein-peptides.csv",index=False)
	os.remove(name+".proteins-ori.xlsx")
	print('--convert pFind results to PEAKS format finished--')

def doMQ(input_path,name):
	pd.set_option('display.unicode.east_asian_width',True)
	dff=pd.read_excel(input_path,usecols=['Sequence','Length','Modifications','Modified sequence','Leading razor protein','Raw file','Charge','m/z','Mass','Retention time','PEP','MS/MS count','MS/MS IDs','Score','Delta score','Intensity'])
	dff.columns=['Sequence','Length','PTM','Peptide','Accession','File','z','m/z','Mass','RT','-10LgP','#Spec','Scan','Score','Delta score','Intensity']
	dff['GN']='0'
	dff['Gene']=''
	dff['Area']='1'
	order=['Accession','GN','Gene','File','Sequence','Peptide','-10LgP','Mass','Length','m/z','z','RT','Area','Intensity','Scan','#Spec','PTM','Score','Delta score']
	dff=dff[order]
	dff['GN']='0'
	dff['Intensity']=dff['Intensity'].fillna('0')
	dff['Peptide']=dff['Peptide'].str.replace('_','')
	dff['Peptide']=dff['Peptide'].str.replace('_','')
	dff['Peptide']=dff['Peptide'].str.replace('(P)','P',regex=False)
	dff['Peptide']=dff['Peptide'].str.replace('(NQ)','NQ',regex=False)
	dff['Peptide']=dff['Peptide'].str.replace('(M)','M',regex=False)
	dff['Peptide']=dff['Peptide'].str.replace('(W)','W',regex=False)
	dff['Peptide']=dff['Peptide'].str.replace('(K)','K',regex=False)
	dff['Peptide']=dff['Peptide'].str.replace('(Y)','Y',regex=False)
	dff['Peptide']=dff['Peptide'].str.replace('(STY)','STY',regex=False)
	dff['-10LgP']=dff['-10LgP'].fillna('1')
	dff['z']=dff['z'].fillna('0')
	dff['RT']=dff['RT'].fillna('100')
	dff['Score']=dff['Score'].fillna('0')
	dff=dff.dropna(subset=['#Spec']) # type: ignore
	dff['#Spec']=dff['#Spec'].astype(int)
	dff['GN']=dff['Accession'].str.split('-').str[0] 
	dff['GN']=dff['GN'].str.split('_').str[0] 
	dff['GN']=dff['GN'].astype(str)
	dff.loc[dff['GN'].str.contains('sp',na=True),['GN']]=dff['Accession'].str.extract(r'\.(\w+)\.', expand=False)
	dff.loc[dff['GN'].str.contains('tr',na=True),['GN']]=dff['Accession'].str.extract(r'\.(\w+)\.', expand=False)
	dff.loc[dff['GN'].str.contains('REV_',na=True),['GN']]=''
	dff.loc[dff['GN']=='CK820_G0031061',['GN']]='MMP20'
	dff.loc[dff['GN']=='CK820_G0035821',['GN']]='AMELX'
	dff.loc[dff['GN']=='A5JJS7',['GN']]='AMELX'
	dff.loc[dff['GN']=='CK820_G0053051',['GN']]='AMELY'
	dff.loc[dff['GN']=='B4DPP6',['GN']]='ALB'
	dff.loc[dff['GN']=='A0A2J8UTQ6',['GN']]='AMTN'
	dff.loc[dff['GN']=='H2Q4M8',['GN']]='MMP20'
	dff.loc[dff['GN']=='A0A2I3SAE8',['GN']]='AMELX'
	dff.loc[dff['GN']=='E1U7Q5',['GN']]='AHSG'
	dff.loc[dff['GN']=='P02458',['GN']]='COL2A1'
	dff.loc[dff['GN']=='P08123',['GN']]='COL1A2'
	dff.loc[dff['GN']=='P02452',['GN']]='COL1A1'

	delete_string='REV_'
	dff=dff[dff['Accession'].str.contains(delete_string)==False]
	delete_string1='CON_'
	dff=dff[dff['Accession'].str.contains(delete_string1,na=False)==False]
	dff['Gene']=dff['GN']
	dff.to_csv(name+".protein-peptides.csv",index=False)
	print('--convert MaxQuant results to PEAKS format finished--')

def doNN(input_path,name):
	pd.set_option('display.unicode.east_asian_width',True)
	dff=pd.read_excel(input_path,usecols=['Protein.Group','Modified.Sequence','Q.Value','Stripped.Sequence','RT','Ms1.Area','Precursor.Normalised','Run','Ms2.Scan']) 
	dff=dff.dropna()
	dff.columns=['Source File','Peptide','Sequence','Accession','RT','Area','Intensity','-10LgP','Scan']
	dff['#Spec']='1'
	dff['Length']='0'
	dff['GN']='0'
	dff['PTM']='0'
	dff=dff.dropna(subset=['#Spec'])
	dff['#Spec']=dff['#Spec'].astype(int)
	dff['GN']=dff['Accession'].str.split('-').str[0] 
	dff['GN']=dff['GN'].str.split('_').str[0] 
	dff['GN']=dff['GN'].astype(str)
	dff.loc[dff['GN'].str.contains('sp',na=True),['GN']]=dff['Accession'].str.extract(r'|(\w+)|', expand=False)
	dff.loc[dff['GN'].str.contains('tr',na=True),['GN']]=dff['Accession'].str.extract(r'|(\w+)|', expand=False)
	dff.loc[dff['GN']=='CK820_G0031061',['GN']]='MMP20'
	dff.loc[dff['GN']=='CK820_G0035821',['GN']]='AMELX'
	dff.loc[dff['GN']=='CK820_G0053051',['GN']]='AMELY'
	dff=dff.dropna(subset=['GN'])
	dff['Length']=dff['Sequence'].str.len()
	dff['Gene']=dff['GN']
	dff['PTM']=dff['Peptide'].str.findall(r'\((.*?)\)')
	dff.to_csv(name+".protein-peptides.csv",index=False)
	print('--convert DIA-NN results to PEAKS format finished--')

def main():
	parser = argparse.ArgumentParser(prog='2PEAKSFormat.py',description='Build PEAKS output from pFind or MaxQuant or DIA-NN reuslt')
	parser.add_argument('input_path', help='The path of input file.')
	parser.add_argument('type', help='The type of input db file.')
	parser.add_argument('prefix', help='The prefix of output csv files.')
	args = parser.parse_args()
	input_path = args.input_path
	type = args.type
	name = args.prefix
	if type == "pFind":
		data = pd.read_csv(input_path, sep='\t',header=None,names=range(19),na_filter=False, skipinitialspace=True)
		data.to_excel(name+".PEAKS.xlsx",index=False)
		doPF(name+".PEAKS.xlsx", name)
	elif type == "MaxQuant":
		data = pd.read_csv(input_path, sep='\t')
		data.to_excel(name+".PEAKS.xlsx",index=False)
		doMQ(name+".PEAKS.xlsx", name)
	elif type == "DIA-NN":
		data = pd.read_parquet(input_path)
		data.to_excel(name+".PEAKS.xlsx", index=False)
		doNN(name+".PEAKS.xlsx", name)
	else:
		print("Please input correct type: pFind or MaxQuant or DIA-NN")

if __name__ == '__main__':
	main()
