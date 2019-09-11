#structure as function / file operator
try:
    import pandas as pd
except ModuleNotFoundError:
    print('This script requires pandas. Please install pandas via conda in with the command conda install -c anaconda pandas')


import warnings
warnings.filterwarnings('ignore')


while True:
    print ('This script uses the eLife pCRISPRi-v2.1 file to give you the guides for genes you selected.')
    print('As a prerequisite, please save the eLife pCRISPRi-v2.1 file and the file that specifies the genes you want guides for in the same folder that this script is in, or provide their filepaths as the input')
    print('Please save your eLife pCRISPRi-v2.1 file as a comma-separated file (.csv)')
    elife_name=input('Please provide the name/filepath to your eLife file: ').strip()
    guides_name=input('Please provide the name of the file that contains your genes of interest: ').strip()
    number=int(input('Please provide the number of guides you want per gene: '))
    break

elife=pd.read_csv(elife_name)
# elife.drop(0, inplace=True)
# elife.dropna(how='all', axis=1, inplace=True)
# elife.dropna(how='all', inplace=True)
elife['selection rank']=elife['selection rank'].astype('int32')
elife=elife[elife['selection rank'].isin(list(range(1, number+1)))]



with open(guides_name) as f:
    guides=f.readlines()
guides=[x.strip('\n') for x in guides]
guides=[x.replace("'", "") for x in guides]



guides_filtered=elife[elife['gene'].isin(guides)]
# guides_filtered.to_csv('{}.csv'.format(guides_name))


def sequence_generator(sequence):
    FW_5='ttg'
    FW_3='gtttaagagc'
    RV_5='ttagctcttaaac'
    RV_3='caacaag'
    basedict={'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    bases=list(sequence)
    reversed_bases=reversed([basedict.get(base,base) for base in bases])
    reversed_bases=''.join(reversed_bases)
    complete_forward=FW_5+sequence+FW_3
    complete_reverse=RV_5+reversed_bases+RV_3
    return complete_forward, complete_reverse

ordering_frame=guides_filtered[['gene','transcript', 'selection rank', 'protospacer sequence']]

ordering_frame['oligo sequence FW'], ordering_frame['oligo sequence RV']=zip(*ordering_frame['protospacer sequence'].map(sequence_generator))
savename='{}_guides.csv'.format(guides_name.replace('.txt',''))
ordering_frame.to_csv(savename, index=False)

not_present=set(guides).difference(set(ordering_frame['gene']))

if len(not_present) != 0:
    print('Unable to find guides for: {}'.format(not_present))

print('Done. Files saved to ', savename)
