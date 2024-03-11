# -*- coding: utf-8 -*-

# Import packages
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
import warnings
warnings.filterwarnings("ignore")
# Declare positions of Spike gene in genome
start = 21562
end = 25383
# Open file with mutations' description and pack it to the dictionary
with open('mutations.txt', 'r') as document:
    mut = {}
    for line in document:
        line = line.split()
        if not line: continue
        mut[line[0]] = line[1:]
muts = list(mut.values())
# Open alignment and check amount of sequences
alignment = AlignIO.read("alignment.fa", "fasta")
num = len([1 for line in open("alignment.fa") if line.startswith(">")])
# Open reference sequence
record = SeqIO.read("reference.fa", "fasta")
s = list(record.seq)
rows = []
# Make table with INDELs
for j in range(1, num):
    sequence = list(alignment[j])
    n = 0
    for z in range(start, end):
        if 'n' in sequence[z]: n += 1
    nCount = "{:.0%}".format(1 - n / (end - start))
    conclusion = "0"
    rows_list = []
    for i in range(0, len(alignment[0])):
        if alignment[j][i] != alignment[0][i]:
            dict1 = {'REF': alignment[0][i], 'POS': i, 'ALT': alignment[j][i]}
            rows_list.append(dict1)
    df0 = pd.DataFrame(rows_list)
    df0.ALT = df0.ALT.replace('-', 'del')
    df0.REF = df0.REF.replace('-', 'ins')
    k = 0
    ins = 0
    # Search mutations and change positions if ins or del
    for ii in range(0, df0.ALT.count()):
        if df0.ALT.iloc[ii] == "del": k += 1
        if df0.REF.iloc[ii] == "ins":
            df0.POS.iloc[ii] = df0.POS.iloc[ii] - k + 2
            ins += 1
        if (df0.ALT.iloc[ii] != "del"): df0.POS.iloc[ii] = df0.POS.iloc[ii] - ins + 1
    # Leave nucleotide deletions for searching amino acids deletions
    dd1 = df0[df0.ALT == 'del']
    # Search deletions in Spike protein
    dd2 = dd1[(dd1.POS >= start) & (dd1.POS <= end)]
    delt = []
    delt = (dd2.POS - start) / 3
    delt = [int(x) for x in delt]
    delt = list(set(delt))
    delt = sorted(delt)  # позиция выпавшей аминокислоты в геноме
    df0 = df0.drop(df0[df0.ALT == 'del'].index)
    df0 = df0.drop(df0[df0.REF == 'ins'].index)
    df0 = df0.dropna()
    # Open pseudo-reference to compare with master reference
    probe = list(alignment[0].seq)
    for i in range(0, len(df0)):
        if (df0.POS.iloc[i] - 1) < len(probe): probe[df0.POS.iloc[i] - 1] = df0.ALT.iloc[i]  # substitute mutations in "reference"
    spike1 = s[start:end]
    spike1 = ''.join(spike1)  # translate the Spike gene sequence into a string form for the reference
    spike2 = probe[start:end]
    spike2 = ''.join(spike2)  # translate the Spike gene sequence into a string form for the pseudo-reference
    translation = Seq(spike1).translate()
    delet = []
    # Make array of amino acid's deletions
    for n in delt:
        if len(translation) > n:
            delet.append(str(str(translation[n]) + str(n + 1) + 'del'))
        else:
            continue
    # Function for creating a table of translated sequences with mutatations
    def trans(seq, probe):
        trans1 = Seq(seq).translate()
        trans2 = Seq(probe).translate()
        data_tuples = list(zip(trans1, trans2))
        dff = pd.DataFrame(data_tuples, columns=['Ref', 'Probe'])
        arr = []
        new_df = dff[dff['Ref'] != dff['Probe']]
        new_df['index'] = new_df.index
        new_df['index'].loc[::] += 1
        new_df = new_df.drop(new_df[new_df.Probe == 'X'].index)
        df = new_df[['Ref', 'index', 'Probe']]
        if df.empty:
            reps = ''
        else:
            x = df.to_string(header=False, index=False, index_names=False).split('\n')
            reps = [''.join(ele.split()) for ele in x]
            arr.append(reps)
        return arr
    # Get string with all mutations in our sample
    if not delet:
        sp = str(trans(spike1, spike2))
    else:
        sp = str(delet) + ', ' + str(trans(spike1, spike2))
    # Delete unnecessary characters in the list of mutations
    chars = ["[", "]", "'"]
    for c in chars: sp = sp.replace(c, "")
    # Make table: Sample-Mutations
    dict0 = {'Sample': [1], 'Conclusion': [2], 'Spike_coverage': [3], 'Spike mutations': [4]}
    dd1 = pd.DataFrame(columns=dict0.keys())
    dd2 = pd.DataFrame(columns=mut.keys())
    result = pd.concat([dd1, dd2])
    # Fill the table
    result['Sample'] = [alignment[j].name]
    result['Conclusion'] = [conclusion]
    result['Spike_coverage'] = [nCount]
    result['Spike mutations'] = [sp]
    # Declare list of unique mutations for genetic lines
    list_omXBB15 = ['F486P']
    list_omXBB = ['Y144del','V83A', 'H146Q','Q183E', 'V213E', 'G252V', 'L368I', 'V445P','F490S']
    list_om_2_75_2 = ['D1199N']
    list_om_2_75 = ['K147E', 'W152R', 'F157L', 'I210V', 'G257S']
    list_omBQ1 = ['K444T']
    list_omBF7 = ['R346T']
    list_om_5 = ['F486V']
    list_om_2 = ['X666X'] # no unique mutations
    kolvo = []
    # Check presense of mutatinos for each line
    for k in range(0, (len(result.columns) - 4)):
        my_list = result['Spike mutations'].loc[0].split(", ")
        same_values = set(my_list) & set(muts[k])
        # Conclude on the number of mutations found for each line
        if len(same_values) == 0:
            result.iloc[0][k + 4] = '-'
        elif len(same_values) == len(muts[k]):
            result.iloc[0][k + 4] = "100%"
        else:
            result.iloc[0][k + 4] = str(len(same_values)) + "/" + str(len(muts[k])) + " " + "{:.0%}".format(
                (len(same_values) / len(muts[k])))
        # Make a conclusion about the line based on the unique mutations of the lines found
        if any(t in x for x in my_list for t in list_omXBB15):
            result['Conclusion'].iloc[0] = "Omicron(XBB.1.5)"
            
        elif any(t in x for x in my_list for t in list_omXBB):
            result['Conclusion'].iloc[0] = "Omicron(XBB)"  
            
        elif any(t in x for x in my_list for t in list_om_2_75_2):
            result['Conclusion'].iloc[0] = "Omicron(BA.2.75.2)"
            
        elif any(t in x for x in my_list for t in list_om_2_75):
            result['Conclusion'].iloc[0] = "Omicron(BA.2.75)"
            
        elif any(t in x for x in my_list for t in list_omBQ1):
                result['Conclusion'].iloc[0] = "Omicron(BQ.1)" 
                
        elif any(t in x for x in my_list for t in list_omBF7):
            result['Conclusion'].iloc[0] = "Omicron(BF.7)"
            
        elif any(t in x for x in my_list for t in list_om_5):
            result['Conclusion'].iloc[0] = "Omicron(BA.4/BA.5)"             
            
        elif any(t in x for x in my_list for t in list_om_2):
            result['Conclusion'].iloc[0] = "Omicron(BA.2)"            
                
        else:
            result['Conclusion'] = "undefined"
        # Write information in a table
    rows.append(result)
    table = pd.concat(rows, ignore_index=True)
    table.reset_index(drop=True, inplace=True)
# Make resulting table and save in a file
from datetime import datetime
date = datetime.now().strftime('%Y-%m-%d_%H-%-M')
filename = f'Result_{date}.csv'
table.to_csv(filename, sep="\t", index=False)
