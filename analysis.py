# -*- coding: utf-8 -*-

# Импортируем требуемые для анализа библиотеки из предустановленного пакета biopython
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
import warnings
warnings.filterwarnings("ignore")
# Объявляем позиции гена в геноме
start = 21562
end = 25383
# Открываем файл с мутациями и переводим его в вид словаря
with open('mutations.txt', 'r') as document:
    mut = {}
    for line in document:
        line = line.split()
        if not line: continue
        mut[line[0]] = line[1:]
muts = list(mut.values())
# Открываем выравненную последовательность и вычисляем количество сиквенсов
alignment = AlignIO.read("alignment.fa", "fasta")
num = len([1 for line in open("alignment.fa") if line.startswith(">")])
# Открываем референсную последовательность
record = SeqIO.read("reference.fa", "fasta")
s = list(record.seq)
rows = []
# Делаем таблицу сиквенса референса и конкретной пробы, обозначаем делеции и вставки
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
    # Ищем мутации, сдвигая нумерацию позиций, если они найдены
    for ii in range(0, df0.ALT.count()):
        if df0.ALT.iloc[ii] == "del": k += 1
        if df0.REF.iloc[ii] == "ins":
            df0.POS.iloc[ii] = df0.POS.iloc[ii] - k + 2
            ins += 1
        if (df0.ALT.iloc[ii] != "del"): df0.POS.iloc[ii] = df0.POS.iloc[ii] - ins + 1
    # Оставляем делеции для поиска выпавших аминокислот
    dd1 = df0[df0.ALT == 'del']
    # Ищем делеции в Spike белке (позиции в геноме 21563:25384)
    dd2 = dd1[(dd1.POS >= start) & (dd1.POS <= end)]
    delt = []
    delt = (dd2.POS - start) / 3
    delt = [int(x) for x in delt]
    delt = list(set(delt))
    delt = sorted(delt)  # позиция выпавшей аминокислоты в геноме
    df0 = df0.drop(df0[df0.ALT == 'del'].index)
    df0 = df0.drop(df0[df0.REF == 'ins'].index)
    df0 = df0.dropna()
    # Открываем "подставной" референс, который будет равен мутировавшему штамму и сравнен с неизмененным референсом
    probe = list(alignment[0].seq)
    for i in range(0, len(df0)):
        if (df0.POS.iloc[i] - 1) < len(probe): probe[df0.POS.iloc[i] - 1] = df0.ALT.iloc[
            i]  # подставляем мутации в "референс"
    spike1 = s[start:end]
    spike1 = ''.join(spike1)  # переводим последовательность Spike гена в строчный вид для настоящего референса
    spike2 = probe[start:end]
    spike2 = ''.join(spike2)  # переводим последовательность Spike гена в строчный вид для подставного референса
    translation = Seq(spike1).translate()
    delet = []
    # Складываем в массив делеции в аминокислотах
    for n in delt:
        if len(translation) > n:
            delet.append(str(str(translation[n]) + str(n + 1) + 'del'))
        else:
            continue
    # Функция, создающая таблицу транслированных последовательностей пробы и референса, и выводящая замены
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
    # Выводим строку со всеми мутациями в анализируемой пробе
    if not delet:
        sp = str(trans(spike1, spike2))
    else:
        sp = str(delet) + ', ' + str(trans(spike1, spike2))
    # Удаляем лишние знаки в списке мутаций
    chars = ["[", "]", "'"]
    for c in chars: sp = sp.replace(c, "")
    # Создаем таблицу: Проба-Мутации
    dict0 = {'Sample': [1], 'Conclusion': [2], 'Spike_coverage': [3], 'Spike mutations': [4]}
    dd1 = pd.DataFrame(columns=dict0.keys())
    dd2 = pd.DataFrame(columns=mut.keys())
    result = pd.concat([dd1, dd2])
    # Заполняем таблицу данными
    result['Sample'] = [alignment[j].name]
    result['Conclusion'] = [conclusion]
    result['Spike_coverage'] = [nCount]
    result['Spike mutations'] = [sp]
    # Объявляем список уникальных мутаций для заданных линий
    list_omXBB15 = ['F486P']
    list_omXBB = ['Y144del','V83A', 'H146Q','Q183E', 'V213E', 'G252V', 'L368I', 'V445P','F490S']
    list_om_2_75_2 = ['D1199N']
    list_om_2_75 = ['K147E', 'W152R', 'F157L', 'I210V', 'G257S']
    list_omBQ1 = ['K444T']
    list_omBF7 = ['R346T']
    list_om_5 = ['F486V']
    list_om_2 = ['X666X']
    kolvo = []
    # Проверяем, есть ли мутации для каждой разновидности вируса
    for k in range(0, (len(result.columns) - 4)):
        my_list = result['Spike mutations'].loc[0].split(", ")
        same_values = set(my_list) & set(muts[k])
        # Делаем заключение о количестве найденных мутаций для каждой линии
        if len(same_values) == 0:
            result.iloc[0][k + 4] = '-'
        elif len(same_values) == len(muts[k]):
            result.iloc[0][k + 4] = "100%"
        else:
            result.iloc[0][k + 4] = str(len(same_values)) + "/" + str(len(muts[k])) + " " + "{:.0%}".format(
                (len(same_values) / len(muts[k])))
        # Делаем заключение о линии по найденным уникальным мутациям линий
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
        # Записываем все в таблицу
    rows.append(result)
    table = pd.concat(rows, ignore_index=True)
    table.reset_index(drop=True, inplace=True)
# Записываем всё в итоговую таблицу и сохраняем на диск
from datetime import datetime
date = datetime.now().strftime('%Y-%m-%d_%H-%-M')
filename = f'Result_{date}.csv'
table.to_csv(filename, index=False)
