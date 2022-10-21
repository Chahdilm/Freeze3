from script.path_variable import *

import pandas as pd

import script.b_script_get_gene.get_gene_from_case as NT_p # import all info for phenopacket nodetype df

import script.b_script_get_gene.get_gene_from_orpha as NT_orpha # import all info for phenopacket nodetype df



####################################
## node type : ORPHA
df_pd1_filtered = NT_orpha.df_pd1_filtered
df_product = NT_orpha.df_product

# STEP one : put all orpha from product one into the df nodetype
df_pd1_filtered_dict = df_pd1_filtered.to_dict('index')
all_interractions = list()
for value in df_pd1_filtered_dict.values():
    oneORPHA = value['ORPHACode']
    onetype = value['DisorderGroup_1']
    all_interractions.append((oneORPHA,onetype))

# STEP two : chck on the df product which containt the orpha child from the product spe if all orpha child are on the df nodetype
#               if not we suppose that their type is  "Subtype of disorder"
df_product_dict = df_product.to_dict('index')
for value in df_product_dict.values():
    oneORPHA_c = value['ORPHAcode_child']
    if oneORPHA_c not in df_pd1_filtered['ORPHACode'].tolist():
        all_interractions.append((oneORPHA_c,"Subtype of disorder"))


df_type_orpha = pd.DataFrame(all_interractions, columns=['ORDO','type'])
df_type_orpha = df_type_orpha.drop_duplicates()

# verification if all orpha script child are on the df
# len(set(df_type_orpha['ORDO']).intersection(set(df_product['ORPHAcode_child'])))


####################################
############### node type gene id
df_pd6_filtered = NT_orpha.df_pd6_filtered

# STEP  : build df with all genes available on orphanet db and all id from other db related
list_pd6 = []
symbol =  source = reference = id_g = ''
with open(PATH_INPUT_PRODUCT_PD6, 'r') as f:
    ligne = f.readline()

    while ligne != "":
        if "<JDBOR" in ligne:
            print("PD6:\t", ligne)
        elif "<DisorderList" in ligne:
            print("PD6:\t", ligne)

        elif '<Gene id=' in ligne:
            # print('rentre symbol')
            new_ligne = ligne.replace('<Gene id="', '').strip()
            new_ligne = new_ligne.replace('">', '').strip()
            id_g = new_ligne
        elif '<Symbol>' in ligne:
            # print('rentre symbol')
            new_ligne = ligne.replace('<Symbol>', '').strip()
            new_ligne = new_ligne.replace('</Symbol>', '').strip()
            symbol = new_ligne
        elif '<ExternalReference id=' in ligne:
            ligne = f.readline()
            new_ligne = ligne.replace('<Source>', '').strip()
            new_ligne = new_ligne.replace('</Source>', '').strip()
            source = new_ligne

            ligne = f.readline()

            new_ligne = ""
            new_ligne = ligne.replace('<Reference>', '').strip()
            new_ligne = new_ligne.replace('</Reference>', '').strip()
            reference = new_ligne

        list_pd6.append((source,reference,symbol,id_g))
        # list_pd6.append(('symbol',symbol))
        ligne = f.readline()

df_type_gene = pd.DataFrame(list_pd6, columns=['source',"ref",'symbol','id'])
df_type_gene= df_type_gene.drop_duplicates()

#df_gene_filtered = df_gene[(df_gene["source"]  =="symbol") |(df_gene["source"]  =="Ensembl") | (df_gene["source"]  =="HGNC") | (df_gene["source"]  =="Reactome")].drop_duplicates()




####################################
## node type : solved/unsolved
all_pheno_list = NT_p.all_pheno_list_f # all json phenopacket
df_excels = NT_p.df_final_5_hpo_gene_afterC     # all solved case associater with one gene

all_cases = set()

for onepheno in all_pheno_list:
    if onepheno in df_excels['phenopacket'].tolist():
        # it s a solved one
        all_cases.add((onepheno,'SOLVED'))
    else:
        all_cases.add((onepheno,'UNSOLVED'))
df_type_p = pd.DataFrame(all_cases, columns=['case','type'])


####################################
## node ERN :

df_input_ERN = pd.read_csv(PATH_INPUT_NODE_TYPE_ERN ,sep=',')

all_interractions=set()
dict_ERN_df = df_input_ERN.to_dict('index')
for value in dict_ERN_df.values():
    oneERN=value['ERNs']
    oneERN_split = oneERN.split('-')
    if 'ERN-EURO-NMD' == oneERN:
        final_ERN= oneERN_split[2]
        all_interractions.add((value['PS ID'],final_ERN))
    elif oneERN_split[0] == 'Not available':
        pass
    else:
        final_ERN= oneERN_split[1]
        all_interractions.add((value['PS ID'],final_ERN))



df_ERN = pd.DataFrame(all_interractions, columns=['element','type'])

df_ERN_filtered= df_ERN[~df_ERN['type'].str.contains("Spain")]



####################################
print('export the df NODE type orpha')
df_type_orpha.to_csv(PATH_OUPUT_NODE_TYPE_ORPHA, sep='\t',index=False)

print('export the df NODE type genes')
df_type_gene.to_csv(PATH_OUPUT_NODE_TYPE_GENE, sep='\t',index=False)


print('export the df NODE type phenopacket')
df_type_p.to_csv(PATH_OUPUT_NODE_TYPE_PHENO, sep='\t',index=False)

print('export the df NODE type ERN')
df_ERN_filtered.to_csv(PATH_OUPUT_NODE_TYPE_ERN, sep='\t',index=False)

####################################