#!usr/bin/python

import sys
import json
import vcf
#from collections import Counter
import re


sample_input = sys.argv
sample_input_1 = sample_input[1]
print(sample_input_1)

vcf_reader = vcf.Reader(open('hg96-hg97.chr21.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf', 'r'))
##counting the number of variants and indels
snp1 = ""
snp2 = ""
snp3 = ""
snp4 = ""
snp5 = ""
count = 0
x = []
y = []
y1 = []
count_1 = 0
count_2 = 0
count_3 = 0
AFR_AF = []
ASN_AF = []
AMR_AF = []
EUR_AF = []
count_4 = 0
count_5 = 0
count_6 = 0
count_7 = 0
origin = []

## record grabs all the column and information within it form the vcf file
sample_name = vcf_reader.samples
for record in vcf_reader:
    ###parsing for certain sample Id 
    call = record.genotype(sample_input_1)
    call_site = call.site
    ###grabbing information for certain sample only form the whole list 
    if call.sample == sample_input_1:
        
        ###############Bonus Allele Frequency##########################################
        ##Allele Freq value for African
        AF1 = record.INFO.get('AFR_AF')
        if AF1 != None:
            African_a = AF1
        #Allele Freq value for Asian 
        AF2 = record.INFO.get('ASN_AF')
        if AF2 != None:
            Asian_a = AF2
        #allele frequency value for American
        AF3 = record.INFO.get('AMR_AF')
        if AF3 != None:
            American_a = AF3
        #Allele Frequency value for European
        AF4 = record.INFO.get('EUR_AF')
        if AF4 != None:
            European_a = AF4
        
        #Getting the highest Allele Frequency value for each record
        high = max(African_a, American_a, Asian_a, European_a)
        #print ("max value", high)
        
        ##couting the total number of snps based on the highest allele freq values 
        if high == African_a:
            #print ("snp african")
            count_4 = count_4 + 1
        if high == Asian_a:
            #print ("snp asian")
            count_5 = count_5 + 1
        if high == American_a:
            #print ("snp american")
            count_6 = count_6 + 1
        if high == European_a:
            #print ("snp European")
            count_7 = count_7 + 1
        
        ##to get the variants for a sample ID
        alternate = record.var_type 
        x.append(alternate)
        
        ##counting indels
        if record.is_indel == True:
            count = count + 1
        
        #counting total snps
        if record.is_snp == True:
         count_1 = count_1 + 1
        
        #counting transition
        if record.is_transition == True:
            count_2 = count_2 + 1
    
        ###couting for each snps in a sample
        if record.var_type == 'snp':
             ##couting snp for A
            if call_site.REF == 'A':
                snp_A = call_site.ALT
                snp1 += str(snp_A)
                ##couting for C
            if call_site.REF == 'C':
                snp_C = call_site.ALT
                snp2 += str(snp_C)
                ##counting for T
            if record.REF == 'T':
                snp_T = call_site.ALT
                snp3 += str(snp_T)
                ##counting for G
            if call_site.REF == 'G':
                snp_G = record.ALT
                snp4 += str(snp_G)
                #counting for N
            if call_site.REF == 'N':
                snp_N = call_site.ALT
                snp5 += str(snp_N)

#########################finding origins#################################################

###finding the highest number of snps for finding the origin
high1 = max(count_4, count_5, count_6, count_7)

if high1 == count_4:
    print ("Origins African since highest humber of snps found for African", count_4)
    origin = 'African'

if high1 == count_5:
    print ("Origins Asian since highest humber of snps found for Asian", count_5)
    origin = 'Asian'

if high1 == count_6:
    print ("Origins American since highest humber of snps found for American", count_6)
    origin = 'American'

if high1 == count_7:
    print ("Origins European since highest humber of snps found for European", count_7)
    origin = 'European'

#gives sample name 
#print (call.sample)
#total number of variants in the file 
#print("count of variants",len(x))

###snp count for A
#print ("snp count for A")
count_A =  snp1.count('A')
count_G =  snp1.count('G')
count_C =  snp1.count('C')
count_T =  snp1.count('T')
count_N =  snp1.count('N')

#print ("G",count_G)
#print ("C", count_C)
#print ("T", count_T)
#print ("N", count_N)
#print ("A", count_A)

##snp count for C
#print ("snp count for C")

count1_A =  snp2.count('A')
count1_G =  snp2.count('G')
count1_C =  snp2.count('C')
count1_T =  snp2.count('T')
count1_N =  snp2.count('N')

#print ("G", count1_G)
#print ("C", count1_C)
#print ("T", count1_T)
#print ("N", count1_N)
#print ("A", count1_A)

##snp count for T
#print ("snp count for T")
count2_A =  snp3.count('A')
count2_G =  snp3.count('G')
count2_C =  snp3.count('C')
count2_T =  snp3.count('T')
count2_N =  snp3.count('N')

#print ("G", count2_G)
#print ("C", count2_C)
#print ("T", count2_T)
#print ("N", count2_N)
#print ("A", count2_A)

##snp count for G
#print ("snp count for G")
count3_A =  snp4.count('A')
count3_G =  snp4.count('G')
count3_C =  snp4.count('C')
count3_T =  snp4.count('T')
count3_N =  snp4.count('N')

#print ("G", count3_G)
#print ("C", count3_C)
#print ("T", count3_T)
#print ("N", count3_N)
#print ("A", count3_A)

##snp count for N
#print ("snp count for N")
count4_A =  snp5.count('A')
count4_G =  snp5.count('G')
count4_C =  snp5.count('C')
count4_T =  snp5.count('T')
count4_N =  snp5.count('N')

#print ("G", count4_G)
#print ("C", count4_C)
#print ("T", count4_T)
#print ("N", count4_N)
#print ("A", count4_A)

##counting transversions; A->C, T->G, A->T & C->G; are classified for transversions 

transversions = ((count_C) + (count2_G) + (count_T) + (count1_G))
#print (transversions)
###couting Ti/Tv rates
rate = (count_2/transversions)
#print ("Ti/Tv rates",rate)
        
#print("all variant types",x)
#print("count of variants",len(x))
#print (x[0])
sv_set = set(x)
sv_count = len(sv_set) #total number of SVs on file
##print ("Structural Variants",sv_count)
#counts for number of indel
#print ("number of indel",count)
#prints for snps total
#print ("number of snps",count_1)
#count for transitions 
#print ("count for Transitions",count_2)
#counting snp

###########################################JSON OUTPUT################################################

data_A = {
    "snps": {
      "A": {
      "A": count_A,
      "C": count_C,
      "T": count_T,
      "G": count_G,
      "N": count_N
    },
    "C": {
      "A": count1_A,
      "C": count1_C,
      "T": count1_T,
      "G": count1_G,
      "N": count1_N
    }
        }
}
data_B = {
    "T": {
      "A": count2_A,
      "C": count2_C,
      "T": count2_T,
      "G": count2_G,
      "N": count2_N
    },
    "G": {
      "A": count3_A,
      "C": count3_C,
      "T": count3_T,
      "G": count3_G,
      "N": count3_N
    }
}
data_C = {
    "N": {
    "A": count4_A,
    "C": count4_C,
    "T": count4_T,
    "G": count4_G,
    "N": count4_N
    }
 },{
    "sample": call.sample,
    "ti-tv": rate,
    "variant_count": len(x),
    "indel_count": count,
    "sv_count": sv_count      
}
data_D = {
    "Origin for the individual is ": origin,
    "since the snps have the highest count of": high1
}
json_str = json.dumps(data_A)
json_str1 = json.dumps(data_B)
json_str2 = json.dumps(data_C)
file_name = ''.join([str(sample_input_1),'.json'])
print (file_name)
# Writing JSON data
with open(file_name, 'w') as f:
     json.dump(data_A, f, indent=2)
     json.dump(data_B,f, indent=2)
     json.dump(data_C,f, indent=2)
     json.dump(data_D,f, indent=2)
print ("done!")

