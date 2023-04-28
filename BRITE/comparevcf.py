#!/usr/bin/env python
import sys
import pysam
import vcfpy

##TO RUN through terminal: make sure vcfpy is installed in environment, module load miniconda, conda activate vcfpy
## python3 comparevcf.py VCF1/PATH VCF2/PATH

vcf1 = sys.argv[1]
vcf2 = sys.argv[2]

#Read in VCF (use reader1 and reader2 for simplicity - but NOT actually reader objects)
#Consist of record object, accessible via index
reader1 = []
for record in vcfpy.Reader.from_path(vcf1): # will iterate through all of this one using basic for loop
    reader1.append(record)

reader2 = []
for record in vcfpy.Reader.from_path(vcf2): #better to read in longer file here
    reader2.append(record) 

# #Read in indist file
# indistinguishables=open('/projectnb/vntrseek/sfiller/hg38_rl150_iterable.indist','r')
# for i in range(10):
#     print(indistinguishables.readline())
# indistinguishables.close()


#FILTERS, set true to turn on
#Chrom, Pos, ID - always printed (can change later if not desired)
#if (bool) {print specific way}
write_to_file=False
#write_data_to_csv=False

print_all=False
print_identical=False
print_different = True

##Feature Selection
print_singletons=True #if True, add another layer of filtering before printing output
print_missing_TRIDs=False
print_missing_singletonTRIDs=False #True

if(print_all): #is true, then identical and different must both be printed
    print_identical=True
    print_different=True


def num_TRIDs(reader): #Count num of records per reader
    i = 0
    for record in reader:
        i+=1
    return str(i)

def num_reads(vcf_path): #pass either vcf1 or vcf2
    header = vcfpy.Reader.from_path(vcf_path).header.copy()
    numreads=header.lines[10].serialize().split('=')[1]
    return numreads

def num_readTRs(vcf_path): #pass either vcf1 or vcf2
    header = vcfpy.Reader.from_path(vcf_path).header.copy()
    numreadTRs=header.lines[11].serialize().split('=')[1]
    return numreadTRs

def num_VNTRs(vcf_path): #pass either vcf1 or vcf2
    header = vcfpy.Reader.from_path(vcf_path).header.copy()
    numVNTRs=header.lines[12].serialize().split('=')[1]
    return numVNTRs

def num_TRs_with_support(vcf_path): #pass either vcf1 or vcf2
    header = vcfpy.Reader.from_path(vcf_path).header.copy()
    numTRsWithSupport=header.lines[13].serialize().split('=')[1]
    return numTRsWithSupport

def get_int_ID(TRID): #pass record.ID to
    TRID = str(TRID)
    stripped = int(TRID.replace("td",""))
    return stripped

def prevent_nonetype(my_list): #makes sure that list is all int values)
    for i, element in enumerate(my_list):
        if (str(element)== "None"):
            my_list[i] = 0 
    return my_list
        
def print_match(record1,record2):
    line = [record1.CHROM, record1.POS, get_int_ID(record1.ID[0])]
    line += [record1.calls[0].data.get('GT'),record2.calls[0].data.get('GT'), record1.calls[0].data.get('SP'),record2.calls[0].data.get('SP'),record1.calls[0].data.get('CGL'),record2.calls[0].data.get('CGL'),int(GT), int(SP), int(CGL)]
    line='\t'.join(map(str,line))
    return line

def count_singletons(reader):
    singletons=0
    for record in reader:
        if(record.FILTER[0]=='PASS'):
            singletons+=1
    return str(singletons)

def print_and_write(output_string):
    print(output_string)
    if (write_to_file):
        new_file.write(output_string+'\n')
    return

# Build and print header & open file to write to
#header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'INFO', 'FILTER', 'FORMAT'] + reader.header.samples.names
vcf1_name = vcfpy.Reader.from_path(vcf1).header.samples.names[0]
vcf2_name = vcfpy.Reader.from_path(vcf2).header.samples.names[0]

if (write_to_file):
    file_name=str(vcf1_name+"_vs_"+vcf2_name+".tsv")
    #if (str(vcf1) and str(vcf2)) #if it contains 18545 in in name put it into this folder, else other genome
    #if it contains all with support put in all with support, else just regular folder
    if(str(vcf1).find('allwithsupport.span2')!= -1 and str(vcf2).find('allwithsupport.span2')!= -1 ):
        complete_path="/projectnb/vntrseek/sfiller/NYGC_na18545/OUTPUTcomparevcf/allTRswithsupport/"+file_name
    else:
        complete_path="/projectnb/vntrseek/sfiller/NYGC_na18545/OUTPUTcomparevcf/"+file_name
    new_file=open(complete_path,"w")

header=[]
files=["VCF1: " + vcf1_name] + ["VCF2: " +vcf2_name]
header+=["##"+'\t'.join(files)]
numTRIDs=[num_TRIDs(reader1),num_TRIDs(reader2)]
header+=["##numTRIDs="+'\t'.join(numTRIDs)]
#PRINTED STRAIGHT FROM VCF HEADER
numreads=[num_reads(vcf1),num_reads(vcf2)]
header+=["##numreads="+'\t'.join(numreads)]
numreadTRs=[num_readTRs(vcf1),num_readTRs(vcf2)]
header+=["##numreadTRs="+'\t'.join(numreadTRs)]
numVNTRs=[num_VNTRs(vcf1),num_VNTRs(vcf2)]
header+=["##numVNTRs="+'\t'.join(numVNTRs)]
numTRsWithSupport=[num_TRs_with_support(vcf1),num_TRs_with_support(vcf2)]
header+=["##numTRsWithSupport="+'\t'.join(numTRsWithSupport)] ## end of VCF header extracts
numsingletons=[count_singletons(reader1),count_singletons(reader2)]
header+=["##numsingletons="+'\t'.join(numsingletons)]
analysis_header = ['#CHROM', 'POS', 'ID', 'GT1', 'GT2','SP1','SP2','CGL1','CGL2'] + ['GT','SP','CGL']
header+=['\t'.join(analysis_header)]

for line in header:
    print_and_write(line)
records_singletons_sp_changes={}
num_same=0 #num with same TRID out of both files
identical=0 #num with exact same GT/SP/CGL for matching TRID
different=0 #num with at least one field of GT/SP/CGL different
num_missing_reader2=0 #num of TRIDs missing from reader2 that are found in reader1 (changes based on which file is which)
missing_reader2=[]
num_missing_reader1=0 #num of TRIDs missing from reader1 that are in reader2
missing_reader1=[]
gt_changes=0
sp_changes=0
cgl_changes=0
sp_unexpected=0
sp_expected=0
singleton_count=0
singletons_different=0
singletons_identical=0
num_singletons_missing_reader2=0 #out of total missing from reader 2, just the number that are singletons
singletons_missing_reader2=[]
num_singletons_missing_reader1=0 #out of total missing from reader 1, just the number that are singletons
singletons_missing_reader1=[]
singleton_gt_changes=0
singleton_sp_changes=0
singleton_cgl_changes=0
singleton_sp_unexpected=0
singleton_sp_expected=0

j=0
for record in reader1: #will iterate all the way through reader1, use index info to go through reader2
    id1 = get_int_ID(record.ID[0]) #only one ID in the list, so specify 0 index to remove list formatting, only returns string
    try: #can raise an exception if shorter than reader1
        id2 = get_int_ID(reader2[j].ID[0])
        while (id1 > id2):
            j+=1
            id2 = get_int_ID(reader2[j].ID[0])
    except: #for when there are more records in file we are iterating through, that need to be counted as missing from the other
        num_missing_reader2+=1
        missing_reader2+=[record.ID[0]]
        if (record.FILTER[0]=='PASS'):
            num_singletons_missing_reader2+=1
            singletons_missing_reader2+=[record]
    else:
        if (id1 == id2): 
            GT = False; SP = False; CGL = False #initially set to false, if different set to true (so that int value will display 1 if there is a difference)
            reader1SPlist = record.calls[0].data.get('SP')
            reader2SPlist = reader2[j].calls[0].data.get('SP')
            if (record.calls[0].data.get('GT') != reader2[j].calls[0].data.get('GT')):
                #print("different GT")
                GT = True
                gt_changes+=1
                if (record.FILTER[0]=='PASS'):
                    singleton_gt_changes+=1
            elif(record.calls[0].data.get('GT') == reader2[j].calls[0].data.get('GT')): #if GT is the same, check to see if SP is same, if not
                if (str(reader1SPlist) != str(reader2SPlist)):
                    sp_changes+=1
                    if (record.FILTER[0]=='PASS'):
                        singleton_sp_changes+=1
                        records_singletons_sp_changes[record.ID[0]]=reader1SPlist,reader2SPlist
                    reader1SPlist = prevent_nonetype(reader1SPlist) #converts 'None' to zero, ok because GT is known to be same, and list 
                    reader2SPlist = prevent_nonetype(reader2SPlist)
                    one_detected = False #only want to increase sp_unexpected once for each set of matching TRIDs (vs once for each time in list it changes)
                    for i in range(len(reader1SPlist)):
                        if ((reader1SPlist[i] > reader2SPlist[i]) and one_detected!=True):
                            sp_unexpected+=1 #num out of total sp_changes when support was higher in first record and lower in second record
                            #if reader1 is restricted and reader2 is full, num out of total sp_changes that are unexpected
                            #print("UNEXPECTED SP change:\t", reader1SPlist[i],reader2SPlist[i])
                            if (record.FILTER[0]=='PASS'):
                                singleton_sp_unexpected+=1
                            one_detected=True
                            
                    one_detected=False
                    for i in range(len(reader1SPlist)):
                        if ((reader1SPlist[i] < reader2SPlist[i]) and one_detected!=True):
                            sp_expected+=1 #num out of total sp_changes when support was lower in first record and higher in second record
                            #if reader1 is restricted and reader2 is full, num out of total sp_changes that are expected
                            #print("EXPECTED SP change:\t", reader1SPlist[i],reader2SPlist[i])
                            if (record.FILTER[0]=='PASS'):
                                singleton_sp_expected+=1
                            one_detected=True
                if (str(record.calls[0].data.get('CGL')) != str(reader2[j].calls[0].data.get('CGL'))):
                    cgl_changes+=1
                    if (record.FILTER[0]=='PASS'):
                        singleton_cgl_changes+=1

            if (str(reader1SPlist) != str(reader2SPlist)):
                #print("different SP")
                SP = True
            if (str(record.calls[0].data.get('CGL')) != str(reader2[j].calls[0].data.get('CGL'))):
                #print("different CGL")
                #add a feature to see if it's rare, if so add here
                CGL = True

            #COUNT SINGLETONS (total of each should add up to test)
            if(record.FILTER[0]=='PASS'and reader2[j].FILTER[0]=='PASS'): #don't need to use both, same result either way
                singleton_count+=1
                if(GT==False and SP==False and CGL==False):
                    singletons_identical+=1
                if(GT or SP or CGL):
                    singletons_different+=1

            if(GT==False and SP==False and CGL==False): #if all false print the line (there's no change in any of the fields)
                identical+=1
                if(print_identical):
                    if(print_singletons):
                        if(record.FILTER[0]=='PASS' and reader2[j].FILTER[0]=='PASS'):
                            print_and_write(print_match(record,reader2[j]))
                    else:
                        print_and_write(print_match(record,reader2[j]))

            if (GT or SP or CGL): #if one is true(different) print
                different+=1
                if(print_different):
                    if(print_singletons):
                        if(record.FILTER[0]=='PASS' and reader2[j].FILTER[0]=='PASS'):
                            print_and_write(print_match(record,reader2[j]))
                    else:
                        print_and_write(print_match(record,reader2[j]))
            j+=1
            num_same+=1
    finally:
        if(id2 > id1):
            num_missing_reader2+=1
            missing_reader2+=[record.ID[0]]
            if (record.FILTER[0]=='PASS'):
                num_singletons_missing_reader2+=1
                singletons_missing_reader2+=[record]

#When record 1 is shorter than 

#To count number of records missing from reader1
m=0
for record in reader2: #will iterate all the way through reader2, use index info to go through reader1
    id2 = get_int_ID(record.ID[0])
    try: #can raise an exception if shorter than reader2
        id1 = get_int_ID(reader1[m].ID[0])
        while (id2 > id1):
            m+=1
            id1 = get_int_ID(reader1[m].ID[0])
    except: #for when there are more records in file we are iterating through, that need to be counted as missing from the other
        num_missing_reader1+=1
        missing_reader1+=[record.ID[0]]
        if (record.FILTER[0]=='PASS'):
            num_singletons_missing_reader1+=1  #out of total missing from reader1
            singletons_missing_reader1+=[record]
    else:
        if (id1 == id2): 
            m+=1
    finally:
        if(id1 > id2): #(same output)
            num_missing_reader1+=1
            missing_reader1+=[record.ID[0]]
            if (record.FILTER[0]=='PASS'):
                num_singletons_missing_reader1+=1  #out of total missing from reader1
                singletons_missing_reader1+=[record]
 
summary=[]
summary+=[str(num_same)+"\t matching TRID entries"]
summary+=[str(identical)+"\t matching TRID and GT/SP/CGL"]
summary+=[str(different)+"\t matching TRID but at least one field of GT/SP/CGL is different"]
#if one has TRID and other is missing
summary+=[str(num_missing_reader2)+"\t number of entries missing from file2 that are in file1"]
summary+=[str(num_missing_reader1)+"\t number of entries missing from file1 that are in file2"]#
summary+=[str(gt_changes)+"\t num times genotype changed"]
summary+=["When genotype was the same:"]
summary+=[str(sp_changes)+"\t num times supporting reads changes when genotypes were the same"]
summary+=[str(cgl_changes)+"\t num times copies gained or lost changed when genotypes were the same"] #rare
#only unexpected when first record is restricted and second record is more full
summary+=[str(sp_unexpected)+"\t num times when support was higher in the first record  and lower in second record"]
summary+=[str(sp_expected)+"\t num times when support was lower in the first record and higher in second record"]
if (print_singletons):
    summary+=["The following numbers are out of all the numbers above; If only different singletons are printed, third line and following describe the output:"]
    summary+=["\t"+ str(singleton_count)+"\t num of singletons out of all matching TRIDs"]
    summary+=["\t"+ str(singletons_identical)+"\t num of singletons out of matching TRID and GT/SP/CGL"]
    summary+=["\t"+ str(singletons_different)+"\t num of singletons with matching TRID but at least one field of GT/SP/CGL is different"]
    summary+=["\t"+ str(num_singletons_missing_reader2)+"\t number of singletons out of all entries missing from file2 that are in file1"]
    summary+=["\t"+ str(num_singletons_missing_reader1)+"\t number of singletons out of all entries missing from file1 that are in file2"]#
    summary+=["\t"+ str(singleton_gt_changes)+"\t num times genotype changed within singletons"]
    summary+=["\t"+ "When genotype of singletons was the same:"]
    summary+=["\t"+ str(singleton_sp_changes)+"\t num times supporting reads changes within singletons when genotypes were the same"]
    summary+=["\t"+ str(singleton_cgl_changes)+"\t num times copies gained or lost changed within singletons when genotypes were the same"] #rare
    #only unexpected when first record is restricted and second record is more full
    summary+=["\t"+ str(singleton_sp_unexpected)+"\t num times when support was higher in the first record  and lower in second record within singletons"]
    summary+=["\t"+ str(singleton_sp_expected)+"\t num times when support was lower in the first record and higher in second record within singletons"]


for line in summary:
    print_and_write(line)

if (print_missing_TRIDs): #can print just numbers using print_and_write(get_int_ID(trid)) #returned as a string
    print_and_write("TRIDs (VNTRs) present in file1 that are missing from file2:")
    for trid in missing_reader2:
        print_and_write(str(trid))
    print_and_write("TRIDs (VNTRs) present in file2 that are missing from file1:")
    for trid in missing_reader1:
        print_and_write(str(trid))

if (print_missing_singletonTRIDs): #can print just numbers using print_and_write(get_int_ID(trid)) #returned as a string
    print_and_write("singleton TRIDs (VNTRs) present in file1 that are missing from file2:")
    for record in singletons_missing_reader2:
        print_and_write(str(record.ID[0])+ "\t"+ str(record.calls[0].data.get('SP')))  #wanted to look at support for missing reads
    print_and_write("singleton TRIDs (VNTRs) present in file2 that are missing from file1:")
    for record in singletons_missing_reader1:
        print_and_write(str(record.ID[0])+ "\t"+ str(record.calls[0].data.get('SP')))

if (write_to_file):
    new_file.close()


#alternative idea: load each record into a dictionary, with key being TRID, and value being CG/SP/CGL

# #Pretty repitive output, but good example of loading into dictionary
# print("TRIDS of singletons when genotype was the same but SP changed, shown are sp changes for each file")
# for key, value in records_singletons_sp_changes.items():
#     print(key, ' : ', value)