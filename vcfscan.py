import bisect
import re
import sys
from tqdm import tqdm


# the information for each mutation is stored in an instance of the Mutation class
class Mutation:

    def __init__(self, name, mut_type, mutation, sample_count, total_count):
        self.name = name
        self.mut_type = mut_type
        self.mutation = mutation
        self.sample_count = sample_count
        self.total_count = total_count
        self.gene_freq = 0
        self.rank = 0
        self.total = 0

# handling invalid input
if len(sys.argv) < 3:
    print("Usage: python vcfclone.py [cosmic_file1.vcf] [cosmic_file2.vcf] [output_file]")
    sys.exit()

# for efficiency, VCF files are treated as plain text and regex is used to extract
# useful data from each line
    
def re_capture(line):
    # defining capture patterns
    name_pattern =  re.compile(r'GENE=(.*?);')
    mut_pattern = re.compile(r'AA=(.*?);')
    sample_count_pattern = re.compile(r'SAMPLE_COUNT=(.*?);')
    can_stat_pattern = re.compile(r'IS_CANONICAL=(.*?);')
    # matching patterns in text
    name_match = name_pattern.search(line)
    mut_match = mut_pattern.search(line)
    sample_count_match = sample_count_pattern.search(line)
    can_stat_match = can_stat_pattern.search(line)
    # extracting data into variables
    name = name_match.group(1)
    mut = mut_match.group(1)
    sample_count = int(sample_count_match.group(1))
    can_stat = can_stat_match.group(1)

    # determining if a given gene was encountered and saving it in a dict if not
    # saving the number of counts for a mutation in a given gene
    target_index = bisect.bisect_left(gene_mutation_names, name)
    if target_index >= len(gene_mutation_names) or gene_mutation_names[target_index] != name:
        bisect.insort(gene_mutation_names, name)
        gene_mutation_counts[name] = [sample_count, []]
        rank_counts[name] = {sample_count: 1}
    else:
        gene_mutation_counts[name][0] += sample_count
        if sample_count in rank_counts[name]:
            rank_counts[name][sample_count] += 1
        else:
            rank_counts[name][sample_count] = 1
        
    
    # using MOAR regex to determine if a mutation is a loss of cysteine
    if re.search("p\.C", mut) != None and can_stat == "y":
        info = "UniProt info"
        if re.search("p\.\w*del\w*", mut) != None or re.search("p\.\w*ins\w*", mut) != None:
            mut_type = "indel"
        elif re.search("p\.\w*fs\w*", mut) != None:
            mut_type = "frameshift"
        elif re.search("p\.\w*dup\w*", mut) != None:
            mut_type = "duplication"
        else:
            mut_type = "point"
        new_mutation = Mutation(name = name, mut_type = mut_type, mutation = mut, sample_count = sample_count, total_count = 0)
        return new_mutation # if not a loss of cysteine, function implicitly returns None

# saving command-line input
filepath1 = sys.argv[1]
filepath2 = sys.argv[2]
output_filename = sys.argv[3]
output_file = open(output_filename, "w")

# defining file lengths (21 header lines) - use wc -l [filename] to find number of lines
FILE_LENGTH1 = (780292 - 21)
FILE_LENGTH2 = (39698736 - 21)

# defining data structures
mutations_list = []
gene_list = []
gene_mutation_names = []
gene_mutation_counts = {}
rank_counts = {}

print("Starting...")
# creating a progress bar
progress_bar = tqdm(total=FILE_LENGTH2, desc="  Parsing genomic mutation screen... ", unit=" mutations")
# opening file and feeding every line into above function
with open(filepath2, 'r') as mut_file2:
    for _ in range(21):
        header = mut_file2.readline()
    for line in mut_file2:
        new_mutation = re_capture(line)
        if new_mutation != None:
            mutations_list.append(new_mutation)
            gene_mutation_counts[new_mutation.name][1].append(new_mutation)
        progress_bar.update(1)

# closing progress bar
progress_bar.close()
sys.stdout.flush()
print("Scanning " + filepath2 + " done.")

progress_bar = tqdm(total=FILE_LENGTH1, desc="  Parsing genomic mutation screen... ", unit=" mutations")
with open(filepath1, 'r') as mut_file2:
    for _ in range(21):
        header = mut_file2.readline()
    for line in mut_file2:
        new_mutation = re_capture(line)
        if new_mutation != None:
            mutations_list.append(new_mutation)
            gene_mutation_counts[new_mutation.name][1].append(new_mutation)
        progress_bar.update(1)

progress_bar.close()
sys.stdout.flush()
print("Scanning " + filepath1 + " done.")

mutations_list_new = []

progress_bar = tqdm(total=len(mutations_list), desc="  Parsing genomic mutation screen... ", unit=" mutations")
# obtaining mutation frequencies
for gene in gene_mutation_names:
    if len(gene_mutation_counts[gene][1]) == 0:
        continue
    for mutation in gene_mutation_counts[gene][1]:
        search_pattern = mutation.mutation[:-1]
        for other_mutation in gene_mutation_counts[gene][1]:
            if gene_mutation_counts[gene][1].index(mutation) != gene_mutation_counts[gene][1].index(other_mutation) and search_pattern == other_mutation.mutation[:-1]:
                mutation.sample_count += other_mutation.sample_count
                gene_mutation_counts[gene][1].pop(gene_mutation_counts[gene][1].index(other_mutation))

        mutation.gene_freq = (mutation.sample_count / (gene_mutation_counts[gene][0])) * 100
        m_before = 0
        for key in rank_counts[mutation.name]:
            value = rank_counts[mutation.name][key]
            mutation.total += value
            if key > mutation.sample_count:
                m_before += value
            elif key == mutation.sample_count:
                if value % 2 == 0:
                    m_before += value / 2
                else:
                    m_before += (value - 1) / 2
        mutation.rank = m_before + 1
        mutations_list_new.append(mutation)
        progress_bar.update(1)

progress_bar.close()
# sorting mutations based on frequency
mutations_list_sorted = sorted(mutations_list_new, key = lambda x : x.gene_freq)

# identifying maximum length of each field for formatting purposes
max_name_len = max(len(mutation_entry.name) for mutation_entry in mutations_list_sorted)
max_type_len = max(len(mutation_entry.mut_type) for mutation_entry in mutations_list_sorted)
max_mutation_len = max(len(mutation_entry.mutation) for mutation_entry in mutations_list_sorted)

# formatting all the data for printing and writing to output file line by line
output_file.write("GENE\tMUTATION\tSAMPLE COUNT\tFREQUENCY\tRANK")
for i in range(len(mutations_list_sorted)):
    line = "{:<{}} {:<{}} {:<{}} {:<{}}\t {:<{}}%\t {:<{}}\n".format(str(mutations_list_sorted[i].name), max_name_len, mutations_list_sorted[i].mut_type, max_type_len, str(mutations_list_sorted[i].mutation), max_mutation_len, str(mutations_list_sorted[i].sample_count), 2, str(round(mutations_list_sorted[i].gene_freq, 3)), 6, str(int(mutations_list_sorted[i].rank)) + " / " + str(mutations_list_sorted[i].total), 8)
    output_file.write(line)

print("Output written to " + output_filename + ".")

# create a graphic explanation of how the program works (diagrams of data structures + examples)
# group cysteines into a single C position mutation (e.g, p.C49H p.C49R p.C49A >> p.C49x)
'''
 OUTPUT FORMAT
 GENE   MUTATION   COUNT   FREQ    RANK
 GAPDH  C53A       12      0.5%    1/250
'''