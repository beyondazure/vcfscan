import re
import sys
import re
from tqdm import tqdm
import bisect


class Mutation:

    def __init__(self, name, mut_type, mutation, sample_count, total_count):
        self.name = name
        self.mut_type = mut_type
        self.mutation = mutation
        self.sample_count = sample_count
        self.total_count = total_count
        self.gene_freq = 0

if len(sys.argv) < 3:
    print("Usage: python vcfclone.py <filepath 1> <filepath 2> <output file>")
    sys.exit()

def re_capture(line):
    name_pattern =  re.compile(r'GENE=(.*?);')
    mut_pattern = re.compile(r'AA=(.*?);')
    sample_count_pattern = re.compile(r'SAMPLE_COUNT=(.*?);')
    can_stat_pattern = re.compile(r'IS_CANONICAL=(.*?);')
    name_match = name_pattern.search(line)
    mut_match = mut_pattern.search(line)
    sample_count_match = sample_count_pattern.search(line)
    can_stat_match = can_stat_pattern.search(line)
    name = name_match.group(1)
    mut = mut_match.group(1)
    sample_count = int(sample_count_match.group(1))
    can_stat = can_stat_match.group(1)
    target_index = bisect.bisect_left(gene_mutation_names, name)
    if target_index >= len(gene_mutation_names) or gene_mutation_names[target_index] != name:
        bisect.insort(gene_mutation_names, name)
        gene_mutation_counts[name] = sample_count
    else:
        gene_mutation_counts[name] += sample_count
    

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
        return new_mutation

filepath1 = sys.argv[1]
filepath2 = sys.argv[2]
output_filename = sys.argv[3]
output_file = open(output_filename, "w")

FILE_LENGTH1 = (780292 - 21)
FILE_LENGTH2 = (39698736 - 21)

mutations_list = []
id_list = []
gene_mutation_names = []
gene_mutation_counts = {}

print("Starting...")
progress_bar = tqdm(total=FILE_LENGTH2, desc="  Parsing genomic mutation screen... ", unit=" mutations")
with open(filepath2, 'r') as mut_file2:
    for _ in range(21):
        header = mut_file2.readline()
    for line in mut_file2:
        new_mutation = re_capture(line)
        if new_mutation != None:
            mutations_list.append(new_mutation)
        progress_bar.update(1)

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
        progress_bar.update(1)

progress_bar.close()
sys.stdout.flush()
print("Scanning " + filepath1 + " done.")

for mutation in mutations_list:
    mutation.gene_freq = (mutation.sample_count / (gene_mutation_counts[mutation.name])) * 100

mutations_list_sorted = sorted(mutations_list, key = lambda x : x.gene_freq)

max_name_len = max(len(mutation_entry.name) for mutation_entry in mutations_list_sorted)
max_type_len = max(len(mutation_entry.mut_type) for mutation_entry in mutations_list_sorted)
max_mutation_len = max(len(mutation_entry.mutation) for mutation_entry in mutations_list_sorted)

for i in range(len(mutations_list_sorted)):
    line = "{:<{}} {:<{}} {:<{}} {:<{}} {:<{}}%\n".format(str(mutations_list_sorted[i].name), max_name_len, mutations_list_sorted[i].mut_type, max_type_len, str(mutations_list_sorted[i].mutation), max_mutation_len, str(mutations_list_sorted[i].sample_count), 2, str(round(mutations_list_sorted[i].gene_freq, 3)), 6)
    output_file.write(line)

print("Output written to " + output_filename + ".")

    