import sys
from Bio import SeqIO
from Bio.Seq import Seq
import math
class Individual:
    def __init__(self, sample_id, date, phenotype, sequence):
        self.sample_id = sample_id
        self.date = date
        self.phenotype = phenotype
        self.sequence = sequence

def read_fasta(file_path):
    individuals = []
    for record in SeqIO.parse(file_path, "fasta"):
        header_parts = record.description.split()
        sample_id, date = header_parts[0].split('_')
        sequence = str(record.seq)
        protein_sequence = Seq(sequence).translate()
        phenotype = "blue" if protein_sequence[3] == "S" else "orange"
        individuals.append(Individual(sample_id, date, phenotype, sequence))
    return individuals

def calculate_probability(k, n, p):
    probability = (math.factorial(n) / (math.factorial(n - k) * math.factorial(k))) * (p ** k) * ((1 - p) ** (n - k))
    return probability

def main():
    if len(sys.argv) != 4:
        print("Usage: python popgen.py Aubie.fasta 0.3 results.txt")
        sys.exit(1)
    file_path = sys.argv[1]
    frequency = float(sys.argv[2])
    output_file = sys.argv[3]
    individuals = read_fasta(file_path)
    n = len(individuals)
    k = sum(1 for ind in individuals if ind.phenotype == "orange")
    probability = calculate_probability(k, n, frequency)

    with open(output_file, 'w') as result_file:
        result_file.write("Results" + "\n\n")
        result_file.write('p (the frequency of "orange" in the population) = '+str(frequency)+ "\n")
        result_file.write(f"n (the number of sampled individuals) = {n}\n")
        result_file.write(f'k (the number of "orange" individuals in the sample set) = {k}\n\n')
        result_file.write(f'Probability of collecting 32 individuals with 5 being "orange" (given a population frequency of 0.3) = {probability}\n')

if __name__ == "__main__":
    main()
