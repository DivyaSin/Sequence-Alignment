from Bio import SeqIO
from random import (choice, random, randint)

class Chromosome:
    target_sequence = "GTAAATTCATACCGGAAATTTTACCAAATGGCGATTTCTTAATTGCTGAGGTGGCCAGCAGAAATCGTCTTTTCATTATTCTGGAATCAAAACACATTCTTTGAATTGTTCACTTTTCTGTTGCCTTGAAATCTTGGTCTTCTTAGTTGACTGTTTCATCAAGGTTGCTCCAAATTCTTTGTGATTTATTGGTAAACTCGGGCATTTTATTGAGT"
    print "Length of target_sequence"
    print len(target_sequence)
    print "Target sequence:", target_sequence

    def __init__(self, gene):
        # Initialize the gene with its fitness
        self.gene = gene
        self.fitnessValue = Chromosome.get_fitness(gene)
    
    def mate(self, mate):
        # Returns new chromosomes after mating
        i = randint(0, len(self.gene) - 1)
        gene1 = self.gene[:i] + mate.gene[i:]
        gene2 = mate.gene[:i] + self.gene[i:]    
        return Chromosome(gene1), Chromosome(gene2)
    
    def mutate(self):
        # Mutate the gene by changing a random character 
        gene = list(self.gene)
        small_change = randint(32, 121)
        i = randint(0, len(gene) - 1)
        gene[i] = chr((ord(gene[i]) + small_change) % 122)
        # where chr returns string of one character on the basis of given ASCII code
        
        return Chromosome(''.join(gene))

    @staticmethod            
    def get_fitness(gene):
        # Returns the fitness of each chromosome
        fitnessValue = 0
        for a, b in zip(gene, Chromosome.target_sequence):
            fitnessValue += abs(ord(a) - ord(b))         
        # print "gene: ", gene
        # print "Chromosome.target_sequence: ", Chromosome.target_sequence
        # print fitness_value  
        return fitnessValue
        
class Population:    
    selectionSize = 3
    def __init__(self, crossover=0.7, mutation=0.3):
        print "Initializing population ...."
        self.mutation = mutation
        self.crossover = crossover
        fasta_sequences = SeqIO.parse(open("Ca.fa"), 'fasta')
        with open("output_file.txt") as out_file:
            sequences = []
            for fasta in fasta_sequences:
                name, sequence = fasta.id, str(fasta.seq) 
                # print sequence
                # write_fasta(out_file)
            # 215 is length of target sequence
            n = 215
            for i in range(0, len(sequence), n):
                sequences.append(Chromosome(sequence[i:i+n]))
            # sequences = [Chromosome(sequence[i:i+n]) for i in range(0, len(sequence), n)]
            self.population = list(sorted(sequences, key=lambda k: k.fitnessValue))

    # def insert_gaps(lmax, len_seq_record, intial_population_size):
    #     length = 1.2*lmax
    #     num_of_gaps = length - len_seq_record
    #     pos1 = randint(0, len_input_sequence1 - 1)  # pick random position to insert gap
    #     query_sequence = "".join((input_sequence1[:pos], "-", input_sequence1[pos:]))  # insert char at gap
    #     pos2 = randint(0, len_seq_record - 1)
    #     target_sequence = "".join((seq_record[:pos], "-", seq_record[pos:]))
    #     return query_sequence, target_sequence

    # query_sequence, target_sequence = insert_gaps(lmax, len_seq_record, intial_population_size)

    def select_chromosome(self):
        # Select a random chromosome using tournament selection 
        b = choice(self.population)
        for i in range(Population.selectionSize):
            c = choice(self.population)
            if (c.fitnessValue < b.fitnessValue):
                b = c                    
        return b

    def select_parents(self):      
        return (self.select_chromosome(), self.select_chromosome())
        
    def evolve(self):
        # print "Generating new population ...."
        # Generate new population using selection, mating and mutation operators
        size = len(self.population)
        # print size
        i = int(round(size * 0.1))
        sequences = self.population[:i]    
        while (i < size):
            # print "Crossover ....."
            if random() <= self.crossover:
                # print "Selecting parents ....."
                (parent1, parent2) = self.select_parents()
                # print "Mating ....."
                children = parent1.mate(parent2)
                # Mutate the children
                for child in children:
                    # print "Mutate and add back to population"
                    if random() <= self.mutation:
                        sequences.append(child.mutate())
                    else:
                        sequences.append(child)
                i = i+2
            else:
                # print "Mutate and add back to population"
                if random() <= self.mutation:
                    sequences.append(self.population[i].mutate())
                else:
                    sequences.append(self.population[i])
                i = i+1    
        self.population = list(sorted(sequences[:size], key=lambda k: k.fitnessValue))

# Main function
if __name__ == "__main__":
    maximumGenerations = 100
    print
    print "Maximum no. of Generations:", maximumGenerations
    print
    print "Fitness value will be decreasing, when it is zero then sequence aligns with our target sequence."
    print
    p = Population(crossover=0.6, mutation=0.3)
    minimumFitness = p.population[0].fitnessValue
    for i in range(1, maximumGenerations + 1):
        print
        print "-------------------------------------------------------------------------------------------"
        print("Generation %d: %s, %d" % (i, p.population[0].gene, p.population[0].fitnessValue))
        if p.population[0].fitnessValue == 0:break

        else:
            print
            print "Evolve using selection, crossover and mutation operators after a generation and add back to population"
            p.evolve()
    else:
        print "Maximum generations reached"
        # print "Fittest String in the given maximum no. of generations:", fittestString, + "with fitness" + p.population[0].fitnessValue
