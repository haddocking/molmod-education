# Python3 program to create target string, starting from
# random string using Genetic Algorithm
# Taken from: https://www.geeksforgeeks.org/genetic-algorithms/
import random

# Target string to be generated
TARGET = "Structural bioinformatics is the best master course!"

# Number of individuals in each generation
POPULATION_SIZE = 200

# Random seed for reproducibility
SEED = 100

# Per-gene mutation probability during mating.
# (The original GeeksforGeeks code hard-coded 0.10, which converges very slowly --
#  a long tail near the optimum; 0.05 converges in a few hundred generations.)
MUTATION_RATE = 0.05

# Optional swap-mutation probability per child: exchange two gene positions.
# Swap only re-orders genes that are already present, so it can NEVER converge on its
# own (it cannot introduce a missing character); it only makes sense in combination
# with point mutation. 0 = off.
SWAP_RATE = 0.0

# Valid genes (single string -- must NOT contain a line break)
GENES = (
    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOP"
    'QRSTUVWXYZ 1234567890, .-;:_!"#%&/()=?@${[]}'
)


class Individual(object):
    """
    Class representing individual in population
    """

    def __init__(self, chromosome):
        self.chromosome = chromosome
        self.fitness = self.cal_fitness()

    @classmethod
    def mutated_genes(cls):
        """
        create random genes for mutation
        """
        return random.choice(GENES)

    @classmethod
    def create_gnome(cls):
        """
        create chromosome or string of genes
        """
        gnome_len = len(TARGET)
        return [cls.mutated_genes() for _ in range(gnome_len)]

    def mate(self, par2):
        """
        Perform mating and produce new offspring
        """

        # each child gene comes from parent 1, parent 2, or a random mutation:
        # mutation happens with probability MUTATION_RATE; the parents split the rest equally
        parent_cut = (1.0 - MUTATION_RATE) / 2.0

        child_chromosome = []
        for gp1, gp2 in zip(self.chromosome, par2.chromosome):
            prob = random.random()

            # gene from parent 1
            if prob < parent_cut:
                child_chromosome.append(gp1)

            # gene from parent 2
            elif prob < 2 * parent_cut:
                child_chromosome.append(gp2)

            # otherwise mutate (random gene), for maintaining diversity
            else:
                child_chromosome.append(self.mutated_genes())

        # optional swap mutation: exchange two gene positions. This only re-orders the
        # genes already present (it cannot introduce a missing character), so it never
        # converges alone -- it only helps in combination with point mutation.
        if SWAP_RATE > 0 and random.random() < SWAP_RATE:
            i = random.randrange(len(child_chromosome))
            j = random.randrange(len(child_chromosome))
            child_chromosome[i], child_chromosome[j] = (
                child_chromosome[j],
                child_chromosome[i],
            )

        # create new Individual(offspring) using
        # generated chromosome for offspring
        return Individual(child_chromosome)

    def cal_fitness(self):
        """
        Calculate fittness score, it is the number of
        characters in string which differ from target
        string.
        """
        global TARGET
        fitness = 0
        for gs, gt in zip(self.chromosome, TARGET):
            if gs != gt:
                fitness += 1
        return fitness


# Driver code


def main():
    global POPULATION_SIZE

    # seed the random-number generator for reproducibility
    random.seed(SEED)

    # current generation
    generation = 1

    found = False
    population = []

    # create initial population
    for _ in range(POPULATION_SIZE):
        gnome = Individual.create_gnome()
        population.append(Individual(gnome))

    while not found:
        # sort the population in increasing order of fitness score
        population = sorted(population, key=lambda x: x.fitness)

        # if the individual having lowest fitness score ie.
        # 0 then we know that we have reached to the target
        # and break the loop
        if population[0].fitness <= 0:
            found = True
            break

        # Otherwise generate new offsprings for new generation
        new_generation = []

        # Perform Elitism, that mean 10% of fittest population
        # goes to the next generation
        s = int((10 * POPULATION_SIZE) / 100)
        new_generation.extend(population[:s])

        # The remaining 90% are offspring of two parents
        # drawn from the 50 fittest individuals
        s = int((90 * POPULATION_SIZE) / 100)
        for _ in range(s):
            parent1 = random.choice(population[:50])
            parent2 = random.choice(population[:50])
            child = parent1.mate(parent2)
            new_generation.append(child)

        population = new_generation

        print(
            "Generation: {}\tString: {}\tFitness: {}".format(
                generation, "".join(population[0].chromosome), population[0].fitness
            )
        )

        generation += 1

    print(
        "Generation: {}\tString: {}\tFitness: {}".format(
            generation, "".join(population[0].chromosome), population[0].fitness
        )
    )


if __name__ == "__main__":
    main()
