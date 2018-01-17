/**************************************************************************

    Copyright (c) 2011-2017

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

/**************************************************************************

    This program attempts to reproduce some of Shakespeare's writings
    from randomly generated strings. It does so using techniques that
    mimick natural selection by creating a large population of strings
    and allowing the most "fit" to survive and reproduce.

    This program should depend only on standard C library operations
    for ease of porting and building.

Possible improvements:

    Speed:

    Don't create/destroy each creature. Store them on a freelist.

    Don't do a random number for each mutated character, do mutations by
    randomly deciding how many there will be for this creature, then
    picking random characters to hit.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <conio.h>
// Disable MS double to float truncation warning
#pragma warning (disable : 4305)
// Float to int...
#pragma warning (disable : 4244)
// Use posix names for console io functions
#pragma warning (disable : 4996)

#ifdef _DEBUG
#define ASSERT(condition) if(!(condition)) { printf("Assert failure: %s Line %d\n", __FILE__, __LINE__); _getch(); exit(-1); }
#else
#define ASSERT(condition)
#endif

/**************************************************************************

        The Lemurs are the creatures that carry the genetic material
        and reproduce.
*/
class Lemur
{

public:

    Lemur(int gene_size);   // Produce an initial random gene
    Lemur(char *gene);      // Produce a creature with inherited genetics
    ~Lemur();               // Old must make way for new!

    Lemur *mate_with(const Lemur *mate) const;

    int fitness;    // Number of characters that match
    char *gene;     // My version of Shakespeare. Random in the first generation
                    // Produced by crossover with mutation for subsequent generations
};

/*
    Possible target strings...
*/
extern char *ToBeOrNot;         // Simple test
extern char *FirstSonnet;       // Still pretty simple, but ridiculously unlikely statistically
extern char *LoversComplaint;   // Not quite Hamlet, but a pretty convincing demonstration
static int gene_length;        // Same as the length of the target text

#if 0
/*
    Let's say we're willing to use a gigabyte for our population.
    Let's round off the size of LoversComplaint to 16k. Then our
    total population size will be about 1 gig divided by 16k
*/
    const int POP_SIZE = 0x40000000 / 0x00004000;
    const int NUM_SURVIVORS = POP_SIZE / 2; // How many survive in each generation
    char *shakespeare = LoversComplaint;
#else
    // This was the first set of parameters that worked pretty well for LoversComplaint
    const int POP_SIZE = 1024 * 8;
    const int NUM_SURVIVORS = POP_SIZE / 2;
    //char *shakespeare = ToBeOrNot;
    //char *shakespeare = FirstSonnet;
    char *shakespeare = LoversComplaint;
#endif

static float mutation_rate;     // Probability of a mutation per-character
static char character_set[256]; // The allowed set of characters in the text
                                // array is sized to accommodate all possible bytes if desired
static int character_set_len;   // Actual number of characters in the character set

class Lemur *population[POP_SIZE];
class Lemur *survivors[NUM_SURVIVORS];

static unsigned maxrandom_uint = RAND_MAX;
static float average_fitness_normalized; // Average fitness normalized to the 0-1.0 range
static int current_max_fitness;
static int all_time_max_fitness;
static bool new_all_time_best;  // true if we found a new all-time best in this generation

static char *all_time_best_result;  // Copy of the gene of the best creature so far
static char *current_best_result;
static int num_lemurs;              // Used to detect memory leaks

static unsigned random_uint();
static char random_char();
static float random_fraction();
static float new_mutation_rate(int gen_num);

static void create_random_population();
static void evaluate_fitness();
static void choose_survivors();
static void reproduce();

/**************************************************************************/
int main(int argc, char *argv[])
{
    //printf("<%s>\n", shakespeare);
    gene_length = strlen(shakespeare);
    printf("Target document length: %d\n", gene_length);
    current_best_result = new char[gene_length + 1];
    all_time_best_result = new char[gene_length + 1];


    /*
        Build the array of allowed characters. All the characters in the target
        text must exist in this array. 
    */
    int n = 0;
    character_set[n++] = (char)0x08;    // Tab
    character_set[n++] = (char)0x0A;    // Linefeed
    character_set[n++] = (char)0x0D;    // Carriage Return
    for(char c = 0x20; c < 0x7E; c++)   // All printable characters
    {
        character_set[n++] = c;
    }
    character_set_len = strlen(character_set);

    /*
        Verify that the text doesn't have any non-printing characters
    */
    for(int n = 0; n < gene_length; n++)
    {
        if(strchr(character_set, shakespeare[n]) == NULL)
        {
            printf("Text contains non-printing character: (char)%d\n", shakespeare[n]);
            printf("Hit any key to exit\n");
            getch();
            exit(-1);
        }
    }


    /*
        Seed the random number generator with a constant so we can reproduce
        our results.
    */
    srand(777);
    printf("Creating random population...\n");
    create_random_population();
    int generation = 0;
    /*
        Run until somebody hits the 'q' key to quit.
    */
    int keystroke = 0;
    while(keystroke != 'q' && keystroke != 'Q')
    {
        printf("Starting generation %d...\n", generation);
        mutation_rate = new_mutation_rate(generation);
        new_all_time_best = false;
        evaluate_fitness();
        choose_survivors();
        reproduce();

        /*
            Uncomment this to get less frequent output
        */
        //if(generation % 10 == 0)
        //{
            printf("Generation %d done, current best fitness: %d of a possible %d\n", generation, current_max_fitness, gene_length);
            printf("All-time best %d of a possible %d\n", all_time_max_fitness, gene_length);
        //}

        if(new_all_time_best)
        {
            printf("New all-time best result: \n\n%s\n\n", all_time_best_result);
            printf("New all-time best fitness: %d of a possible %d\n", all_time_max_fitness, gene_length);
            if(all_time_max_fitness == gene_length)
            {
                printf("TaDa!\n");
                break;
            }
        }
        generation++;

        if(kbhit())
        {
            keystroke = getch();
        }
    }

    printf("Hit any key to exit\n");
    getch();
}
/**************************************************************************

    Construct a Lemur with a gene of the given size, which will be
    interpreted as a character string. Internally the size is one
    larger to accomodate a terminator so we can print it.
*/
Lemur::Lemur(int gene_size)
{
    num_lemurs++;
    fitness = 0;
    gene = new char[gene_size + 1];
    gene[gene_size - 1] = 0; // Terminate the string so we can print it;
    for(int n = 0; n < gene_size - 1; n++)
    {
        gene[n] = random_char();
    }
}
/**************************************************************************
    Construct a Lemur from an existing gene, normally the crossed-over genes
    of the parents.
*/
Lemur::Lemur(char *g)
{
    num_lemurs++;
    fitness = 0;
    gene = g;
}

/**************************************************************************

*/
Lemur::~Lemur()
{
    num_lemurs--;
    delete[] gene;
}
/**************************************************************************
    Create a new gene by crossing the genes of this with the mate.
    Introduce random mutations at the current mutation rate.
*/
Lemur *Lemur::mate_with(const Lemur *mate) const
{

    char *moms_gene = gene;
    char *dads_gene = mate->gene;
    char *childs_gene = new char[gene_length + 1];
    childs_gene[gene_length] = 0;

    /*
        Split at two points so that we get half the genetic
        information from each parent.
    */
    int crossover1 = random_fraction() * gene_length / 2;
    int crossover2 = crossover1 + gene_length / 2;

    for(int n = 0; n < crossover1; n++)
    {
        if(random_fraction() < mutation_rate)
            childs_gene[n] = random_char();
        else
            childs_gene[n] = moms_gene[n];
    }

    for(int n = crossover1; n < crossover2; n++)
    {
        if(random_fraction() < mutation_rate)
            childs_gene[n] = random_char();
        else
            childs_gene[n] = dads_gene[n];
    }


    for(int n = crossover2; n < gene_length; n++)
    {
        if(random_fraction() < mutation_rate)
            childs_gene[n] = random_char();
        else
            childs_gene[n] = moms_gene[n];
    }

    Lemur *child = new Lemur(childs_gene);
    return child;

}
/**************************************************************************
    Called to create the original primordial soup of genetic material
*/
static void create_random_population()
{
    for(int n = 0; n < POP_SIZE; n++)
    {
        population[n] = new Lemur(gene_length);
    }
}
/**************************************************************************

    For each creature, calculate its fitness defined as the number of
    characters that match the target string.

    Also calculates the current average and best fitness.

*/
static void evaluate_fitness()
{
    int total_fitness = 0;
    current_max_fitness = 0;

    for(int n = 0; n < POP_SIZE; n++)
    {
        Lemur *lemur = population[n];
        char *gene = lemur->gene;
        int fitness = 0;
        for(int nc = 0; nc < gene_length; nc++)
        {
            if(gene[nc] == shakespeare[nc])
            {
                fitness++;
            }
        }
        if(fitness > current_max_fitness)
        {
            current_max_fitness = fitness;
        }
        lemur->fitness = fitness;
        total_fitness += fitness;
    }

    average_fitness_normalized = float(total_fitness) / POP_SIZE;
    average_fitness_normalized /= gene_length;
}
/**************************************************************************

    Function required by quicksort to sort the population.
*/
static int compare_fn(const void *lemur1, const void *lemur2)
{
    int fit1 = (*(Lemur **)lemur1)->fitness;
    int fit2 = (*(Lemur**)lemur2)->fitness;

    if(fit1 < fit2)
        return 1;
    else if(fit1 == fit2)
        return 0;
    else
        return -1;

}
/**************************************************************************

      Sort the population by fitness, then kill off some fraction of
      them and create a survivors array that will be mated to create
      the next generation.
*/
static void choose_survivors()
{

    qsort(population, POP_SIZE, sizeof(void *), compare_fn);
    strcpy(current_best_result, population[0]->gene);

    if(current_max_fitness > all_time_max_fitness)
    {
        new_all_time_best = true;
        all_time_max_fitness = current_max_fitness;
        strcpy(all_time_best_result, population[0]->gene);
    }

    for(int n = 0; n < NUM_SURVIVORS; n++)
    {
        survivors[n] = population[n];
    }
    /*
        Kill off the previous generation.
    */
    for(int n = NUM_SURVIVORS; n < POP_SIZE; n++)
    {
        delete population[n];
        population[n] = NULL;
    }
}
/**************************************************************************

    Go through the survivors and allow the lucky couples to mate.

*/
static void reproduce()
{

    /*
        Randomize the new population so we won't risk odd selection effects
        associated with mating creatures with identical fitness.
    */
    for(int n = 0; n < NUM_SURVIVORS; n++)
    {
        int n2 = random_uint() % (NUM_SURVIVORS);
        Lemur *temporary = survivors[n];
        survivors[n] = survivors[n2];
        survivors[n2] = temporary;
    }

    int num_new_lemurs = 0;
    for(int nparent = 0; nparent < NUM_SURVIVORS - 1; nparent += 2)
    {

        Lemur *mom = survivors[nparent];
        Lemur *dad = survivors[nparent + 1];
        /*
            We need multiple children from each couple to fully repopulate.
        */
        population[num_new_lemurs++] = mom->mate_with(dad);
        population[num_new_lemurs++] = mom->mate_with(dad);
        population[num_new_lemurs++] = mom->mate_with(dad);
        population[num_new_lemurs++] = mom->mate_with(dad);
    }

    ASSERT(num_new_lemurs == POP_SIZE);

    /*
        Make way for the new!
    */
    for(int n = 0; n < NUM_SURVIVORS; n++)
    {
        delete survivors[n];
    }
    ASSERT(num_lemurs == POP_SIZE);
}
/*****************************************************************************

    Produce a mutation rate that changes over time.

    This is probably realistic in the sense that in the course of
    evolution on the Earth the environmental factors have changed
    (probably less radioactivity now) and the creatures themselves
    will have become more or less subject to mutation as the
    molecular biology of reproduction changed. The chemicals in the
    primoridal soup probably had much larger mutation rates than
    current creatures have.
*/
static float new_mutation_rate(int gen_num)
{
    if(gen_num < 50)
        return 0.1;
    else if(gen_num < 100)
        return 0.05;
    else if(gen_num < 200)
        return 0.02;
    else if(gen_num < 300)
        return 0.01;
    else if(gen_num < 1000)
        return 0.001;
    else
        return 0.0001;
}
/**************************************************************************
    Create a random unsigned integer between 0 and maxrandom_uint
*/
static unsigned random_uint()
{
#if 0
    /*
       A pretty good fast random number generator. Linear congruential.
       See Knuth.
    */
    static unsigned maxrandom_uint = 1771875;
    static unsigned a = 2416;
    static unsigned b = 374441;
    static unsigned rndseed = 777;
    rndseed =  (a * rndseed + b ) % maxrandom_uint;
    return rndseed;
#else
    /*
        Use the standard library random number generator. The
        downside of this is that 1) it may be a crappy generator and
        2) results may not be reproducible on other compilers.
    */
    return rand();
#endif
}
/**************************************************************************

    Produce a random character by generating a random integer, then
    using it modulo the number of characters to index into the
    character set.
*/
static char random_char()
{
    unsigned rnd = random_uint();
    char c = character_set[rnd % character_set_len];
    return c;
}
/**************************************************************************
    Return a random float between 0.0 and 1.0
*/
static float random_fraction()
{
    float frac = float(random_uint()) / maxrandom_uint;
    return frac;
}

