
/*****************************************************************************\

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

\*****************************************************************************/

#include <stdio.h>
#include <conio.h>
#include <future>

#include "utilities.h"
#include "simparameters.h"
#include "lemur.h"
#include "simfunctions.h"


/*****************************************************************************\ 
 
    This module contains the actual meat of the simulation.
 
\*****************************************************************************/

SimParameters simParameters;
SimResult simResult;

/***********************************************\
    Tiny class to track the most-fit creature.
\***********************************************/
class BestFitness
{
public:
    bool operator>(const BestFitness& rhs) const { return fitness > rhs.fitness; }
    BestFitness& operator=(const BestFitness&) = default;

    int fitness = 0;
    Lemur* theMostFitCreature = nullptr;
};

class Parents
{
public:
    Lemur* mom = nullptr;
    Lemur* dad = nullptr;
};

static int numThreads;

// Give each thread a different generator to avoid the need
// for synchronization
static std::vector<LRandom*> randomNumberGenerators;

static std::vector<Lemur*> population;
static std::vector<Lemur*> survivors;
static std::vector<Parents> parents;

static void buildRandomNumberGenerators();
static void createRandomPopulation();
static int createRandomPopulationThread(int threadnum);

static void evaluateFitness();
static BestFitness evaluateFitnessThread(int threadnum);

static void chooseSurvivors();
static void reproduce();
static int reproduceThread(int threadnum);

static void mutate();
static int individualMutateThread(int threadnum);
static int populationMutateThread(int threadnum);

static void evaluateAverageStats();
static int numUnmatchedLocations();
static void range(int threadnum, int num_elements, int& from, int& to);

/*****************************************************************************\
 
    Set up our data structures according to the parameters which are assumed
    to have been poked into simParameters.
 
\*****************************************************************************/
void initSimulation()
{

    // Our algorithm operates on even boundaries for simplicity. This enforces it.
    Assert((simParameters.getNumSurvivors() & 1) == 0);
    Assert((simParameters.getNumCouples() & 1) == 0);

    numThreads = simParameters.numThreads;
    Lemur::numGenes = static_cast<int>(simParameters.definitionOfFitness.length());

    buildRandomNumberGenerators();

    printf("Creating random population...\n");
    Stopwatch stopwatch;
    stopwatch.start();
    createRandomPopulation();

    simResult.populationCreationTime = stopwatch.elapsedSeconds();
    simResult.totalTime += simResult.populationCreationTime;
}


/*****************************************************************************\ 
 
    If we build a random number generator for each thread then each thread is
    completely independent, and we don't require synchronization except for the
    freelist operations.
 
\*****************************************************************************/

static void buildRandomNumberGenerators()
{
    For(n, numThreads)
    {
        /*
            Seed the random number generators with constants rather than time,
            so that the results are reproducible. Give each thread a different
            generator with a different seed.
        */
        uint32_t seed = (n + 1) * 77777;  // I like 7's. They're lucky.
        LRandom* lrand = LRandom::createRandomizer(simParameters.randomAlgorithm, seed);

        lrand->initDistribution(LRandom::POPULATION_SIZE_RNDDIST, 0, simParameters.populationSize);
        lrand->initDistribution(LRandom::NUM_SURVIVORS_RNDDIST, 0, simParameters.getNumSurvivors());
        lrand->initDistribution(LRandom::NUM_COUPLES_RNDDIST, 0, simParameters.getNumCouples());
        lrand->initDistribution(LRandom::GENE_LOCUS_RNDDIST, 0, Lemur::numGenes);

        /* 
            Cover the ascii character set. The characters before 9 are old
            teletype control characters. We'll still pick up some silly stuff
            like device control characters, which will give us more trash in
            the genome, but the algorithm should be robust enough to handle
            that, it just strengthens the argument for the effectiveness of
            selection.
             
        */ 
        lrand->initDistribution(LRandom::ALLELE_RNDDIST, 9, 127);

        /* 
            For simplicity in the table lookup algorithm, just implement
            doubles via an integer and a division. This line has no effect
            for the std:: generators.
        */ 
        lrand->initDistribution(LRandom::DOUBLE_TABLE_RNDDIST, 0, 1024 * 1024);

        randomNumberGenerators.push_back(lrand);
    }
}
/******************************************************************************\

    Called to create the original primordial soup of genetic material.
 
    Initialization is expensive for large populations, so we run in
    multiple threads.

\******************************************************************************/
static void createRandomPopulation()
{
    population.resize(simParameters.populationSize);
    survivors.resize(simParameters.getNumSurvivors());
    parents.resize(simParameters.getNumCouples());

    std::vector<std::future<int>> results;
    For(t, numThreads)
    {
        results.push_back(std::async(std::launch::async, createRandomPopulationThread, t));
    }
    For(t, numThreads)
    {
        int actual = results[t].get();
        Assert(actual == t);
    }

    // Evaluate and record the initial population stats
    evaluateFitness();
    evaluateAverageStats();
}

/***********************************************\
    Run a thread to handle a subset of the
    population
\***********************************************/
static int createRandomPopulationThread(int threadnum)
{
    LRandom* generator = randomNumberGenerators[threadnum];
    int from;
    int to;
    range(threadnum, simParameters.populationSize, from, to);
    for(int n = from; n < to; n++)
    {
        population[n] = Lemur::acquire();
        population[n]->randomize(generator);
    }

    return threadnum;
}

/*****************************************************************************\
    Run the simulation for one generation.
\*****************************************************************************/

void runOneGeneration()
{
    if (simResult.allTimeMaxFitness == Lemur::numGenes)
    {
        printf("Perfection Achieved!\n");
        return;
    }

    printf("Generation %d...\n", simResult.numGenerations);

    Stopwatch stopwatch;

    stopwatch.start();

    simResult.numGenerations++;

    prn("Choosing survivors...\n");
    chooseSurvivors();
    prn("Reproducing the next generation...\n");
    reproduce();
    prn("Mutating...\n");
    mutate();

    prn("Evaluating fitness...\n");
    evaluateFitness();
    evaluateAverageStats();

    simResult.lastGenerationTime = stopwatch.elapsedSeconds();
    simResult.totalTime += simResult.lastGenerationTime;

    prn("Generation %d: Max Fitness %d/%d, execution time =  %g seconds\n",
        simResult.numGenerations, simResult.currentMaxFitness, Lemur::numGenes, simResult.lastGenerationTime);
    prn("All-time best %d of a possible %d\n", simResult.allTimeMaxFitness, Lemur::numGenes);

}

/*****************************************************************************\
    Compute some statistical information for reporting.
\*****************************************************************************/
static void evaluateAverageStats() 
{
    double normalized_ave_fitness = 0;
    double ave_mutrate = 0;
    for (auto creature : population)
    {
        normalized_ave_fitness += double(creature->fitness) / Lemur::numGenes;
        ave_mutrate += creature->mutationRate;
    }

    normalized_ave_fitness = normalized_ave_fitness / simParameters.populationSize;
    ave_mutrate = ave_mutrate / simParameters.populationSize;
    simResult.normalizedAverageFitness.push_back(normalized_ave_fitness);
    simResult.averageMutationRate.push_back(ave_mutrate);
    simResult.bestFitness.push_back(simResult.currentMaxFitness);
    simResult.normalizedBestFitness.push_back(double(simResult.currentMaxFitness) / Lemur::numGenes);
}

/**************************************************************************\ 
 
    For each creature, calculate its fitness defined as the number of
    characters that match the definitionOfFitness string.

\**************************************************************************/
static void evaluateFitness()
{
    simResult.currentMaxFitness = 0;
    simResult.foundNewAllTimeBest = false;

    std::vector<std::future<BestFitness>> results;

    For(t, numThreads)
    {
        results.push_back(std::async(std::launch::async, evaluateFitnessThread, t));
    }

    BestFitness best;

    For(t, numThreads)
    {
        BestFitness bf;
        bf = results[t].get();
        if(bf > best)
        {
            best = bf;
        }
    }

    simResult.currentMaxFitness = best.fitness;

    if(simResult.currentMaxFitness > simResult.allTimeMaxFitness)
    {
        simResult.foundNewAllTimeBest = true;
        simResult.allTimeMaxFitness = simResult.currentMaxFitness;
        simResult.allTimeBestResult = best.theMostFitCreature->chromosome;
    }

}

/***********************************************\
   Run a thread to handle a subset of the
   population
\***********************************************/
static BestFitness evaluateFitnessThread(int threadnum)
{
    // Range of individuals this thread operates on...
    int from;
    int to;
    range(threadnum, simParameters.populationSize, from, to);

    BestFitness best;

    // Doing reads via c_str() has a significant performance advantage.
    const char* perfection = simParameters.definitionOfFitness.c_str();
    for(int n = from; n < to; n++)
    {
        Lemur* lemur = population[n];
        Assert(lemur);

        const char* chrom = lemur->chromosome.c_str();
        int fitness = 0;
        for(int nc = 0; nc < Lemur::numGenes; nc++)
        {
            if(chrom[nc] == perfection[nc])
            {
                fitness++;
            }
        }
        lemur->fitness = fitness;

        if(fitness > best.fitness)
        {
            best.fitness = fitness;
            best.theMostFitCreature = lemur;
        }
    }

    return best;
}

/*****************************************************************************\
    
    Sort the population by fitness, then kill off some fraction of them and
    create a survivors array that will be mated to create the next generation.

\*****************************************************************************/
static void chooseSurvivors()
{
    int popsize = simParameters.populationSize;
    int numsurvivors = simParameters.getNumSurvivors();

    std::sort(population.begin(), population.end(), [](auto l0, auto l1) { return l0->fitness > l1->fitness; });
    Assert(population[0]->fitness >= population[numsurvivors]->fitness);

    // Take the best as survivors
    for(int n = 0; n < numsurvivors; n++)
    {
        survivors[n] = population[n];
    }

    //  Kill off the previous generation
    for(int n = numsurvivors; n < popsize; n++)
    {
        Lemur::release(&population[n]);
    }

    // Do this to allow assertions in reproduction routines...
    For(n, popsize)
    {
        population[n] = nullptr;
    }
}

/**************************************************************************\

    Go through the survivors and allow the lucky couples to mate.

\**************************************************************************/
static void reproduce()
{
    /*
        If matesAreSelective is true, then pairs will be joined with mates of
        similar fitness. Otherwise we randomize the survivors so there is
        no correlation between fitness of mates.
    */
    if(!simParameters.matesAreSelective)
    {
        /*
            We don't use std::shuffle here because our random number algorithm
            is not necessarily an std engine.
        */
        LRandom* generator = randomNumberGenerators[0];
        for(int n = 0; n < simParameters.getNumSurvivors(); n++)
        {
            int n2 = generator->randomInt(LRandom::NUM_SURVIVORS_RNDDIST);
            Swap<Lemur*>(survivors[n], survivors[n2]);
        }
    }

    for(int n = 0; n < simParameters.getNumCouples(); n++)
    {
        Assert((n * 2 + 1) < simParameters.getNumSurvivors());
        parents[n].mom = survivors[n * 2];
        parents[n].dad = survivors[n * 2 + 1];
    }

    std::vector<std::future<int>> results;
    For(t, numThreads)
    {
        results.push_back(std::async(std::launch::async, reproduceThread, t));
    }

    For(t, numThreads)
    {
        /* 
            Empirically, if we don't force a result from the async thread, then
            the async mechanism doesn't work as expected. That's why we have
            some code here to prevent the optimizer from messing this up.
        */ 
        int actual = results[t].get();
        Assert(actual == t);
    }

    //  Make way for the new!
    for(int n = 0; n < simParameters.getNumSurvivors(); n++)
    {
        Lemur::release(&survivors[n]);
    }
}

/***********************************************\
    Run a thread to handle a subset of the
    population
\***********************************************/
static int reproduceThread(int threadnum)
{
    LRandom* generator = randomNumberGenerators[threadnum];

    // The range of parent couples that we will process
    int from;
    int to;
    range(threadnum, simParameters.getNumCouples(), from, to);

    // The range of slots in the population vector that this
    // thread will fill in
    int popfrom;
    int popto;
    range(threadnum, simParameters.populationSize, popfrom, popto);
    int popslot = popfrom;
    for(int n = from; n < to; n++)
    {
        Lemur* mom = parents[n].mom;
        Lemur* dad = parents[n].dad;

        // Each couple produces 4 offspring
        Assert(population[popslot] == nullptr);
        population[popslot++] = mom->mateWith(dad, generator);

        Assert(population[popslot] == nullptr);
        population[popslot++] = mom->mateWith(dad, generator);

        // Reverse the order of mom and dad so we won't be biased  in terms of
        // which parent comes first in the mateWith algorithm
        Assert(population[popslot] == nullptr);
        population[popslot++] = dad->mateWith(mom, generator);

        Assert(population[popslot] == nullptr);
        population[popslot++] = dad->mateWith(mom, generator);
    }
    Assert(popslot == popto);
    return threadnum;
}


/*****************************************************************************\
 
    Select a mutation algorithm, whole-population or individual, based on
    the SimulationParameters

\*****************************************************************************/
static void mutate()
{
    std::vector<std::future<int>> results;
    For(t, numThreads)
    {
        if(simParameters.evolveMutationRate)
        {
            results.push_back(std::async(std::launch::async, individualMutateThread, t));
        }else{
            results.push_back(std::async(std::launch::async, populationMutateThread, t));
        }
    }
    For(t, numThreads)
    {
        int actual = results[t].get();
        Assert(actual == t);
    }
}

/***********************************************\ 
 
  Use each creature's individual mutation rate. This
  allows the mutation rate itself to be subject
  to natural selection. It also means our performance
  is lower, because of the need to run the random
  number generator for every individual.
 
\***********************************************/
static int individualMutateThread(int threadnum)
{

    // Range of individuals this thread operates on...
    int from;
    int to;
    range(threadnum, simParameters.populationSize, from, to);
    LRandom* generator = randomNumberGenerators[threadnum];

    for (int n = from; n < to; n++)
    {
        population[n]->mutate(generator);
    }
    return threadnum;
}
/***********************************************\ 
 
  Use a single mutation rate for the entire population.
 
  Compute the expected number of mutations given the current mutation rate and
  number of creatures, then randomly change that number of "genes", where a
  gene is a single character in the chromosome.
 
\***********************************************/
static int populationMutateThread(int threadnum)
{
    // If the rate is exactly zero, then we are testing the total absence
    // of mutation. This is used only for test, to verify that a population
    // with unmatched genes will never converge without mutation.
    if (simParameters.populationMutationRate == 0)
    {
        return threadnum;
    }

    int from;
    int to;
    range(threadnum, simParameters.populationSize, from, to);
    LRandom* generator = randomNumberGenerators[threadnum];

    int num_lemurs = to - from;
    double d_num_alleles = static_cast<double>(num_lemurs) * Lemur::numGenes;
    double d_num_mutations = d_num_alleles * simParameters.populationMutationRate;
    int64_t num_mutations = static_cast<int64_t>(d_num_mutations);
    /*
        For very small populations and low mutation rates, this can round off
        to zero, which is not realistic. There will always be SOME finite amount
        of mutation.
    */
    num_mutations = Max(num_mutations, 1);

//    if(threadnum == 0)
//    {
        // Generate some sanity-check output
        //prn("Thread %d Generating %g mutations in %d individuals and %g alleles\n",
        //       threadnum, d_num_mutations, num_lemurs, d_num_alleles);
//    }

    generator->initDistribution(LRandom::PER_THREAD_POPULATION_SUBSET_RNDDIST, from, to);
    for(int64_t n = 0; n < num_mutations; n++)
    {
        int which_individual = generator->randomInt(LRandom::PER_THREAD_POPULATION_SUBSET_RNDDIST);
        Assert(Between(from, which_individual, to - 1));
        Lemur* lemur = population[which_individual];

        int which_gene = generator->randomInt(LRandom::GENE_LOCUS_RNDDIST);
        lemur->chromosome[which_gene] = generator->randomInt(LRandom::ALLELE_RNDDIST);
    }
    return threadnum;
}

/*****************************************************************************\
 
    Return a range of indices on which the given thread can operate without 
    colliding with other threads.
 
\*****************************************************************************/

static void range(int threadnum, int num_elements, int& from, int& to)
{
    int span = num_elements / numThreads;
    from = threadnum * span;
    if(threadnum == numThreads - 1)
    {
        to = num_elements;
    } else
    {
        to = (threadnum + 1) * span;
    }
}

/*****************************************************************************\

    The checksum is a simple one for use in regression testing during
    refactoring, to verify we haven't changed the algorithm.

\*****************************************************************************/
int populationCheckSum()
{
    int csum = 0;
    for(int n = 0; n < simParameters.populationSize; ++n)
    {
        csum += checksum(population[n]->chromosome);
    }
    return csum;
}

/*****************************************************************************\
       
   Find the number of character locations at which no member of the
   population has a match. If this number is non-zero, then it will
   be impossible to achieve perfection without mutation.

   This is actually a rather expensive calculation, probably because it blows
   the cache in the inner loop. It was only used to verify the fact that
   mutation is necessary in populations with unmatched locations.

\*****************************************************************************/
static int numUnmatchedLocations()
{
    int num_unmatched = 0;

    For(n, Lemur::numGenes)
    {
        char target_char = simParameters.definitionOfFitness[n];
        bool matched = false;
        for (Lemur* lemur : population)
        {
            if (lemur->chromosome[n] == target_char)
            {
                matched = true;
                break;
            }
        }

        if (!matched)
        {
            ++num_unmatched;
        }
    }
    return num_unmatched;
}

