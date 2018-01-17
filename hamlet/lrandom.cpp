
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

#include <stdlib.h>
#include <set>
#include <algorithm>
#include <vector>

#include "utilities.h"
#include "lrandom.h"

// In theory, as long as this is large enough to hold one copy of 
// each value, then it's big enough. We guarantee that below.
static int minRandomIntegerTableSize = 2 * 1024 * 1024;

/*****************************************************************************\
    
    Create a distribution that will generate integers on the half-open
    interval [min, max)

\*****************************************************************************/
void StdGenerator::initDistribution(DistributionNumber which, int minvalue, int maxvalue)
{
    if(range[which].second > 0)
    {
        // Be sure they haven't changed the range for this distribution
        Assert(range[which].first == minvalue);
        Assert(range[which].second == maxvalue);
    } else
    {
        // WARNING: std generates integers on the closed interval. That's why
        // we use max-1 here
        mDistributions[which] = std::uniform_int_distribution<>(minvalue, maxvalue - 1);
        range[which].first = minvalue;
        range[which].second = maxvalue;
    }
}


/*****************************************************************************\
                        Linear Congruential
\*****************************************************************************/
int StdLinearCongruential::randomInt(DistributionNumber which)
{
    return (mDistributions[which])(mEngine);
}

double StdLinearCongruential::randomDouble()
{
    return realDistribution(mEngine);
}

/*****************************************************************************\
                           Mersenne Twister
\*****************************************************************************/
int StdMersenneTwister32::randomInt(DistributionNumber which)
{
    return (mDistributions[which])(mEngine);
}

double StdMersenneTwister32::randomDouble()
{
    return realDistribution(mEngine);
}

/*****************************************************************************\
                          Table Lookup.
\*****************************************************************************/
TableLookupRandom::TableLookupRandom(int seed)
{
    // The table index is the next place we'll go to look up a random number.
    // Each distribution maintains it's own index, to avoid allowing the number
    // of calls to one distribution biasing the values found by a different
    // distribution.
    for(int i = 0; i < NUM_DISTRIBUTIONS; ++i)
    {
        mTableIndex[i] = seed % minRandomIntegerTableSize;
    }
}

/*****************************************************************************\
    
    Build a table of random numbers for this distribution. It should contain
    each possible value the same number of times as every other value.

\*****************************************************************************/
void TableLookupRandom::initDistribution(DistributionNumber which, int minvalue, int maxvalue)
{
    if(range[which].second > 0)
    {
        // Be sure they haven't changed the range for this distribution
        Assert(range[which].first == minvalue);
        Assert(range[which].second == maxvalue);
        return;
    }

    if(which ==  DOUBLE_TABLE_RNDDIST)
    {
        // This allows us to use a multiplication rather than the slightly
        // more expensive division.
        maxDoubleTableValueInverse = 1.0 / double(maxvalue + 1);
    }

    range[which].first = minvalue;
    range[which].second = maxvalue;

    //////////////////////////////////////////////////
    //  We want every integer in the range to occur with equal frequency,
    //  so we'll size the array to be a multiple of the number of values,
    //  but at least as big as the min table size.

    // The user provides us with an open interval, so there is no +1 here
    int num_values = maxvalue - minvalue;

    Assert(minRandomIntegerTableSize > num_values);

    // Ensure that each integer occurs the same number of times in our table,
    int num_copies = minRandomIntegerTableSize / num_values;
    // Add 1 to be sure we satisfy our minimum size requirement.
    num_copies += 1;
    int table_size = num_copies * num_values;

    Assert(table_size >= minRandomIntegerTableSize);
    Assert((table_size % num_values) == 0);

    std::vector<int>& table = randomIntegerTables[which];
    table.resize(table_size);
    randomIntegerTableSizes[which] = table_size;

    /*
       To guarantee that all values occur and occur with equal
       probability, we fill the array with all integers then
       shuffle it.
    */
    int the_random_int = 0;
    for(int i = 0; i < table_size; ++i)
    {
        the_random_int = minvalue + (i % num_values);
        table.at(i) = the_random_int;
    }
    // Verify that we have generated the complete set of values
    Assert(the_random_int == maxvalue - 1);

    // Shuffle the table with a 64 bit Mersenne twister
    std::mt19937_64 engine; 
    std::shuffle(table.begin(), table.end(), engine);
}

/***********************************************\
      
\***********************************************/
int TableLookupRandom::randomInt(DistributionNumber which)
{
    int idx = mTableIndex[which];
    idx = (idx == randomIntegerTableSizes[which] - 1) ? 0 : idx + 1;
    mTableIndex[which] = idx;

    int value = randomIntegerTables[which][idx];

    return value;
}

/***********************************************\
      
\***********************************************/
double TableLookupRandom::randomDouble()
{
    int value = randomInt(DOUBLE_TABLE_RNDDIST);

    return value * maxDoubleTableValueInverse;
    
}


/*****************************************************************************\
    Randomizer Factory Function
\*****************************************************************************/
LRandom* LRandom::createRandomizer(RandomAlgorithm algorithm, int seed)
{
    LRandom* result;
    switch(algorithm)
    {
        case LRandom::LINEAR_CONGRUENTIAL:
            result = new StdLinearCongruential(seed);
            break;
        case MERSENNE_TWISTER:
            result = new StdMersenneTwister32(seed);
            break;
        case MERSENNE_TABLE_LOOKUP:
            result = new TableLookupRandom(seed);
            break;
        default:
            Assert(0);
            return nullptr;
    }

    result->algorithm = algorithm;

    return result;
}

/*****************************************************************************\
        Generate a textual label for the given algorithm
\*****************************************************************************/
const char* LRandom::enumLabel(RandomAlgorithm e)
{
    switch(e)
    {
        case LINEAR_CONGRUENTIAL:
            return "Linear Congruential";
        case MERSENNE_TWISTER:
            return "Mersenne Twister";
        case MERSENNE_TABLE_LOOKUP:
            return "Mersenne Table Lookup";
        default:
            Assert(0);
            return "Unknown";
    }
}

/*****************************************************************************\
        Perform a sanity check on the random number generation.
\*****************************************************************************/
void LRandom::testRandomInt(DistributionNumber which)
{

    int minvalue = range[which].first;
    int maxvalue = range[which].second;
    int num_values = maxvalue - minvalue;

    std::vector<double> num_hits(maxvalue);

    int num_calls = 10000000;
    For(c, num_calls)
    {
        int nr = randomInt(which);
        Assert(Between(minvalue, nr, maxvalue - 1));
        num_hits[nr]++;
    }

    double expected_frac = 1 / double(num_values);
    printf("Expected fraction %g\n", expected_frac);

    double max_variation = 0;
    for (int n = minvalue; n < maxvalue; n++)
    {
        double frac = num_hits[n] / num_calls;
        max_variation = Max(fractionalDifference(frac, expected_frac), max_variation);
        if(n == 0)
            printf("Sample value %d : fraction %g\n", n, frac);
    }
    printf("Fractional variation in expected value: %g\n", max_variation);
}
void LRandom::testRandomDouble()
{
    printf("Testing randomDouble...\n");

    const int numbins = 100;
    int num_hits[numbins] = {0};

    int num_calls = 10000000;
    For(c, num_calls)
    {
        double rd = randomDouble();
        Assert(Between(0, rd, 1.0));
        /*
            The uniform_real_distribution returns numbers on the open interval,
            unlike the int distributions. That's why we don't have to worry
            about the bin number exceeding numbins-1.
         */
        int bin = int(rd * numbins);
        num_hits[bin]++;
    }

    double expected_frac = 1 / double(numbins);
    printf("Expected fraction %g\n", expected_frac);

    double max_variation = 0;
    for (int n = 0; n < numbins; n++)
    {
        double frac = double(num_hits[n]) / num_calls;
        max_variation = Max(fractionalDifference(frac, expected_frac), max_variation);
        if (n == 0)
            printf("Sample value %d : fraction %g\n", n, frac);
    }
    printf("Fractional variation in expected value: %g\n", max_variation);
}
/*****************************************************************************\
 
    Spits out a bunch of information on the random numbers. Inspect visually.
 
\*****************************************************************************/
void LRandom::test()
{

    for(int i = 0; i < NUM_RANDOM_ALGORITHMS; ++i)
    {
        printf("\n---------------------------------\nTesting algorithm %s\n", enumLabel(RandomAlgorithm(i)));
        LRandom* generator;
        generator = createRandomizer(RandomAlgorithm(i), 7777);
        generator->initDistribution(GENE_LOCUS_RNDDIST, 0, 1000);
        generator->initDistribution(PER_THREAD_POPULATION_SUBSET_RNDDIST, 4096, 8500);
        generator->initDistribution(NUM_COUPLES_RNDDIST, 0, 32);
        generator->initDistribution(ALLELE_RNDDIST, 7, 35);
        // This one is strictly speaking only needed by the table lookup generator
        generator->initDistribution(LRandom::DOUBLE_TABLE_RNDDIST, 0, 1024 * 1024);

        generator->testRandomInt(GENE_LOCUS_RNDDIST);
        generator->testRandomInt(PER_THREAD_POPULATION_SUBSET_RNDDIST);
        generator->testRandomInt(NUM_COUPLES_RNDDIST);
        generator->testRandomInt(ALLELE_RNDDIST);
        generator->testRandomDouble();

    }

}
