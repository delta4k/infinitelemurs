#pragma once

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

#include <random>


/*****************************************************************************\ 
 
    We provide a choice of random number algorithms. LRandom is the abstract
    base class for all random number generation algorithms.

\*****************************************************************************/
class LRandom
{

public:

    enum RandomAlgorithm
    {
        LINEAR_CONGRUENTIAL,
        MERSENNE_TWISTER,
        MERSENNE_TABLE_LOOKUP,
        NUM_RANDOM_ALGORITHMS
    };

    /* 
        We store a different distribution for every range of random numbers
        that is needed for the simulation. In this way we can avoid using
        hacks, like modulo division that reduces randomness, to generate
        integers in each specific range.
    */ 
    enum DistributionNumber
    {
        POPULATION_SIZE_RNDDIST,
        NUM_SURVIVORS_RNDDIST,
        NUM_COUPLES_RNDDIST,
        GENE_LOCUS_RNDDIST, // Location of a gene on the chromosome
        PER_THREAD_POPULATION_SUBSET_RNDDIST,
        ALLELE_RNDDIST, // Possible alleles for a gene
        DOUBLE_TABLE_RNDDIST, // Only used internally by table lookup algorithm
        NUM_DISTRIBUTIONS
    };

    LRandom() : range {}, algorithm {MERSENNE_TWISTER} {}
    virtual ~LRandom() {}

    virtual void initDistribution(DistributionNumber which, int minvalue, int maxvalue) = 0;
    virtual int randomInt(DistributionNumber which) = 0;
    virtual double randomDouble() = 0;
    bool randomBool() { return randomDouble() > 0.5; }

    // For debugging
    static void test();
    static const char* enumLabel(LRandom::RandomAlgorithm e);

    // Factory function
    static LRandom* createRandomizer(RandomAlgorithm algorithm, int seed);

    RandomAlgorithm algorithm;
    
protected:

    void LRandom::testRandomInt(DistributionNumber which);
    void testRandomDouble();

    // This is used primarily for error checking, to enforce our assumption
    // that the simulation isn't changing the distribution after it's first
    // created.
    std::pair<int, int> range[NUM_DISTRIBUTIONS];

};

/*****************************************************************************\
    Generator based on the standard C++ library

    TODO: If we decide to add more variations of engines, it would probably
    make sense to make this a template class to avoid the subclassing.

\*****************************************************************************/
class StdGenerator  : public LRandom
{
public:

    StdGenerator() {};
    void initDistribution(DistributionNumber which, int minvalue, int maxvalue) override;

protected:

    std::uniform_int_distribution<> mDistributions[NUM_DISTRIBUTIONS];
    std::uniform_real_distribution<> realDistribution { 0, 1 };

};

/*****************************************************************************\
    Use a minimal standard linear congruential engine
\*****************************************************************************/
class StdLinearCongruential : public StdGenerator
{
public:

    explicit StdLinearCongruential(int seed) : mEngine(seed) {}

    int randomInt(DistributionNumber which) override;
    double randomDouble() override;

protected:

    std::minstd_rand mEngine;

};

/*****************************************************************************\
    Use a Mersenne Twister engine
\*****************************************************************************/
class StdMersenneTwister32 : public StdGenerator
{
public:

    explicit StdMersenneTwister32(int seed) : mEngine(seed) {}

    int randomInt(DistributionNumber which) override;
    double randomDouble() override;

protected:

    std::mt19937 mEngine; 
};
/*****************************************************************************\
 
    Create a huge array of random numbers, and just loop through it. This
    generator is optimized for speed.

\*****************************************************************************/
class TableLookupRandom : public LRandom
{

public:

    explicit TableLookupRandom(int seed);

    void initDistribution(DistributionNumber which, int minvalue, int maxvalue) override;
    int randomInt(DistributionNumber which) override;
    double randomDouble() override;

private:

    std::vector<int> randomIntegerTables[LRandom::NUM_DISTRIBUTIONS];
    int randomIntegerTableSizes[LRandom::NUM_DISTRIBUTIONS];

    int mTableIndex[NUM_DISTRIBUTIONS];

    double maxDoubleTableValueInverse = 0;
};






