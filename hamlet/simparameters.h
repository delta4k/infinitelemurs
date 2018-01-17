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
#include "lrandom.h"

/*****************************************************************************\

    Details of the simulation method.
    
    The parameters set as defaults are not necessarily optimal, but they match
    the values used in an early successful experiment, which achieved perfection
    in 784368 generations.

\*****************************************************************************/
class SimParameters
{
public:

    std::string definitionOfFitness;

    // You'll should make this value a multiple of 4 because of the
    // way we manage the population.
    int populationSize = 1024;
    LRandom::RandomAlgorithm randomAlgorithm = LRandom::MERSENNE_TWISTER;
    int numThreads = 4;
    int maxGenerations = INT_MAX - 1;


    /* 
        Mutation rate is the probability that single gene (i.e. a character)
        will be randomized in a new generation.
     
        If evolveMutationRate is true, then the mutation rate itself is subject
        to selection forces. Otherwise, populationMutationRate is applied to
        all creatures.
    */
    bool evolveMutationRate = false;

    double populationMutationRate = 0.00001;

    // If this is true, then individuals will mate with partners
    // of similar fitness.
    bool matesAreSelective = false;

    // WARNING: To simplify the algorithm, we assume that the number
    // of survivors and the number of couples satisfy these constraints.
    // You'll hit assertions otherwise.
    int getNumSurvivors() const { return populationSize / 2; };
    int getNumCouples() const { return populationSize / 4; }
};

/*****************************************************************************\
    Output of the simulation    
\*****************************************************************************/
class SimResult
{
public:
    // All times are in seconds
    double populationCreationTime = 0;// Time it took to create the initial population
    double lastGenerationTime = 0; // computation time for the most recent generation
    double totalTime = 0;
    int numGenerations = 0;

    bool foundNewAllTimeBest = false; // Did the latest generation find a new best fitness?
    int currentMaxFitness = 0; // Best fitness in the latest generation
    int allTimeMaxFitness = 0;
    std::string allTimeBestResult; // Copy of the gene of the best creature so far

    // We store the stats of evolution over time here...
    std::vector<double> normalizedAverageFitness;
    std::vector<double> averageMutationRate; // Only useful if evolveMutationRate is true
    std::vector<int> bestFitness;
    std::vector<double> normalizedBestFitness;

    void printReport(bool to_file, bool brief, const SimParameters& parms);
    void generateSpreadsheet();
    std::string printStats(int generation, const SimParameters& parms);
};
