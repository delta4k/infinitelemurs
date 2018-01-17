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

#include "utilities.h"
#include "simparameters.h"
#include "lemur.h"

/*****************************************************************************\
 
    Print summary information about the current state of the simulation
    to the console.

    If to_file is true, open and print the report to a file named "report.txt"
    in the current working directory.

    If brief is true, we generate an abbreviated report, used for the benchmarking
    mode.
 
\*****************************************************************************/
void SimResult::printReport(bool to_file, bool brief, const SimParameters& parms)
{
    std::vector<std::string> report;

    if(brief)
    {
        // Generate a brief report just showing the timing. Useful for documenting
        // experiments.
        report.push_back(formatString("\n// Random Number Algorithm: %s\n", LRandom::enumLabel(parms.randomAlgorithm)));
        report.push_back(formatString("Population Creation Time: %g\n", populationCreationTime));
        report.push_back(formatString("Last Generation Time: %g\n", lastGenerationTime));
        report.push_back(formatString("Total Time: %g\n", totalTime));

    } else
    {
        report.push_back(formatString("\n\nSummary Report\n"));
        report.push_back(formatString("--------------\n"));
        report.push_back(formatString("Parameters...\n"));
        report.push_back(formatString("Population Size: %d\n", parms.populationSize));
        report.push_back(formatString("Selective Mates: %s\n", parms.matesAreSelective ? "true" : "false"));
        report.push_back(formatString("Text Length: %d\n", Lemur::numGenes));
        report.push_back(formatString("Random Number Algorithm: %s\n", LRandom::enumLabel(parms.randomAlgorithm)));
        report.push_back(formatString("Number of Threads: %d\n", parms.numThreads));

        if(parms.evolveMutationRate)
        {
            report.push_back(formatString("Evolving Mutation Rate"));
        } else
        {
            report.push_back(formatString("Population Mutation Rate: %g\n", parms.populationMutationRate));
        }

        //////////////////////////////////////////////////
        // Results

        report.push_back(formatString("\nResults...\n"));
        report.push_back(formatString("Num Generations: %d\n", numGenerations));
        report.push_back(formatString("Population Creation Time: %g\n", populationCreationTime));
        report.push_back(formatString("Last Generation Time: %g\n", lastGenerationTime));
        report.push_back(formatString("Total Time: %g\n", totalTime));
        report.push_back(formatString("All Time Max Fitness: %d of a possible %d\n", allTimeMaxFitness, Lemur::numGenes));
        // We have averages for the zeroth generation, plus an average at the end of
        // every completed generation.
        Assert(normalizedAverageFitness.size() == numGenerations + 1);

        // Don't print a ridiculous number of lines
        int maxlines = 100;
        int step = numGenerations / maxlines;
        if(step == 0) step = 1;

        int gen;
        for(gen = 0; gen <= numGenerations; gen += step)
        {
            report.push_back(printStats(gen, parms));
        }

        // Be sure we print out the last datapoint
        report.push_back("\nFinal...\n");
        report.push_back(printStats(numGenerations, parms));
    }

    FILE* fp = nullptr;
    if(to_file)
    {
        fp = fopen("report.txt", "wt");
    }

    for(std::string& line : report)
    {
        printf("%s", line.c_str());
        if(fp)
        {
            fprintf(fp, "%s", line.c_str());
        }
    }

    if(fp)
    {
        fclose(fp);
    }
}

/*****************************************************************************\ 
 
        Generate a csv file of the results
 
\*****************************************************************************/
void SimResult::generateSpreadsheet()
{
    std::vector<std::string> rows;

    // Don't print a ridiculous number of lines
    int maxlines = 100;
    int step = numGenerations / maxlines;
    if(step == 0) step = 1;

    // Add a row of column labels
    rows.push_back(formatString("Generation#,Avg. Fitness,Best Fitness\n"));
    int gen;
    bool printed_last_datapoint = false;
    for(gen = 0; gen <= numGenerations; gen += step)
    {
        if (gen == numGenerations) printed_last_datapoint = true;
        rows.push_back(formatString("%d,%g,%g\n", gen, normalizedAverageFitness[gen], normalizedBestFitness[gen]));
    }
    // Be sure we print out the last datapoint
    if (!printed_last_datapoint)
    {
        gen = numGenerations;
        rows.push_back(formatString("%d,%g,%g\n", gen, normalizedAverageFitness[gen], normalizedBestFitness[gen]));
    }

    FILE* fp = fopen("report.csv", "wt");

    if(fp)
    {
        for(std::string& line : rows)
        {
            fprintf(fp, "%s", line.c_str());
        }
        fclose(fp);
    }
}

/*****************************************************************************\
 
        Print out some summary information about the result of the
        simulation in a given generation.
 
\*****************************************************************************/
std::string SimResult::printStats(int generation, const SimParameters& parms)
{
    std::string result;
    if(parms.evolveMutationRate)
    {
        result = formatString("Gen %d, avg. mutation rate %g, best fitness %d, normalized avg. fitness %g\n",
                        generation,
                        averageMutationRate[generation],
                        bestFitness[generation],
                        normalizedAverageFitness[generation]);
    } else
    {
        result = formatString("Gen %d, best fitness %d, normalized avg. fitness %g\n",
                        generation,
                        bestFitness[generation],
                        normalizedAverageFitness[generation]);
    }

    return result;
}
