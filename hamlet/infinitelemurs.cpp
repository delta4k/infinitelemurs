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
/*****************************************************************************\ 
 
    What does this program do?
    --------------------------
 
    This program reproduces the entire text of Shakespeare's Hamlet from
    randomly generated characters, through a process of natural selection.
 
    The techniques used are described in the classic text:
 
        "Genetic Algorithms in Search, Optimization & Machine Learning"
        by David E. Goldberg, Pub. Addison-Wesley, 1989
 
 
    To quote from the text:
 
        "Genetic algorithms are search algorithms based on the mechanics of
        natural selection and natural genetics. They combine survival of the
        fittest among string structures with a structured yet randomized
        information exchange."
 
    A genetic algorithm requires a genetic code and a definition of "fitness".
    In our case, the genetic code is a string of characters, and fitness is
    defined by the number of characters in the string that match the
    corresponding character in the same location in the text of Hamlet.
 
    It is important to note that every character is generated randomly. The
    only reference to the original text is in the definition of fitness.
 
    The program has a number of free parameters (see the SimParameters
    structure) which allow testing of different details of the algorithm.
 
    What key features of evolution does this simulation demonstrate?
    ----------------------------------------------------------------

    While this is obviously a vastly simplified simulation, it demonstrates
    several key features of evolution...

    1)  Randomization plus selection can yield results that are so unlikely as
        to seem impossible. The probability of producing the phrase "Order from
        chaos!" by randomly combining characters is approximately one in 10 to
        the 33rd power. If you produced a million random trials per second, it
        would take you 10^20 years, or 10 billion times the age of the
        universe, to arrive at this result by pure chance. Calculations like
        this are frequently used as arguments against evolution.
 
    2)  Without mutation, evolution is impossible. A gene pool has a finite
        amount of variation. Mutation introduces new variations that have never
        been seen before. So, for example in our simulation, if no individuals
        in the initial population match a given character of Hamlet in a given
        position, then without mutation none of the descendents would ever
        match that character.

    3)  Without sex, evolution is impossible. You have to mix up the genes to
        produce new variations. There is a caveat in that evolution could
        in theory be based purely on mutation, but that would be
        extremely innefficient (this program could be used to test such a
        method).

    4)  Without death, evolution is impossible, you have to make room for the
        new variations.
 
 
    How does this simulation compare to biological evolution?
    ---------------------------------------------------------
 
    This simulation differs from biological evolution in several ways:
 
    1)  In biological evolution the "fitness function" is survival in the
        environment. There is no specific goal, hence there are many possible
        solutions to the survival problem, each solution being a different
        species.

    2)  Biological evolution is much more "parallel". A modern computer could
        not simulate the complete behaviour of a single atom in real time, much
        less an entire organism, much less a population of organisms interacting
        with their environment.

    3)  Biological evolution operates on much larger populations (e.g. the
        original population consisted of all the molecules in the primaeval ocean).

    4)  Biological evolution operates for a much longer period of time. There
        has been life on earth for something like a third of the entire age
        of the universe.
 
    5)  The details of the crossover algorithm differ, although this is somewhat
        a matter of nomenclature: If you called each character a
        "chromosome" (we call it a "gene") then there would be a closer
        correspondence between our algorithm and the typical biological
        crossover mechanism. However all that really matters is that the
        genetic material of the two parents gets randomly mixed.

\*****************************************************************************/

#include <stdio.h>
#include <conio.h>
#include <list>
#include <iostream>

#include "utilities.h"
#include "lemur.h"
#include "simfunctions.h"

/*****************************************************************************\

    The program can run in several different modes, basically test modes for
    development as opposed to the HAMLET mode which is the final demonstration 
    of the algorithm's capabilities.

\*****************************************************************************/
enum SimMode
{
    HAMLET, // Run the full experiment
    REGRESSION_TEST, // Use when refactoring to verify result hasn't changed
    SCRATCH, // For doing miscellaneous experiments
    BENCHMARK, // Use to test efficiency
    RND_TEST // Use to inspect the random number generator output
};


static SimMode simMode = HAMLET;
//static SimMode simMode = REGRESSION_TEST;
//static SimMode simMode = SCRATCH;
//static SimMode simMode = BENCHMARK;
//static SimMode simMode = RND_TEST;

// If true, run the program under the control of a simple text-mode menu.
static bool useMenu = false;

// Text for a simple, fast regression test...
const char* ToBeOrNot = "To be, or not to be.";

static void runSimToCompletion();
static void runFromMenu();
static void testRandomNumberGenerators();
static void buildSimParameters();

/*****************************************************************************\

    This program runs one complete simulation. There is no way to restart the
    simulation and do another run except by terminating the program and running
    it again. 
    
    The arguments to main() are not used. Experiments are done by changing the
    hardcoded simParameters initialization in buildSimParameters.
 
\*****************************************************************************/
int main(int argc, char* argv[])
{
    buildSimParameters();

    initSimulation();

    if(useMenu)
    {
        runFromMenu();
    } else
    {
        runSimToCompletion();
    }

    // Print a brief report if we are benchmarking. Makes it easier
    // to collect and compare results
    simResult.printReport(true, simMode == BENCHMARK, simParameters);
    simResult.generateSpreadsheet();

    if (simMode == REGRESSION_TEST)
    {
        // This verifies that we haven't broken anything when refactoring.
        int checksum = populationCheckSum();
        Assert(checksum == 0x63B52);
        printf("\nVerified checksum: 0x%X\n", checksum);
    }else if(simMode == BENCHMARK)
    {
        // Assert if we experience a significant regression in performance
        if(simParameters.evolveMutationRate)
        {
            Assert(simResult.lastGenerationTime < 4.0);
        } else
        {
            Assert(simResult.lastGenerationTime < 0.55);
        }

    }
    printf("Hit any key to exit\n");
    getch();
}

/*****************************************************************************\
    Run the program under menu control, rather than just running the simulation
    to completion without user interaction.
\*****************************************************************************/
static void runFromMenu()
{
    while(1)
    {
        printf("\n\n");
        printf("q - %s\n", "Quit");
        printf("1 - %s\n", "Run One Generation");
        printf("2 - %s\n", "Run Continuously");
        printf("3 - %s\n", "Print Current Best Creature");
        printf("4 - %s\n", "Print Report");
        printf("5 - %s\n", "Print Report to file");
        printf("6 - %s\n", verboseOutput ? "Turn Off Verbosity" : "Turn On Verbosity");

        int keystroke = getch();
        switch (keystroke)
        {
        case 'q':
        case 'Q':
            printf("Hit 'q' to confirm exit...\n");
            keystroke = getch();
            if(keystroke == 'q' || keystroke == 'Q')
            {
                exit(0);
            }
            break;
        case '1':
            runOneGeneration();
            break;
        case '2':
            runSimToCompletion();
            break;
        case '3':
            {
                
                printf("\n\nBest Creature So Far:\n");
                // printf may not be able to handle a very long string, so I'll use
                // cout, which I assume won't have such a problem.
                fflush(stdout);
                std::cout << simResult.allTimeBestResult;
                std::cout.flush();
                printf("\nFitness: %d of a possible %d\n", simResult.allTimeMaxFitness, Lemur::numGenes);
            }
            break;
        case '4':
            // Don't generate a report file in this case
            simResult.printReport(false, false, simParameters);
            break;
        case '5':
            simResult.printReport(true, false, simParameters);
            simResult.generateSpreadsheet();
            break;
        case '6':
            verboseOutput = !verboseOutput;
        default:
            break;
        }
    }
}

/*****************************************************************************\
    Create SimulationParameters appropriate for our current SimMode.
\*****************************************************************************/
static void buildSimParameters()
{
    if(simMode == REGRESSION_TEST)
    {
        simParameters.definitionOfFitness = ToBeOrNot;
        simParameters.populationSize = 256;
        simParameters.numThreads = 4;
        simParameters.evolveMutationRate = false;

        // With a population of only 256, we verify that there are unmatched
        // locations and that mutation is necessary.
        // Uncomment this line, and you will never match the definitionOfFitness.
        //simParameters.populationMutationRate = 0;
        useMenu = false;

    } else if(simMode == SCRATCH)
    {
        simParameters.definitionOfFitness = ToBeOrNot;
        simParameters.populationSize = 256;
        simParameters.numThreads = 4;
        simParameters.evolveMutationRate = false;
        simParameters.populationMutationRate = 0;
        useMenu = true;
    } else if(simMode == BENCHMARK)
    {
        const char* fname = "hamlet.txt";
        simParameters.definitionOfFitness = loadText(fname);
        simParameters.populationSize = 1024 * 8;
        simParameters.maxGenerations = 5;

        simParameters.numThreads = 4;
        verboseOutput = false;

        //        simParameters.randomAlgorithm = LRandom::MERSENNE_TABLE_LOOKUP;
    } else if(simMode == HAMLET)
    {
        //////////////////////////////////////////////////
        // Large population sizes gain fitness in fewer generations,
        // but take more computation time. In the real world, populations
        // run in parallel, so the tradeoff isn't there.

        const char* fname = "hamlet.txt";
        simParameters.definitionOfFitness = loadText(fname);

        //////////////////////////////////////////////////
        // Put parameters here that deviate from the original 
        // successful parameters (that are now the default).


        //////////////////////////////////////////////////
        verboseOutput = false;
        useMenu = true;
    } else if(simMode == RND_TEST)
    {
        // This procedure doesn't return
        testRandomNumberGenerators();
    } else
    {
        Assert(0);
    }
}

/******************************************************************************

    Run the simulation until we achieve perfection, or are stopped by the user,
    or until the maximum number of generations expires.

******************************************************************************/
static void runSimToCompletion()
{
    int keystroke = 0;
    while(simResult.numGenerations < simParameters.maxGenerations)
    {
        // Don't do anything if we've already reached perfection
        if(simResult.allTimeMaxFitness == Lemur::numGenes)
        {
            printf("Perfection Achieved!\n");
            break;
        }

        runOneGeneration();

        if(simResult.foundNewAllTimeBest)
        {
            printf("New all-time best fitness: %d of a possible %d\n", simResult.allTimeMaxFitness, Lemur::numGenes);
        }

        // kbhit is windows-specific, but isn't hard to emulate if you want to
        // port to linux or mac.
        if(kbhit())
        {
            // Any keystroke will cause us to pause.
            getch();
            printf("Hit 'q' to break out of run loop, any other key to continue...\n");
            keystroke = getch();

            // If we are running under menu control, this will throw us back
            // into the menu, otherwise the program will terminate after
            // printing a report.
            if(keystroke == 'q' || keystroke == 'Q')
            {
                break;
            }
        }
    }
}


/*****************************************************************************\
\*****************************************************************************/
static void testRandomNumberGenerators()
{
    LRandom::test();

    printf("Hit any key to exit...\n");
    getch();
    exit(0);
}

