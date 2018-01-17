
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

#include <cstdio>
#include <atomic>
#include <mutex>

#include "utilities.h"
#include "lemur.h"

int Lemur::numGenes; // Length of the text

// Use a freelist for memory management. Major performance boost.
static Lemur* freelist;

// The freelist is accessed from multiple threads.
static std::mutex freeListMutex;

/*****************************************************************************\
 
    Retrieve a Lemur from the pool of souls, unless there are none, in which
    case we do a conventional construction.
 
\*****************************************************************************/
Lemur* Lemur::acquire()
{
    Lemur* result;

    {
        std::lock_guard<std::mutex> lock{freeListMutex};
        if(freelist)
        {
            result = freelist;
            freelist = freelist->next;
            result->fitness = 0;
            result->next = nullptr;
        } else
        {
            result = new Lemur();
        }
    }


    return result;
}

/*****************************************************************************\
 
    Release a Lemur back into limbo to await rebirth in the next
    generation.
 
\*****************************************************************************/
void Lemur::release(Lemur** lemur)
{
    std::lock_guard<std::mutex> lock{freeListMutex};
    (*lemur)->next = freelist;
    freelist = *lemur;
    *lemur = nullptr;
}

/*****************************************************************************\
 
    Construct a Lemur with a gene of the given size, which will be
    interpreted as a character string.
 
\*****************************************************************************/
Lemur::Lemur()
{
    next = nullptr;
    mutationRate = 0;
    fitness = 0;

    Assert(numGenes >0);
    chromosome.resize(numGenes);
}

/*****************************************************************************\
 
    When the population is initially created, we give each creature a totally
    random chromosome.
 
\*****************************************************************************/
void Lemur::randomize(LRandom* generator)
{

    for(int n = 0; n < numGenes; n++)
    {
        chromosome[n] = generator->randomInt(LRandom::ALLELE_RNDDIST);
    }
    mutationRate = generator->randomDouble();
}

/*****************************************************************************\

    For the mode where we are letting the mutation rate itself evolve, we
    use this function to implement mutations. Otherwise the mutations are
    applied to the whole population and this function isn't used.

\*****************************************************************************/
void Lemur::mutate(LRandom* generator)
{
    double expected_num_mutations = mutationRate * numGenes;
    int i_num_mutations = int(expected_num_mutations);

    /*
        If the mutation rate gets very low, the above calculation will always
        produce zero, but we know that we will still have a non-zero
        probability of mutation. Say the expected number of mutations is 0.1.
        What that means is that one of of ten creatures will have one mutation.
        We therefore use the fractional part to probabalistically increment the
        number of expected mutations.
     */
    double fractional_num_mutations = expected_num_mutations - i_num_mutations;

    if(generator->randomDouble() < fractional_num_mutations)
    {
        ++i_num_mutations;
    }


    for (int n = 0; n < i_num_mutations; n++)
    {
        unsigned which_gene = generator->randomInt(LRandom::GENE_LOCUS_RNDDIST);
        chromosome[which_gene] = generator->randomInt(LRandom::ALLELE_RNDDIST);
    }

    // Lastly, mutate the mutation rate itself!
    if (generator->randomDouble() < mutationRate)
    {
        mutationRate = generator->randomDouble();
    }

}


/******************************************************************************\

    Create a child by mixing the genes of this with the mate.

\******************************************************************************/
Lemur* Lemur::mateWith(const Lemur* mate, LRandom* generator) const
{
    // Doing reads via c_str() empirically has a significant performance advantage.
    const char* moms_chrom = chromosome.c_str();
    const char* dads_chrom = mate->chromosome.c_str();

    Lemur* child = acquire();
    std::string& childs_chrom = child->chromosome;

    /*
        Split at two points separated by half the length of the chromosome, so
        that we get half the genetic information from each parent. If we only
        split at one point we would not be getting equal numbers of genes from
        each parent.
     
        Randomizing the sections used means that a single couple can produce a
        large number of variations in their offspring (just like in real life).
    */
    int crossover1 = generator->randomInt(LRandom::GENE_LOCUS_RNDDIST) / 2;
    int crossover2 = crossover1 + numGenes / 2;
    Assert(crossover2 < numGenes + 1);

    for(int n = 0; n < crossover1; n++)
    {
        childs_chrom[n] = moms_chrom[n];
    }

    for(int n = crossover1; n < crossover2; n++)
    {
        childs_chrom[n] = dads_chrom[n];
    }

    for(int n = crossover2; n < numGenes; n++)
    {
        childs_chrom[n] = moms_chrom[n];
    }

    if(generator->randomBool())
    {
        child->mutationRate = this->mutationRate;
    }else
    {
        child->mutationRate = mate->mutationRate;
    }

    return child;
}
