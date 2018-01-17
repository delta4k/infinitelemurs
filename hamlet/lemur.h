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
 
    Lemurs are the creatures in our simulation. They carry the genetic material,
    and reproduce.
 
    Why Lemurs? Because "infinitemonkeys.com" was already taken. Besides,
    Lemurs are cuter than monkeys, and Lemurs are a better demonstration of
    the principles of evolution, having evolved in isolation on the island
    of Madagascar.

\*****************************************************************************/
class Lemur
{
private:

    Lemur* next; // Used for freelist
    Lemur();    // Only allow construction by acquire()
    ~Lemur() = delete;  // never delete, just go to freelist via release()


public:

    // Use a freelist rather than construction/destruction. Much faster.
    static Lemur* acquire();
    static void release(Lemur** lemur);

    void randomize(LRandom* generator);
    void mutate(LRandom* generator);

    Lemur* mateWith(const Lemur* mate, LRandom* generator) const;

    // Each character is a "gene". numGenes is the length of the text.
    static int numGenes;

    // See comments on mutation rates in SimParameters
    double mutationRate;

    // Number of characters that match the definitionOfFitness
    int fitness; 

    // The chromosome is the string of text that defines this creature. Random
    // in the first generation, produced by crossover with mutation for
    // subsequent generations
    std::string chromosome;
};

