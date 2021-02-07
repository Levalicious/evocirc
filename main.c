#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <signal.h>

#include "circuit.h"

#define POP (4096)
#define REPCST (1000)
#define REPWHEN (1500)
#define CIRCLN (40)

/* Choose an initial energy such that newborns die off
 * immediately if they don't provide an improvement. If they
 * do manage to be competitive, they may survive, and might also reproduce */
#define INITENERG (800)

#define COMPTHRESH (512)

#define REW (800)
#define FEEDINGS (1024)
#define FDRAT (2)
#define FDMX (2)
#define FDSIZE (3)

#define LMUR (0.0f)
#define MUR (0.3f)
#define TMUT (0.15f)
#define BMUT (0.3f)

#define REPRFREQ (0.0f)
#define REPRSIZE (10)

#define MAXITERS (5000000)

static volatile int keepRunning = 1;

void inthandler(int dummy) {
    keepRunning = 0;
}

void seedr(u64* state, u64 seed) {
    state[0] = 0xb6d47cfacccc53f8LU ^ seed;
    state[1] = 0x30b319a052624be7LU ^ seed;
    state[2] = 0xfbeb173c6d0227d8LU ^ seed;
    state[3] = 0x99cfe60a00bdd4feLU ^ seed;
}

int main() {
    u64 rstate[4];
    seedr(rstate, time(NULL));

    circ** pop = (circ**) malloc(sizeof(circ*) * POP);
    if (pop == NULL) {
        printf("Failed to allocate population array.\n");
        return 0;
    }

    circ** live = (circ**) malloc(sizeof(circ*) * POP);
    if (live == NULL) {
        printf("Failed to allocate living array.\n");
        return 0;
    }

    for (u64 i = 0; i < POP; ++i) {
        pop[i] = initcirc(CIRCLN);
        randcirc(rstate, pop[i]);
        pop[i]->energy = INITENERG;
        pop[i]->born = 0;
    }

    f32* vins = (f32*) calloc((CIRCLN + 5) * 2, sizeof(f32));
    if (vins == NULL) {
        printf("Failed to alloc voltage array.\n");
        return -1;
    }

    sigheap* h = initheap();

    u64 iters = 0;

    u64 alive = 0;
    u64 dead = 0;
    u64 borna = 0;
    u64 best = UINT64_MAX;
    u64 worst = 0;
    f64 avg = 0.0;
    f64 avglvng = 0.0;
    f64 avgenrg = 0.0;
    f64 avglvngenerg = 0.0;
    u64 generated = 0;
    u64 mutants = 0;
    u64 oldest, youngest;
    u64 maxzers = 0;
    u64 predead = 0;
    u64 borns = 0;
    f64 rncst = 0.0;

    while (keepRunning) {
        f32 iternoise = rf(rstate);

        /* Evaluate all circuits */
        alive = 0;
        dead = 0;
        borna = 0;
        borns = 0;
        best = UINT64_MAX;
        worst = 0;
        avg = 0.0;
        avglvng = 0.0;
        avgenrg = 0.0;
        avglvngenerg = 0.0;
        generated = 0;
        mutants = 0;
        predead = 0;
        rncst = 0.0;

        oldest = iters;
        youngest = 0;
        maxzers = 0;

        for (u64 currCirc = 0; currCirc < POP; ++currCirc) {
            u64 tstate[4];
            memcpy(tstate, rstate, sizeof(u64) * 4);
            // seedr(tstate, 5);
            /* Same noise for each circ */
            if (pop[currCirc]->energy == 0) predead++;
            int res = run(h, tstate, pop[currCirc], vins);
            rncst += (res / ((f64) POP));
            if (iters - pop[currCirc]->born > 100) {
                /* 'Old Age' */
                pop[currCirc]->energy /= (iters - pop[currCirc]->born) - 100;
            }
            if (pop[currCirc]->defects < best) best = pop[currCirc]->defects;
            if (pop[currCirc]->defects > worst) worst = pop[currCirc]->defects;

            if (pop[currCirc]->energy == 0) {
                dead++;
            } else {
                live[alive] = pop[currCirc];
                alive++;
                avglvng += ((f64) pop[currCirc]->defects);
                avglvngenerg += ((f64) pop[currCirc]->energy);
                if (pop[currCirc]->born > youngest) youngest = pop[currCirc]->born;
                if (pop[currCirc]->born < oldest) oldest = pop[currCirc]->born;
                if (pop[currCirc]->zeros > maxzers) maxzers = pop[currCirc]->zeros;
            }

            if (pop[currCirc]->zeros == COMPTHRESH) {
                printf("Found solution on iter %lu\n", iters);
                printcircuit(pop[currCirc]);
                return 0;
            }
        }

        if (alive) {
            /* Competition over food */
            for (u64 fdng = 0; fdng < FEEDINGS; ++fdng) {
                circ* fdgroup[FDSIZE];
                /* Pick random sample group */
                for (u64 i = 0; i < FDSIZE; ++i) {
                    fdgroup[i] = live[ru(rstate) % alive];
                }

                /* Bubble sort according to score */
                u8 swapped = 0;
                u64 sortn = FDSIZE;
                do {
                    swapped = 0;
                    for (u64 i = 1; i < sortn; ++i) {
                        if (fdgroup[i - 1]->defects > fdgroup[i]->defects) {
                            circ* temp = fdgroup[i];
                            fdgroup[i] = fdgroup[i - 1];
                            fdgroup[i - 1] = temp;
                            swapped = 1;
                        }
                    }
                    sortn--;
                } while (swapped);
                /* Reward based on configured distribution */
                u64 curew = REW;
                u64 curewind = 0;
                while (curew != 0 && curewind != FDMX) {
                    fdgroup[curewind]->energy += curew;
                    curew /= FDRAT;
                    curewind++;
                }
            }

            /* Competition over reproduction */
            for (u64 currCirc = 0; currCirc < POP; ++currCirc) {
                /* If not dead, don't try to replace. Maybe mutate. */
                if (pop[currCirc]->energy != 0) {
                    if (rf(rstate) < LMUR) {
                        mutcirc(rstate, pop[currCirc], TMUT, BMUT);
                        mutants++;
                        pop[currCirc]->born = iters;
                    }
                    continue;
                }

                circ* reprgroup[REPRSIZE];
                /* Pick random sample group */
                for (u64 i = 0; i < REPRSIZE; ++i) {
                    reprgroup[i] = live[ru(rstate) % alive];
                }

                /* Bubble sort according to energy */
                u8 swapped = 0;
                u64 sortn = FDSIZE;
                do {
                    swapped = 0;
                    for (u64 i = 1; i < sortn; ++i) {
                        if (reprgroup[i - 1]->defects > reprgroup[i]->defects) {
                            circ* temp = reprgroup[i];
                            reprgroup[i] = reprgroup[i - 1];
                            reprgroup[i - 1] = temp;
                            swapped = 1;
                        }
                    }
                    sortn--;
                } while (swapped);

                /* If top two have sufficient energy to reproduce,
                 * breed them and remove energy.
                 * TODO: If sexes are implemented, different energy costs. */
                if (reprgroup[0]->energy >= REPWHEN) { //  && rf(rstate) < REPRFREQ
                    if (rf(rstate) < REPRFREQ && reprgroup[1]->energy >= REPWHEN) { 
                        crosscirc(rstate, pop[currCirc], reprgroup[0], reprgroup[1]);
                        reprgroup[0]->energy -= REPCST;
                        reprgroup[1]->energy -= REPCST;
                        borns++;
                    } else {
                        repcirc(pop[currCirc], reprgroup[0]);
                        pop[currCirc]->energy = INITENERG;
                        reprgroup[0]->energy -= REPCST;
                        borna++;
                    }
                    
                    if (rf(rstate) < MUR) {
                        mutcirc(rstate, pop[currCirc], TMUT, BMUT);
                        mutants++;
                    }
                    pop[currCirc]->born = iters;
                } else {
                    /*
                    if (alive < POP / 3) {
                        randcirc(rstate, pop[currCirc]);
                        pop[currCirc]->energy = INITENERG;
                        pop[currCirc]->born = iters;
                        generated++;
                    }
                    */
                }
            }
        } else {
            printf("Regenerated solution pool.\n");
            for (u64 i = 0; i < POP; ++i) {
                /* If can't reproduce, generate new */
                randcirc(rstate, pop[i]);
                pop[i]->energy = INITENERG;
                pop[i]->born = iters;
                generated++;
            }
        }

        printf("Iteration %8lu : Pop. %lu , %lu deaths, %lu asexual births, %lu sexual births, %lu mutations\n", iters, alive, dead - predead, borna, borns, mutants);
        printf("\tBest circuit: %f\n", best / ((f64) TESTREPS));
        printf("\tRuncost: %f\n", rncst);
        if (alive) {
            printf("\tAvg living circuit: %f\n", ((avglvng) / ((f64) alive)) / ((f64) TESTREPS));
            printf("\tAvg living energy: %f\n", avglvngenerg / ((f64) alive));
            printf("\tAges: %lu to %lu\n", iters - youngest, iters - oldest);
            printf("\tMax zeros: %lu\n", maxzers);
        }
        
        iters++;
        if (iters == MAXITERS) break;
   }

    for (u64 i = 0; i < POP; ++i) {
        free(pop[i]->code);
        free(pop[i]->repcode);
        free(pop[i]);
    }
    freeheap(h);
    free(pop);
    free(live);
    free(vins);
    return 0;
}