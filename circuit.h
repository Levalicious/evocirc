#pragma once

#include "types.h"
#include "heap.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MINDEL (0.5f)
#define MAXDEL (4.f)

typedef struct {
    u64 energy;
    u64* code;
    u64* repcode;
    u64 clen;
    u64 hash;
    u64 defects;
    u64 zeros;
    u64 born;
} circ;

u64 ru(u64* state) {
	u64 *s = state;
    
	u64 const result = (((s[1] * 5) << 7) | ((s[1] * 5) >> (64 - 7))) * 9;
	u64 const t = s[1] << 17;

	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];

	s[2] ^= t;
	s[3] = (s[3] << 45) | (s[3] >> (64 - 45));

	return result;
}

f32 rf(u64 *s) {
    u32* seed = (u32*) s;
    f32 res;
    seed[0] = (i32) (((i32) seed[0]) * 16807);
    *((u32*) &res) = (((seed[0]) >> 9U) | 0x3f800000U);
    return res - 1.f;
}

void hashcirc(circ* c) {
    c->hash = 0;

    for (u64 i = 0; i < c->clen; ++i) {
        c->hash ^= 1LU << (c->code[i] % 61U);
        c->hash ^= (1LU << 61U) << (c->code[i] % 3U); 
    }
}

circ* initcirc(u64 len) {
    circ* out = (circ*) malloc(sizeof(circ));
    if (out == NULL) {
        printf("Failed to init circuit\n");
        exit(-1);
    }
    out->hash = 0;
    out->clen = len + 5;
    out->code = (u64*) calloc(out->clen, sizeof(u64));
    out->repcode = (u64*) calloc(out->clen, sizeof(u64));
    if (out->code == NULL || out->repcode == NULL) {
        printf("Failed to init circuit.\n");
        exit(-1);
    }
    out->defects = UINT64_MAX;
    out->energy = 50;
    out->zeros = 0;
    hashcirc(out);
    return out;
}

void randcirc(u64* state, circ* c) {
    for (u64 i = 0; i < c->clen; ++i) {
        c->code[i] = ru(state);
        /* Clear 0 bit */
        c->code[i] &= ~(1UL);

        /* 50% chance to toggle 0 bit */
        if (ru(state) % 2 == 0) {
            c->code[i] |= 1UL;
        }
        c->repcode[i] = ru(state);
    }
    c->zeros = 0;
    hashcirc(c);
}

void mutcirc(u64* state, circ* c, f32 tmut, f32 bmut) {
    for (u64 i = 0; i < c->clen; ++i) {
        f32 r = rf(state);
        if (r < tmut) {
            c->code[i] ^= 1UL;
        }
        if (r < bmut) {
            c->code[i] ^= (1UL << 1U) << (ru(state) % 63);
        }

        c->repcode[i] ^= (1LU << (ru(state) % 64));
    }

    c->zeros = 0;

    hashcirc(c);
}

void crosscirc(u64* state, circ* c, circ* a, circ* b) {
    u64 xover, afirst;
    afirst = ru(state) & 1LU;
    xover = ru(state) % a->clen;

    if (afirst) {
        memcpy(c->code, a->code, sizeof(u64) * xover);
        memcpy(c->code + xover, b->code + xover, (a->clen - xover) * sizeof(u64));
    } else {
        memcpy(c->code, b->code, sizeof(u64) * xover);
        memcpy(c->code + xover, a->code + xover, (a->clen - xover) * sizeof(u64));
    }
    c->zeros = 0;
}

void repcirc(circ* c, circ* a) {
    memcpy(c->code, a->code, sizeof(u64) * a->clen);
    memcpy(c->repcode, a->repcode, sizeof(u64) * a->clen);
    // for (u64 i = 0; i < a->clen; ++i) c->code[i] ^= c->repcode[i];
    c->zeros = 0;
}

f32 calcdel(u64* state, f32 mindel, f32 maxdel, f32 t, u32 tr) {
    return mindel + (rf(state) * (maxdel - mindel));
}

void test0(sigheap* h, f32* vins, circ* c, u64* seed) {
    u64 dmsk = ((1U << 31U) - 1U);

    f32 currt = 0.f;

    /* LO signal from A */
    u32 a1 = (c->code[0] >> 2U) & dmsk;
    u32 a2 = (c->code[0] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, a1), a1, 0.1f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, a2), a2, 0.1f);
    
    /* LO signal from B */
    u32 b1 = (c->code[1] >> 2U) & dmsk;
    u32 b2 = (c->code[1] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, b1), b1, 0.1f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, b2), b2, 0.1f);

    /* HI signal from t */
    u32 t1 = (c->code[2] >> 2U) & dmsk;
    u32 t2 = (c->code[2] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1), t1, 1.f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2), t2, 1.f);

    /* HI signal from P (power) */
    u32 P1 = (c->code[3] >> 2U) & dmsk;
    u32 P2 = (c->code[3] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, P1), P1, 1.f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, P2), P2, 1.f);

    u32 currind = 0;
    u32 tind = 0;
    f32 currv = 0.f;
    while (remmin(h, &currt, &currind, &currv) == 0) {
        if (c->energy > 0) {
            c->energy--;
        } else {
            break;
        }

        currind %= (c->clen * 2);
        tind = currind % c->clen;

        if (tind < 4) {
            /* Patch in */
            tind += 4;
        }
        if (tind < 6) {
            /* Output nodes. */
            /* Test 0: 0 0 1 -> 0 1 or x 0 */
            vins[tind] = currv;
            /* If 'complete' goes HI */
            if (vins[5] > 0.7f) {
                if (vins[4] > 0.7f) {
                    /* Illegal: 0 0 1 -> 1 1 */
                    c->defects++;
                }
            }

            if (vins[4] > 0.7f) {
                if (vins[5] > 0.7f) {
                    /* Illegal: 0 0 1 -> 1 1 */
                    c->defects++;
                }
            }

            continue;
        }

        vins[currind] = currv;

        t1 = (c->code[tind] >> 2U) & dmsk;
        t2 = (c->code[tind] >> 33U);

        if ((c->code[tind] & 0b10LU)) {
            /* Wire */
            insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind]);
            insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind]);
        } else {
            if ((c->code[tind] & 1LU)) {
                /* 1 : P-type */
                /* LO : Connected */
                /* TODO: Small loss 'a'=1, large loss 'a'=0 */
                if (vins[tind + c->clen] < 0.3f) {
                    /* LO : Connected */
                    /* Forward current value of 'a' to both recipients */
                    if (vins[tind] > 0.7f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.95);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.95);
                    } else if (vins[tind] < 0.3f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.75);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.75);
                    } else {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                    }
                } else if (vins[tind + c->clen] > 0.7f) {
                    /* HI : Disconnected */
                    /* Forward 0.1 to both recipients. LO, but circuit is active. */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                } else {
                    /* ∅ : ∅ */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                }
            } else {
                /* 0 : N-type */
                /* HI : Connected */
                /* TODO: Small loss 'a'=0, large loss 'a'=1 */
                if (vins[tind + c->clen] > 0.7f) {
                    /* HI : Connected */
                    /* Forward current value of 'a' to both recipients */
                    if (vins[tind] < 0.3f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.95);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.95);
                    } else if (vins[tind] > 0.7f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.75);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.75);
                    } else {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                    }
                } else if (vins[tind + c->clen] < 0.3f) {
                    /* LO : Disconnected */
                    /* Forward 0.1 to both recipients. LO, but circuit is active. */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                } else {
                    /* ∅ : ∅ */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                }
            }
        }
    }

    /* Must eventually 'complete' */
    /* 0 0 1 -> 0 1 */
    if (!(vins[5] > 0.7f)) c->defects++;
    memset(vins, 0, sizeof(f32) * c->clen * 2);
    h->n = 0;
}

void test1(sigheap* h, f32* vins, circ* c, u64* seed) {
    u64 dmsk = ((1U << 31U) - 1U);

    f32 currt = 0.f;

    /* LO signal from A */
    u32 a1 = (c->code[0] >> 2U) & dmsk;
    u32 a2 = (c->code[0] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, a1), a1, 0.1f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, a2), a2, 0.1f);
    
    /* HI signal from B */
    u32 b1 = (c->code[1] >> 2U) & dmsk;
    u32 b2 = (c->code[1] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, b1), b1, 1.f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, b2), b2, 1.f);

    /* HI signal from t */
    u32 t1 = (c->code[2] >> 2U) & dmsk;
    u32 t2 = (c->code[2] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1), t1, 1.f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2), t2, 1.f);

    /* HI signal from P (power) */
    u32 P1 = (c->code[3] >> 2U) & dmsk;
    u32 P2 = (c->code[3] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, P1), P1, 1.f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, P2), P2, 1.f);

    u32 currind = 0;
    u32 tind = 0;
    f32 currv = 0.f;
    while (remmin(h, &currt, &currind, &currv) == 0) {
        if (c->energy > 0) {
            c->energy--;
        } else {
            break;
        }

        currind %= (c->clen * 2);
        tind = currind % c->clen;

        if (tind < 4) {
            /* Patch in */
            tind += 4;
        }
        if (tind < 6) {
            /* Output nodes. */
            /* Test 1: 0 1 1 -> 0 1 or x 0 */
            vins[tind] = currv;
            /* If 'complete' goes HI */
            if (vins[5] > 0.7f) {
                if (vins[4] > 0.7f) {
                    /* Illegal: 0 1 1 -> 1 1 */
                    c->defects++;
                }
            }

            if (vins[4] > 0.7f) {
                if (vins[5] > 0.7f) {
                    /* Illegal: 0 1 1 -> 1 1 */
                    c->defects++;
                }
            }

            continue;
        }

        vins[currind] = currv;

        t1 = (c->code[tind] >> 2U) & dmsk;
        t2 = (c->code[tind] >> 33U);

        if ((c->code[tind] & 0b10LU)) {
            /* Wire */
            insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind]);
            insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind]);
        } else {
            if ((c->code[tind] & 1LU)) {
                /* 1 : P-type */
                /* LO : Connected */
                /* TODO: Small loss 'a'=1, large loss 'a'=0 */
                if (vins[tind + c->clen] < 0.3f) {
                    /* LO : Connected */
                    /* Forward current value of 'a' to both recipients */
                    if (vins[tind] > 0.7f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.95);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.95);
                    } else if (vins[tind] < 0.3f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.75);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.75);
                    } else {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                    }
                } else if (vins[tind + c->clen] > 0.7f) {
                    /* HI : Disconnected */
                    /* Forward 0.1 to both recipients. LO, but circuit is active. */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                } else {
                    /* ∅ : ∅ */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                }
            } else {
                /* 0 : N-type */
                /* HI : Connected */
                /* TODO: Small loss 'a'=0, large loss 'a'=1 */
                if (vins[tind + c->clen] > 0.7f) {
                    /* HI : Connected */
                    /* Forward current value of 'a' to both recipients */
                    if (vins[tind] < 0.3f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.95);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.95);
                    } else if (vins[tind] > 0.7f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.75);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.75);
                    } else {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                    }
                } else if (vins[tind + c->clen] < 0.3f) {
                    /* LO : Disconnected */
                    /* Forward 0.1 to both recipients. LO, but circuit is active. */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                } else {
                    /* ∅ : ∅ */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                }
            }
        }
    }

    /* Must eventually 'complete' */
    /* 0 1 1 -> 0 1 */
    if (!(vins[5] > 0.7f)) c->defects++;
    memset(vins, 0, sizeof(f32) * c->clen * 2);
    h->n = 0;
}

void test2(sigheap* h, f32* vins, circ* c, u64* seed) {
    u64 dmsk = ((1U << 31U) - 1U);

    f32 currt = 0.f;

    /* HI signal from A */
    u32 a1 = (c->code[0] >> 2U) & dmsk;
    u32 a2 = (c->code[0] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, a1), a1, 1.f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, a2), a2, 1.f);
    
    /* LO signal from B */
    u32 b1 = (c->code[1] >> 2U) & dmsk;
    u32 b2 = (c->code[1] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, b1), b1, 0.1f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, b2), b2, 0.1f);

    /* HI signal from t */
    u32 t1 = (c->code[2] >> 2U) & dmsk;
    u32 t2 = (c->code[2] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1), t1, 1.f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2), t2, 1.f);

    /* HI signal from P (power) */
    u32 P1 = (c->code[3] >> 2U) & dmsk;
    u32 P2 = (c->code[3] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, P1), P1, 1.f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, P2), P2, 1.f);

    u32 currind = 0;
    u32 tind = 0;
    f32 currv = 0.f;
    while (remmin(h, &currt, &currind, &currv) == 0) {
        if (c->energy > 0) {
            c->energy--;
        } else {
            break;
        }

        currind %= (c->clen * 2);
        tind = currind % c->clen;

        if (tind < 4) {
            /* Patch in */
            tind += 4;
        }
        if (tind < 6) {
            /* Output nodes. */
            /* Test 2: 1 0 1 -> 0 1 or x 0 */
            vins[tind] = currv;
            /* If 'complete' goes HI */
            if (vins[5] > 0.7f) {
                if (vins[4] > 0.7f) {
                    /* Illegal: 1 0 1 -> 1 1 */
                    c->defects++;
                }
            }

            if (vins[4] > 0.7f) {
                if (vins[5] > 0.7f) {
                    /* Illegal: 1 0 1 -> 1 1 */
                    c->defects++;
                }
            }

            continue;
        }

        vins[currind] = currv;

        t1 = (c->code[tind] >> 2U) & dmsk;
        t2 = (c->code[tind] >> 33U);

        if ((c->code[tind] & 0b10LU)) {
            /* Wire */
            insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind]);
            insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind]);
        } else {
            if ((c->code[tind] & 1LU)) {
                /* 1 : P-type */
                /* LO : Connected */
                /* TODO: Small loss 'a'=1, large loss 'a'=0 */
                if (vins[tind + c->clen] < 0.3f) {
                    /* LO : Connected */
                    /* Forward current value of 'a' to both recipients */
                    if (vins[tind] > 0.7f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.95);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.95);
                    } else if (vins[tind] < 0.3f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.75);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.75);
                    } else {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                    }
                } else if (vins[tind + c->clen] > 0.7f) {
                    /* HI : Disconnected */
                    /* Forward 0.1 to both recipients. LO, but circuit is active. */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                } else {
                    /* ∅ : ∅ */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                }
            } else {
                /* 0 : N-type */
                /* HI : Connected */
                /* TODO: Small loss 'a'=0, large loss 'a'=1 */
                if (vins[tind + c->clen] > 0.7f) {
                    /* HI : Connected */
                    /* Forward current value of 'a' to both recipients */
                    if (vins[tind] < 0.3f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.95);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.95);
                    } else if (vins[tind] > 0.7f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.75);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.75);
                    } else {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                    }
                } else if (vins[tind + c->clen] < 0.3f) {
                    /* LO : Disconnected */
                    /* Forward 0.1 to both recipients. LO, but circuit is active. */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                } else {
                    /* ∅ : ∅ */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                }
            }
        }
    }

    /* Must eventually 'complete' */
    /* 1 0 1 -> 0 1 */
    if (!(vins[5] > 0.7f)) c->defects++;
    memset(vins, 0, sizeof(f32) * c->clen * 2);
    h->n = 0;
}

void test3(sigheap* h, f32* vins, circ* c, u64* seed) {
    u64 dmsk = ((1U << 31U) - 1U);

    f32 currt = 0.f;

    /* HI signal from A */
    u32 a1 = (c->code[0] >> 2U) & dmsk;
    u32 a2 = (c->code[0] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, a1), a1, 1.f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, a2), a2, 1.f);
    
    /* HI signal from B */
    u32 b1 = (c->code[1] >> 2U) & dmsk;
    u32 b2 = (c->code[1] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, b1), b1, 1.f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, b2), b2, 1.f);

    /* HI signal from t */
    u32 t1 = (c->code[2] >> 2U) & dmsk;
    u32 t2 = (c->code[2] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1), t1, 1.f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2), t2, 1.f);

    /* HI signal from P (power) */
    u32 P1 = (c->code[3] >> 2U) & dmsk;
    u32 P2 = (c->code[3] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, P1), P1, 1.f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, P2), P2, 1.f);

    u32 currind = 0;
    u32 tind = 0;
    f32 currv = 0.f;
    while (remmin(h, &currt, &currind, &currv) == 0) {
        if (c->energy > 0) {
            c->energy--;
        } else {
            break;
        }

        currind %= (c->clen * 2);
        tind = currind % c->clen;

        if (tind < 4) {
            /* Patch in */
            tind += 4;
        }
        if (tind < 6) {
            /* Output nodes. */
            /* Test 3: 1 1 1 -> 1 1 */
            vins[tind] = currv;
            /* If 'complete' goes HI */
            if (vins[5] > 0.7f) {
                if (!(vins[4] > 0.7f)) {
                    /* Illegal: 
                     * 1 1 1 -> 0 1 
                     * 1 1 1 -> ∅ 1 */
                    c->defects++;
                }
            }

            if (!(vins[4] > 0.7f)) {
                if (vins[5] > 0.7f) {
                    /* Illegal: 
                     * 1 1 1 -> 0 1 
                     * 1 1 1 -> ∅ 1 */
                    c->defects++;
                }
            }

            continue;
        }

        vins[currind] = currv;

        t1 = (c->code[tind] >> 2U) & dmsk;
        t2 = (c->code[tind] >> 33U);

        if ((c->code[tind] & 0b10LU)) {
            /* Wire */
            insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind]);
            insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind]);
        } else {
            if ((c->code[tind] & 1LU)) {
                /* 1 : P-type */
                /* LO : Connected */
                /* TODO: Small loss 'a'=1, large loss 'a'=0 */
                if (vins[tind + c->clen] < 0.3f) {
                    /* LO : Connected */
                    /* Forward current value of 'a' to both recipients */
                    if (vins[tind] > 0.7f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.95);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.95);
                    } else if (vins[tind] < 0.3f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.75);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.75);
                    } else {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                    }
                } else if (vins[tind + c->clen] > 0.7f) {
                    /* HI : Disconnected */
                    /* Forward 0.1 to both recipients. LO, but circuit is active. */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                } else {
                    /* ∅ : ∅ */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                }
            } else {
                /* 0 : N-type */
                /* HI : Connected */
                /* TODO: Small loss 'a'=0, large loss 'a'=1 */
                if (vins[tind + c->clen] > 0.7f) {
                    /* HI : Connected */
                    /* Forward current value of 'a' to both recipients */
                    if (vins[tind] < 0.3f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.95);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.95);
                    } else if (vins[tind] > 0.7f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.75);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.75);
                    } else {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                    }
                } else if (vins[tind + c->clen] < 0.3f) {
                    /* LO : Disconnected */
                    /* Forward 0.1 to both recipients. LO, but circuit is active. */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                } else {
                    /* ∅ : ∅ */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                }
            }
        }
    }

    /* Must eventually 'complete' */
    /* 1 1 1 -> 1 1 */
    if (!(vins[5] > 0.7f)) c->defects++;
    memset(vins, 0, sizeof(f32) * c->clen * 2);
    h->n = 0;
}

void test4(sigheap* h, f32* vins, circ* c, u64* seed) {
    u64 dmsk = ((1U << 31U) - 1U);

    f32 currt = 0.f;

    /* LO signal from A */
    u32 a1 = (c->code[0] >> 2U) & dmsk;
    u32 a2 = (c->code[0] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, a1), a1, 0.1f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, a2), a2, 0.1f);
    
    /* LO signal from B */
    u32 b1 = (c->code[1] >> 2U) & dmsk;
    u32 b2 = (c->code[1] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, b1), b1, 0.1f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, b2), b2, 0.1f);

    /* LO signal from t */
    u32 t1 = (c->code[2] >> 2U) & dmsk;
    u32 t2 = (c->code[2] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1), t1, 0.1f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2), t2, 0.1f);

    /* HI signal from P (power) */
    u32 P1 = (c->code[3] >> 2U) & dmsk;
    u32 P2 = (c->code[3] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, P1), P1, 1.f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, P2), P2, 1.f);

    u32 currind = 0;
    u32 tind = 0;
    f32 currv = 0.f;
    while (remmin(h, &currt, &currind, &currv) == 0) {
        if (c->energy > 0) {
            c->energy--;
        } else {
            break;
        }
        
        currind %= (c->clen * 2);
        tind = currind % c->clen;

        if (tind < 4) {
            /* Patch in */
            tind += 4;
        }
        if (tind < 6) {
            /* Output nodes. */
            /* Test 4: 0 0 0 -> x 0 */
            vins[tind] = currv;
            /* If 'complete' goes HI */
            if (vins[5] > 0.7f) {
                c->defects++;
            }

            continue;
        }

        vins[currind] = currv;

        t1 = (c->code[tind] >> 2U) & dmsk;
        t2 = (c->code[tind] >> 33U);

        if ((c->code[tind] & 0b10LU)) {
            /* Wire */
            insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind]);
            insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind]);
        } else {
            if ((c->code[tind] & 1LU)) {
                /* 1 : P-type */
                /* LO : Connected */
                /* TODO: Small loss 'a'=1, large loss 'a'=0 */
                if (vins[tind + c->clen] < 0.3f) {
                    /* LO : Connected */
                    /* Forward current value of 'a' to both recipients */
                    if (vins[tind] > 0.7f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.95);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.95);
                    } else if (vins[tind] < 0.3f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.75);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.75);
                    } else {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                    }
                } else if (vins[tind + c->clen] > 0.7f) {
                    /* HI : Disconnected */
                    /* Forward 0.1 to both recipients. LO, but circuit is active. */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                } else {
                    /* ∅ : ∅ */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                }
            } else {
                /* 0 : N-type */
                /* HI : Connected */
                /* TODO: Small loss 'a'=0, large loss 'a'=1 */
                if (vins[tind + c->clen] > 0.7f) {
                    /* HI : Connected */
                    /* Forward current value of 'a' to both recipients */
                    if (vins[tind] < 0.3f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.95);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.95);
                    } else if (vins[tind] > 0.7f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.75);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.75);
                    } else {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                    }
                } else if (vins[tind + c->clen] < 0.3f) {
                    /* LO : Disconnected */
                    /* Forward 0.1 to both recipients. LO, but circuit is active. */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                } else {
                    /* ∅ : ∅ */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                }
            }
        }
    }

    memset(vins, 0, sizeof(f32) * c->clen * 2);
    h->n = 0;
}

void test5(sigheap* h, f32* vins, circ* c, u64* seed) {
    u64 dmsk = ((1U << 31U) - 1U);

    f32 currt = 0.f;

    /* LO signal from A */
    u32 a1 = (c->code[0] >> 2U) & dmsk;
    u32 a2 = (c->code[0] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, a1), a1, 0.1f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, a2), a2, 0.1f);
    
    /* HI signal from B */
    u32 b1 = (c->code[1] >> 2U) & dmsk;
    u32 b2 = (c->code[1] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, b1), b1, 1.f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, b2), b2, 1.f);

    /* LO signal from t */
    u32 t1 = (c->code[2] >> 2U) & dmsk;
    u32 t2 = (c->code[2] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1), t1, 0.1f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2), t2, 0.1f);

    /* HI signal from P (power) */
    u32 P1 = (c->code[3] >> 2U) & dmsk;
    u32 P2 = (c->code[3] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, P1), P1, 1.f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, P2), P2, 1.f);

    u32 currind = 0;
    u32 tind = 0;
    f32 currv = 0.f;
    while (remmin(h, &currt, &currind, &currv) == 0) {
        if (c->energy > 0) {
            c->energy--;
        } else {
            break;
        }

        currind %= (c->clen * 2);
        tind = currind % c->clen;

        if (tind < 4) {
            /* Patch in */
            tind += 4;
        }
        if (tind < 6) {
            /* Output nodes. */
            /* Test 5: 0 1 0 -> x 0 */
            vins[tind] = currv;
            /* If 'complete' goes HI */
            if (vins[5] > 0.7f) {
                c->defects++;
            }

            continue;
        }

        vins[currind] = currv;

        t1 = (c->code[tind] >> 2U) & dmsk;
        t2 = (c->code[tind] >> 33U);

        if ((c->code[tind] & 0b10LU)) {
            /* Wire */
            insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind]);
            insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind]);
        } else {
            if ((c->code[tind] & 1LU)) {
                /* 1 : P-type */
                /* LO : Connected */
                /* TODO: Small loss 'a'=1, large loss 'a'=0 */
                if (vins[tind + c->clen] < 0.3f) {
                    /* LO : Connected */
                    /* Forward current value of 'a' to both recipients */
                    if (vins[tind] > 0.7f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.95);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.95);
                    } else if (vins[tind] < 0.3f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.75);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.75);
                    } else {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                    }
                } else if (vins[tind + c->clen] > 0.7f) {
                    /* HI : Disconnected */
                    /* Forward 0.1 to both recipients. LO, but circuit is active. */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                } else {
                    /* ∅ : ∅ */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                }
            } else {
                /* 0 : N-type */
                /* HI : Connected */
                /* TODO: Small loss 'a'=0, large loss 'a'=1 */
                if (vins[tind + c->clen] > 0.7f) {
                    /* HI : Connected */
                    /* Forward current value of 'a' to both recipients */
                    if (vins[tind] < 0.3f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.95);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.95);
                    } else if (vins[tind] > 0.7f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.75);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.75);
                    } else {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                    }
                } else if (vins[tind + c->clen] < 0.3f) {
                    /* LO : Disconnected */
                    /* Forward 0.1 to both recipients. LO, but circuit is active. */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                } else {
                    /* ∅ : ∅ */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                }
            }
        }
    }

    memset(vins, 0, sizeof(f32) * c->clen * 2);
    h->n = 0;
}

void test6(sigheap* h, f32* vins, circ* c, u64* seed) {
    u64 dmsk = ((1U << 31U) - 1U);

    f32 currt = 0.f;

    /* HI signal from A */
    u32 a1 = (c->code[0] >> 2U) & dmsk;
    u32 a2 = (c->code[0] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, a1), a1, 1.f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, a2), a2, 1.f);
    
    /* LO signal from B */
    u32 b1 = (c->code[1] >> 2U) & dmsk;
    u32 b2 = (c->code[1] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, b1), b1, 0.1f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, b2), b2, 0.1f);

    /* LO signal from t */
    u32 t1 = (c->code[2] >> 2U) & dmsk;
    u32 t2 = (c->code[2] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1), t1, 0.1f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2), t2, 0.1f);

    /* HI signal from P (power) */
    u32 P1 = (c->code[3] >> 2U) & dmsk;
    u32 P2 = (c->code[3] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, P1), P1, 1.f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, P2), P2, 1.f);

    u32 currind = 0;
    u32 tind = 0;
    f32 currv = 0.f;
    while (remmin(h, &currt, &currind, &currv) == 0) {
        if (c->energy > 0) {
            c->energy--;
        } else {
            break;
        }

        currind %= (c->clen * 2);
        tind = currind % c->clen;

        if (tind < 4) {
            /* Patch in */
            tind += 4;
        }
        if (tind < 6) {
            /* Output nodes. */
            /* Test 6: 1 0 0 -> x 0 */
            vins[tind] = currv;
            /* If 'complete' goes HI */
            if (vins[5] > 0.7f) {
                c->defects++;
            }

            continue;
        }

        vins[currind] = currv;

        t1 = (c->code[tind] >> 2U) & dmsk;
        t2 = (c->code[tind] >> 33U);

        if ((c->code[tind] & 0b10LU)) {
            /* Wire */
            insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind]);
            insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind]);
        } else {
            if ((c->code[tind] & 1LU)) {
                /* 1 : P-type */
                /* LO : Connected */
                /* TODO: Small loss 'a'=1, large loss 'a'=0 */
                if (vins[tind + c->clen] < 0.3f) {
                    /* LO : Connected */
                    /* Forward current value of 'a' to both recipients */
                    if (vins[tind] > 0.7f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.95);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.95);
                    } else if (vins[tind] < 0.3f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.75);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.75);
                    } else {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                    }
                } else if (vins[tind + c->clen] > 0.7f) {
                    /* HI : Disconnected */
                    /* Forward 0.1 to both recipients. LO, but circuit is active. */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                } else {
                    /* ∅ : ∅ */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                }
            } else {
                /* 0 : N-type */
                /* HI : Connected */
                /* TODO: Small loss 'a'=0, large loss 'a'=1 */
                if (vins[tind + c->clen] > 0.7f) {
                    /* HI : Connected */
                    /* Forward current value of 'a' to both recipients */
                    if (vins[tind] < 0.3f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.95);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.95);
                    } else if (vins[tind] > 0.7f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.75);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.75);
                    } else {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                    }
                } else if (vins[tind + c->clen] < 0.3f) {
                    /* LO : Disconnected */
                    /* Forward 0.1 to both recipients. LO, but circuit is active. */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                } else {
                    /* ∅ : ∅ */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                }
            }
        }
    }

    memset(vins, 0, sizeof(f32) * c->clen * 2);
    h->n = 0;
}

void test7(sigheap* h, f32* vins, circ* c, u64* seed) {
    u64 dmsk = ((1U << 31U) - 1U);

    f32 currt = 0.f;

    /* HI signal from A */
    u32 a1 = (c->code[0] >> 2U) & dmsk;
    u32 a2 = (c->code[0] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, a1), a1, 1.f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, a2), a2, 1.f);
    
    /* HI signal from B */
    u32 b1 = (c->code[1] >> 2U) & dmsk;
    u32 b2 = (c->code[1] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, b1), b1, 1.f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, b2), b2, 1.f);

    /* LO signal from t */
    u32 t1 = (c->code[2] >> 2U) & dmsk;
    u32 t2 = (c->code[2] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1), t1, 0.1f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2), t2, 0.1f);

    /* HI signal from P (power) */
    u32 P1 = (c->code[3] >> 2U) & dmsk;
    u32 P2 = (c->code[3] >> 33U);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, P1), P1, 1.f);
    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, P2), P2, 1.f);

    u32 currind = 0;
    u32 tind = 0;
    f32 currv = 0.f;
    while (remmin(h, &currt, &currind, &currv) == 0) {
        if (c->energy > 0) {
            c->energy--;
        } else {
            break;
        }
        
        currind %= (c->clen * 2);
        tind = currind % c->clen;

        if (tind < 4) {
            /* Patch in */
            tind += 4;
        }
        if (tind < 6) {
            /* Output nodes. */
            /* Test 7: 1 1 0 -> x 0 */
            vins[tind] = currv;
            /* If 'complete' goes HI */
            if (vins[5] > 0.7f) {
                c->defects++;
            }

            continue;
        }

        vins[currind] = currv;

        t1 = (c->code[tind] >> 2U) & dmsk;
        t2 = (c->code[tind] >> 33U);

        if ((c->code[tind] & 0b10LU)) {
            /* Wire */
            insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind]);
            insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind]);
        } else {
            if ((c->code[tind] & 1LU)) {
                /* 1 : P-type */
                /* LO : Connected */
                /* TODO: Small loss 'a'=1, large loss 'a'=0 */
                if (vins[tind + c->clen] < 0.3f) {
                    /* LO : Connected */
                    /* Forward current value of 'a' to both recipients */
                    if (vins[tind] > 0.7f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.95);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.95);
                    } else if (vins[tind] < 0.3f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.75);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.75);
                    } else {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                    }
                } else if (vins[tind + c->clen] > 0.7f) {
                    /* HI : Disconnected */
                    /* Forward 0.1 to both recipients. LO, but circuit is active. */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                } else {
                    /* ∅ : ∅ */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                }
            } else {
                /* 0 : N-type */
                /* HI : Connected */
                /* TODO: Small loss 'a'=0, large loss 'a'=1 */
                if (vins[tind + c->clen] > 0.7f) {
                    /* HI : Connected */
                    /* Forward current value of 'a' to both recipients */
                    if (vins[tind] < 0.3f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.95);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.95);
                    } else if (vins[tind] > 0.7f) {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, vins[tind] * 0.75);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, vins[tind] * 0.75);
                    } else {
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                        insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                    }
                } else if (vins[tind + c->clen] < 0.3f) {
                    /* LO : Disconnected */
                    /* Forward 0.1 to both recipients. LO, but circuit is active. */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                } else {
                    /* ∅ : ∅ */
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t1 ^ currind), t1, 0.1f);
                    insmin(h, currt + calcdel(seed, MINDEL, MAXDEL, currt, t2 ^ currind), t2, 0.1f);
                }
            }
        }
    }

    memset(vins, 0, sizeof(f32) * c->clen * 2);
    h->n = 0;
}

#define TESTREPS (512)

void printcircuit(circ* c) {
    for (u64 i = 0; i < c->clen; ++i) {
        u32 o1, o2;
        char* o1s, *o2s;
        o1 = (c->code[i] >> 2U) & ((1U << 31U) - 1U);
        o2 = (c->code[i] >> 33U);
        o1s = (o1 < c->clen) ? ("a") : ("s");
        o2s = (o2 < c->clen) ? ("a") : ("s");
        o1 %= c->clen;
        o2 %= c->clen;
        if (c->code[i] & 0b10LU) {
            printf("W : %lu -> %s.%u %s.%u\n", i, o1s, o1, o2s, o2);
        } else {
            if (c->code[i] & 1LU) {
                printf("P : %lu -> %s.%u %s.%u\n", i, o1s, o1, o2s, o2);
            } else {
                printf("N : %lu -> %s.%u %s.%u\n", i, o1s, o1, o2s, o2);
            }
        }
        
    }
}

int run(sigheap* h, u64* seednoise, circ* c, f32* vins) {
    /* Lo: 0.0 - 0.3
     *  ∅: 0.3 - 0.7
     * Hi: 0.7 - 1.0 */
    /* Async AND gate
     * 0: 0 0 1 -> 0 1 or x 0
     * 1: 0 1 1 -> 0 1 or x 0
     * 2: 1 0 1 -> 0 1 or x 0
     * 3: 1 1 1 -> 1 1 or x 0
     * 4: 0 0 0 -> x 0
     * 5: 0 1 0 -> x 0
     * 6: 1 0 0 -> x 0
     * 7: 1 1 0 -> x 0 */

    c->defects = 0;
    u64 bnrg, enrg;

    bnrg = c->energy;
    test0(h, vins, c, seednoise);
    test1(h, vins, c, seednoise);
    test2(h, vins, c, seednoise);
    test3(h, vins, c, seednoise);
    test4(h, vins, c, seednoise);
    test5(h, vins, c, seednoise);
    test6(h, vins, c, seednoise);
    test7(h, vins, c, seednoise);
    enrg = c->energy;
    for (u64 i = 0; i < TESTREPS - 1; ++i) {
        test0(h, vins, c, seednoise);
        test1(h, vins, c, seednoise);
        test2(h, vins, c, seednoise);
        test3(h, vins, c, seednoise);
        test4(h, vins, c, seednoise);
        test5(h, vins, c, seednoise);
        test6(h, vins, c, seednoise);
        test7(h, vins, c, seednoise);
    }
    c->energy = enrg;
    if (bnrg == enrg) {
        bnrg = c->energy;
        enrg = 0;
        c->defects += c->energy;
        c->energy = 0;
    }

    if (c->defects == 0) {
        c->zeros++;
    } else {
        c->zeros = 0;
    }
    return bnrg - enrg;
}