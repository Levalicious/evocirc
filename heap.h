#pragma once

#include "types.h"
#include <stdlib.h>
#include <stdio.h>

typedef struct {
    f32* times;
    u32* inds;
    f32* vs;
    u64 n;
    u64 cap;
} sigheap;

#define D (4)

sigheap* initheap() {
    sigheap* out = (sigheap*) malloc(sizeof(sigheap));
    out->n = 0;
    out->cap = 32;
    out->times = (f32*) malloc(sizeof(f32) * out->cap);
    out->inds = (u32*) malloc(sizeof(u32) * out->cap);
    out->vs = (f32*) malloc(sizeof(f32) * out->cap);

    return out;
}

void freeheap(sigheap* heap) {
    free(heap->times);
    free(heap->inds);
    free(heap->vs);
    free(heap);
}

void hfymax(sigheap* heap, u64 i) {
    u64 min = i;
    for (u64 l = (i << D) + 1; l < (i << D) + D + 1; ++l) {
        if (l < heap->n && heap->times[l] > heap->times[min]) min = l;
    }

    if (min != i) {
        f32 tt = heap->times[0];
        u32 ti = heap->inds[0];
        f32 tv = heap->vs[0];
        heap->times[0] = heap->times[min];
        heap->inds[0] = heap->inds[min];
        heap->vs[0] = heap->vs[min];
        heap->times[min] = tt;
        heap->inds[min] = ti;
        heap->vs[min] = tv;

        hfymax(heap, min);
    }
}

void insmax(sigheap* heap, f32 it, u32 iind, f32 v) {
    if (heap->n == heap->cap) {
        heap->cap *= 2;
        heap->times = (f32*) realloc(heap->times, sizeof(f32) * heap->cap);
        if (heap->times == NULL) {
            printf("Failed to realloc heap times.\n");
            exit(0);
        }

        heap->inds = (u32*) realloc(heap->inds, sizeof(u32) * heap->cap);
        if (heap->inds == NULL) {
            printf("Failed to realloc heap inds.\n");
            exit(0);
        }

        heap->vs = (f32*) realloc(heap->vs, sizeof(f32) * heap->cap);
        if (heap->vs == NULL) {
            printf("Failed to realloc heap voltages.\n");
            exit(0);
        }
    }

    u64 i = heap->n;
    u64 pind = (i - 1LU) >> D;

    while (i > 0 && heap->times[pind] > it) {
        heap->times[i] = heap->times[pind];
        heap->inds[i] = heap->inds[pind];
        heap->vs[i] = heap->vs[pind];
        i = pind;
        pind = (i - 1LU) >> D;
    }

    heap->times[i] = it;
    heap->inds[i] = iind;
    heap->vs[i] = v;
    heap->n++;
}

int remmax(sigheap* heap, f32* ot, u32* oind, f32* v) {
    if (heap->n == 0) {
        return -1;
    }

    *ot = heap->times[0];
    *oind = heap->inds[0];
    *v = heap->vs[0];
    heap->n--;
    heap->times[0] = heap->times[heap->n];
    heap->inds[0] = heap->inds[heap->n];
    heap->vs[0] = heap->vs[heap->n];

    hfymax(heap, 0);

    return 0;
}

void hfymin(sigheap* heap, u64 i) {
    u64 min = i;
    for (u64 l = (i << D) + 1; l < (i << D) + D + 1; ++l) {
        if (l < heap->n && heap->times[l] < heap->times[min]) min = l;
    }

    if (min != i) {
        f32 tt = heap->times[0];
        u32 ti = heap->inds[0];
        f32 tv = heap->vs[0];
        heap->times[0] = heap->times[min];
        heap->inds[0] = heap->inds[min];
        heap->vs[0] = heap->vs[min];
        heap->times[min] = tt;
        heap->inds[min] = ti;
        heap->vs[min] = tv;

        hfymax(heap, min);
    }
}

void insmin(sigheap* heap, f32 it, u32 iind, f32 v) {
    if (heap->n == heap->cap) {
        heap->cap *= 2;
        heap->times = (f32*) realloc(heap->times, sizeof(f32) * heap->cap);
        if (heap->times == NULL) {
            printf("Failed to realloc heap times.\n");
            exit(0);
        }

        heap->inds = (u32*) realloc(heap->inds, sizeof(u32) * heap->cap);
        if (heap->inds == NULL) {
            printf("Failed to realloc heap inds.\n");
            exit(0);
        }

        heap->vs = (f32*) realloc(heap->vs, sizeof(f32) * heap->cap);
        if (heap->vs == NULL) {
            printf("Failed to realloc heap voltages.\n");
            exit(0);
        }
    }

    u64 i = heap->n;
    u64 pind = (i - 1LU) >> D;

    while (i > 0 && heap->times[pind] < it) {
        heap->times[i] = heap->times[pind];
        heap->inds[i] = heap->inds[pind];
        heap->vs[i] = heap->vs[pind];
        i = pind;
        pind = (i - 1LU) >> D;
    }

    heap->times[i] = it;
    heap->inds[i] = iind;
    heap->vs[i] = v;
    heap->n++;
}

int remmin(sigheap* heap, f32* ot, u32* oind, f32* v) {
    if (heap->n == 0) {
        return -1;
    }

    *ot = heap->times[0];
    *oind = heap->inds[0];
    *v = heap->vs[0];
    heap->n--;
    heap->times[0] = heap->times[heap->n];
    heap->inds[0] = heap->inds[heap->n];
    heap->vs[0] = heap->vs[heap->n];

    hfymin(heap, 0);

    return 0;
}