#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

typedef struct {
    int* data;
    size_t size;
} vector_int;

vector_int new_vector_int(size_t);

void delete_vector_int(vector_int* const);

void print_vector_int(vector_int const* const);

int find_sum_vector_int(vector_int const* const);

double find_sample_mean_vector_int(vector_int const* const);

int* find_max_value_vector_int(vector_int const* const);

int* find_min_value_vector_int(vector_int const* const);

size_t count_if_vector_int(bool (*compare)(int const, void const* const), 
    vector_int const* const, void const* const);
