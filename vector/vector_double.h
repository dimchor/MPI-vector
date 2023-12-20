#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

typedef struct {
    double* data;
    size_t size;
} vector_double;

vector_double new_vector_double(size_t);

void delete_vector_double(vector_double* const);

void print_vector_double(vector_double const* const);

double find_sum_vector_double(vector_double const* const);

double find_sample_mean_vector_double(vector_double const* const);

double* find_max_value_vector_double(vector_double const* const);

double* find_min_value_vector_double(vector_double const* const);

size_t count_if_vector_double(bool (*compare)(double const, void const* const), 
    vector_double const* const, void const* const);
