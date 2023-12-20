#include "vector_double.h"

vector_double new_vector_double(size_t size) {
    vector_double v;
    v.size = size;
    v.data = (double*) calloc(size, sizeof(double));
    return v;
}

void delete_vector_double(vector_double* const v) {
    free(v->data);
    v->data = NULL;
    v->size = 0;
}

void print_vector_double(vector_double const* const v) {
    double* item = v->data;
    printf("[ ");
    while(item < v->data + v->size) {
        printf("%lf ", *item++);
    }
    printf("]\n");
}

double find_sum_vector_double(vector_double const* const v) {
    double sum = 0, *item = v->data;
    while(item < v->data + v->size) {
        sum += *item++;
    }
    return sum;
}

double find_sample_mean_vector_double(vector_double const* const v) {
    return (double) find_sum_vector_double(v) / v->size;
}

double* find_max_value_vector_double(vector_double const* const v) {
    double* max_value = v->data, *item = v->data;
    while(item < v->data + v->size) {
        if (*item > *max_value) {
            max_value = item;
        }
        ++item;
    }
    return max_value;
}

double* find_min_value_vector_double(vector_double const* const v) {
    double* min_value = v->data, *item = v->data;
    while(item < v->data + v->size) {
        if (*item < *min_value) {
            min_value = item;
        }
        ++item;
    }
    return min_value;
}

size_t count_if_vector_double(
    bool (*compare)(double const, void const* const), // compare function
    vector_double const* const v,                 // vector
    void const* const args                 // any args
) {
    size_t counter = 0;
    {
        double* item = v->data;
        while (item < v->data + v->size) {
            if (compare(*item++, args)) {
                counter++;
            }
        }
    }
    return counter;
}

