#include "vector_int.h"

vector_int new_vector_int(size_t size) {
    vector_int v;
    v.size = size;
    v.data = (int*) calloc(size, sizeof(int));
    return v;
}

void delete_vector_int(vector_int* const v) {
    free(v->data);
    v->data = NULL;
    v->size = 0;
}

void print_vector_int(vector_int const* const v) {
    int* item = v->data;
    printf("[ ");
    while(item < v->data + v->size) {
        printf("%d ", *item++);
    }
    printf("]\n");
}

int find_sum_vector_int(vector_int const* const v) {
    int sum = 0, *item = v->data;
    while(item < v->data + v->size) {
        sum += *item++;
    }
    return sum;
}

double find_sample_mean_vector_int(vector_int const* const v) {
    return (double) find_sum_vector_int(v) / v->size;
}

int* find_max_value_vector_int(vector_int const* const v) {
    int* max_value = v->data, *item = v->data;
    while(item < v->data + v->size) {
        if (*item > *max_value) {
            max_value = item;
        }
        ++item;
    }
    return max_value;
}

int* find_min_value_vector_int(vector_int const* const v) {
    int* min_value = v->data, *item = v->data;
    while(item < v->data + v->size) {
        if (*item < *min_value) {
            min_value = item;
        }
        ++item;
    }
    return min_value;
}

size_t count_if_vector_int(
    bool (*compare)(int const, void const* const), // compare function
    vector_int const* const v,                 // vector_int
    void const* const args                 // any args
) {
    size_t counter = 0;
    {
        int* item = v->data;
        while (item < v->data + v->size) {
            if (compare(*item++, args)) {
                counter++;
            }
        }
    }
    return counter;
}

/*
double find_sum_of_distances_vector(vector_int const* const v) {
    double const m = find_sample_mean_vector(v);
    double sum = 0;
    int* item = v->data;
    while(item < v->data + v->size) {
        sum += pow(*item++ - m, 2);
    }
    return sum;
}

double find_var_vector(vector_int const* const v) {
    return (double) find_sum_of_distances_vector(v) / v->size;
}
*/
