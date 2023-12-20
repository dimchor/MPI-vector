#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <mpich-x86_64/mpi.h>

#include "../vector/vector_int.h"
#include "../vector/vector_double.h"

vector_int make_vector();

void split_array(
    size_t const,        // in:  array size
    unsigned int const,  // in:  divisor
    int** const,         // out: cnts: array of sizes
    int** const          // out: displs: array of offsets
);

void find_min_value(vector_int const* const, int* const);
void find_max_value(vector_int const* const, int* const);
void find_sum_value(vector_int const* const, int* const);

bool is_greater_than_double(int const, void const* const);
bool is_smaller_than_double(int const, void const* const);

void menu(
    vector_int const* const,
    double const,
    int const, 
    int const,
    double const, 
    vector_double const* const,
    size_t const,
    vector_int const* const
);

void print_everything(
    vector_int const* const,
    int const,
    int const,
    double const,
    int const, 
    int const,
    double const, 
    vector_double const* const,
    size_t const,
    vector_int const* const
);

#define ROOT 0

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    vector_int initial_v = make_vector();
    if (mpi_size > initial_v.size) {
        delete_vector_int(&initial_v);
        if (mpi_rank == ROOT) {
            printf("MPI_Comm_size and vector size mismatch. Please run this \
program with less processes or larger vector size.\n");
        }
        MPI_Finalize();
        return 0;
    }
    /*
    // for testing purposes
    if (mpi_rank == ROOT) {
        printf("input vector: ");
        print_vector_int(&initial_v);
    }
    */
    // Send initial_v
    int *sendcnts_initial_v, *displs_initial_v;
    split_array(initial_v.size, mpi_size, &sendcnts_initial_v, 
        &displs_initial_v);

    // make recv_v to store the numbers
    vector_int recv_v = new_vector_int(sendcnts_initial_v[mpi_rank]);
    MPI_Scatterv(initial_v.data, sendcnts_initial_v, displs_initial_v, MPI_INT, 
        recv_v.data, recv_v.size, MPI_INT, ROOT, MPI_COMM_WORLD);
    // Done

    /*
        +---------------------+
        | Find min, max, mean | 
        +---------------------+
    */

    // Final values
    int min_value, max_value, sum_value;

    find_min_value(&recv_v, &min_value);
    find_max_value(&recv_v, &max_value);
    find_sum_value(&recv_v, &sum_value);

    // Final values
    double const mean_value = (double) sum_value / initial_v.size;
    size_t greater_than_mean_counter, smaller_than_mean_counter;

    // Local values
    {
        size_t local_greater_than_mean_counter = count_if_vector_int(
            is_greater_than_double, &recv_v, &mean_value);
        size_t local_smaller_than_mean_counter = count_if_vector_int(
            is_smaller_than_double, &recv_v, &mean_value);

        MPI_Reduce(&local_greater_than_mean_counter, 
            &greater_than_mean_counter, 1, MPI_UNSIGNED_LONG, MPI_SUM, ROOT,
            MPI_COMM_WORLD);
        MPI_Reduce(&local_smaller_than_mean_counter, 
            &smaller_than_mean_counter, 1, MPI_UNSIGNED_LONG, MPI_SUM, ROOT,
            MPI_COMM_WORLD);
    }
    /*
    // for testing purposes
    if (mpi_rank == ROOT) {
        printf(
            "   min value: %d, max value: %d, mean value: %lf\n", 
            min_value, max_value, mean_value
        );
        printf("   numbers greater than mean: %ld\n   numbers smaller than \
mean: %ld\n", greater_than_mean_counter, smaller_than_mean_counter);
    }
    */

    /*  
        +---------------+
        | Find variance |
        +---------------+
    */

    double var_sum;
    {
        double local_var_sum = 0;
        for (int* it = recv_v.data; it < recv_v.data + recv_v.size; ++it) {
            local_var_sum += pow(*it - mean_value, 2);
        }
        MPI_Reduce(&local_var_sum, &var_sum, 1, MPI_DOUBLE, MPI_SUM, ROOT,
            MPI_COMM_WORLD);
    }
    double const var = var_sum / initial_v.size;
    /*
    // for testing purposes
    if (mpi_rank == ROOT) {
        printf("   variance: %lf\n", var);
    }
    */

    /*
        +----------+
        | vector D |
        +----------+
    */

    vector_double d = new_vector_double(initial_v.size);
    size_t max_value_d_position;
    {
        /*
            +---------------+
            | Fill vector d |
            +---------------+
        */
        vector_double local_d = new_vector_double(sendcnts_initial_v[mpi_rank]);

        for (size_t i = 0; i < local_d.size; ++i) {
            local_d.data[i] = (
                (recv_v.data[i] - min_value) / (double) (max_value - min_value)
                ) * 100;
        }

        MPI_Gatherv(local_d.data, local_d.size, MPI_DOUBLE, d.data, 
            sendcnts_initial_v, displs_initial_v, MPI_DOUBLE, ROOT, 
            MPI_COMM_WORLD);

        /*
            +----------------------+
            | Find max of vector d |
            +----------------------+
        */

        {
            double* local_max_d = find_max_value_vector_double(&local_d);
            size_t local_max_position = (local_max_d - local_d.data) + 
                displs_initial_v[mpi_rank];

            MPI_Reduce(&local_max_position, &max_value_d_position , 1, 
                MPI_UNSIGNED_LONG, MPI_MAX, ROOT, MPI_COMM_WORLD);
        }

        delete_vector_double(&local_d);
    }
    /*
    // for testing purposes
    if (mpi_rank == ROOT) {
        printf("vector d: ");
        print_vector_double(&d);
        printf("   max value position: %ld\n   → d[%ld] = %lf\n   → \
in[%ld] = %d\n", 
            max_value_d_position, max_value_d_position, 
            d.data[max_value_d_position], max_value_d_position, 
            initial_v.data[max_value_d_position]);
    }
    */


    /*
        +-------------+
        | Prefix sums |
        +-------------+
    */

    vector_int prefix_sums = new_vector_int(initial_v.size);

    {
        vector_int local_prefix_sums = new_vector_int(
            sendcnts_initial_v[mpi_rank]);

        for (size_t i = 0; i < local_prefix_sums.size; ++i) {
            for (size_t j = 0; j <= i; ++j) {
                local_prefix_sums.data[i] += recv_v.data[j];
            }
        }

        int last_element_prefix_sum;
        MPI_Scan(local_prefix_sums.data + local_prefix_sums.size - 1, 
            &last_element_prefix_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if (mpi_rank < mpi_size - 1) {
            MPI_Send(&last_element_prefix_sum, 1, MPI_INT, mpi_rank + 1, 0, 
                MPI_COMM_WORLD);
        }

        int correction_prefix;
        if (mpi_rank > ROOT) {
            MPI_Status status;
            MPI_Recv(&correction_prefix, 1, MPI_INT, mpi_rank - 1, 0, 
                MPI_COMM_WORLD, &status);

            for (size_t i = 0; i < local_prefix_sums.size; ++i) {
                local_prefix_sums.data[i] += correction_prefix;
            }
        }

        MPI_Gatherv(local_prefix_sums.data, local_prefix_sums.size, MPI_INT, 
            prefix_sums.data, sendcnts_initial_v, displs_initial_v, MPI_INT, 
            ROOT, MPI_COMM_WORLD);
        
        delete_vector_int(&local_prefix_sums);
    }
    /*
    // for testing purposes
    if (mpi_rank == ROOT) {
        printf("Prefix sums: ");
        print_vector_int(&prefix_sums);
    }
    */

    /*  
        +------+
        | Menu |
        +------+
    */

    if (mpi_rank == ROOT) {
        menu(&initial_v, mean_value, greater_than_mean_counter, 
            smaller_than_mean_counter, var, &d, max_value_d_position, 
            &prefix_sums);
        //print_everything(&initial_v, max_value, min_value, mean_value, 
        //    greater_than_mean_counter, smaller_than_mean_counter, var, &d, 
        //    max_value_d_position, &prefix_sums);
    }

    /*
        +----------+
        | Clean-up |
        +----------+
    */

    delete_vector_int(&prefix_sums);
    delete_vector_double(&d);

    // Delete received data
    delete_vector_int(&recv_v);
    free(sendcnts_initial_v); free(displs_initial_v);
    delete_vector_int(&initial_v);
    
    MPI_Finalize();
    return 0;
}

vector_int make_vector() {
    srand(time(NULL));

    vector_int v = new_vector_int(15);

    /*
    // random numbers
    {
        int* it = v.data;
        while (it < v.data + v.size) {
            *it++ = rand() % 20 + 10;
        }
    }
    */

    // ascending numbers
    {
        int* it = v.data;
        while(it < v.data + v.size) {
            *it++ = it - v.data + 1;
        }
    }

    return v;
}

void split_array(
    size_t const size,          // in:  array size
    unsigned int const divisor, // in:  divisor
    int** const cnts,           // out: cnts: array of sizes
    int** const displs          // out: displs: array of offsets
) {
    unsigned int const division = size / divisor;

    *displs = (int*) calloc(divisor, sizeof(int));
    *cnts = (int*) calloc(divisor, sizeof(int));
    for (size_t i = 0; i < divisor; ++i) {
        *(*displs + i) = i * division;
        *(*cnts + i) = division;
    }

    if (size % divisor == 0) {
        return;
    }

    *(*cnts + divisor - 1) += size - (division * divisor); 
}

void find_min_value(vector_int const* const v, int* const min_value) {
    int local_min_value = *find_min_value_vector_int(v);
    MPI_Allreduce(&local_min_value, min_value, 1, MPI_INT, MPI_MIN, 
        MPI_COMM_WORLD);
}

void find_max_value(vector_int const* const v, int* const max_value) {
    int local_max_value = *find_max_value_vector_int(v);
    MPI_Allreduce(&local_max_value, max_value, 1, MPI_INT, MPI_MAX, 
        MPI_COMM_WORLD);
}

void find_sum_value(vector_int const* const v, int* const sum_value) {
    int local_sum_value = find_sum_vector_int(v);
    MPI_Allreduce(&local_sum_value, sum_value, 1, MPI_INT, MPI_SUM, 
        MPI_COMM_WORLD);
}

bool is_greater_than_double(int const item, void const* const value) {
    if (item > *(double*) value) {
        return true;
    }
    return false;
}

bool is_smaller_than_double(int const item, void const* const value) {
    if (item < *(double*) value) {
        return true;
    }
    return false;
}

void menu(
    vector_int const* const initial_v,
    double const mean_value,
    int const greater_than_mean_counter, 
    int const smaller_than_mean_counter,
    double const variance, 
    vector_double const* const vector_d,
    size_t const max_value_d_position,
    vector_int const* const prefix_sums
) {
    char option;
    do {
        printf("\n--= Option menu =--\n  1. Show how many numbers are greater \
and smaller than the mean value.\n  2. Show the variance of the vector.\n  3. \
Show vector D\n  4. Largest value found in vector D.\n  5. Show prefix sums.\n\
  6. Exit.\nPlease select an option: ");
        option = getchar();
        while(getchar() != '\n');
        switch (option) {
        case '1':
            printf("Input vector: ");
            print_vector_int(initial_v);
            printf("   mean value: %lf\n   numbers greater than mean: %ld\n   \
numbers smaller than mean: %ld\n", mean_value, greater_than_mean_counter, 
                smaller_than_mean_counter);
            break;
        case '2':
            printf("Input vector: ");
            print_vector_int(initial_v);
            printf("   variance: %lf\n", variance);
            break;
        case '3':
            printf("Input vector: ");
            print_vector_int(initial_v);
            printf("Vector d: ");
            print_vector_double(vector_d);
            break;
        case '4':
            printf("Vector d: ");
            print_vector_double(vector_d);
            printf("   max value position: %ld\n   → d[%ld] = %lf\n   → \
input[%ld] = %d\n", 
                max_value_d_position, max_value_d_position, 
                vector_d->data[max_value_d_position], max_value_d_position, 
                initial_v->data[max_value_d_position]);
            break;
        case '5':
            printf("Prefix sums: ");
            print_vector_int(prefix_sums);
            break;
        case '6':
            break;
        default:
            printf("\nInvalid option. Press enter to continue...");
            while(getchar() != '\n');
            break;
        }
    } while (option != '6');
}

void print_everything(
    vector_int const* const initial_v,
    int const max_value,
    int const min_value,
    double const mean_value,
    int const greater_than_mean_counter, 
    int const smaller_than_mean_counter,
    double const variance, 
    vector_double const* const vector_d,
    size_t const max_value_d_position,
    vector_int const* const prefix_sums
) {
    // for testing purposes
    printf("Input vector: ");
    print_vector_int(initial_v);
    printf("   max value: %d\n   min value: %d\n", max_value, min_value);
    printf("   mean value: %lf\n   numbers greater than mean: %ld\n   \
numbers smaller than mean: %ld\n", mean_value, greater_than_mean_counter, 
        smaller_than_mean_counter);
    printf("   variance: %lf\n", variance);
    printf("Vector d: ");
    print_vector_double(vector_d);
    printf("   max value position: %ld\n   → d[%ld] = %lf\n   → \
input[%ld] = %d\n", 
        max_value_d_position, max_value_d_position, 
        vector_d->data[max_value_d_position], max_value_d_position, 
        initial_v->data[max_value_d_position]);
    printf("Prefix sums: ");
    print_vector_int(prefix_sums);
}

