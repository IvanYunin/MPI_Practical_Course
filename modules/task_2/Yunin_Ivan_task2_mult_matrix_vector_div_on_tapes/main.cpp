// Copyright  2018 Ivan Yunin
#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <string>
#include <ctime>
#include <cstdlib>
#include <iostream>


/* Description
Task 2 "Matrix vector multiplication (Division on tapes)"
The task is to multiply matrix by vector
Matrix has a dim [n x m]
Vector has a dim [m x 1]
Result of multiply is vector with dim [n x 1]

(*) Each element of result vector is scalar multiplication vector by row in matrix 

    Parallelism
    
Quantity processes = number_process

Matrix divided on tapes 
Quantity rows in tape is  k = n / number_process;    (rounded up)
Each process receives part of matrix with dim [k x m]
Each process calculate  multiplication vector by  PART of matrix
Result of work process is part of result vector

    Variables
col_num -> m; row_num -> n; sub_row_num-> k
proc_num -> number_process
flag_out -> if true - show serial result vector and parallel result vector

    Functions
check() - Checking for vector equality
scal_mult() - Scalar multiplication vector by other vector 
*/

int check(double* first, double* second, int size) {
    int res = 0;
    for (int i = 0; i < size; i++) {
        if (fabs(first[i]-second[i]) > 0.0000001) {
            res = 1;
            break;
        }
    }
    return res;
}


double scal_mult(double* a, double* b, int size) {
    double res = 0;
    for (int i = 0; i < size; i++) {
        res+=a[i]*b[i];
    }
    return res;
}


int main(int argc, char*argv[]) {
    int col_num = 100, row_num = 100, sub_row_num;
    int proc_num, proc_id, flag;
    int flag_out = 0;
    double *vector, *matrix = nullptr, *sub_matrix;
    double *serial_res = nullptr, *parallel_res = nullptr, *sub_parallel_res;
    double serial_time =0.0, parallel_time = 0.0;
    std::srand(static_cast<int>(time(0)));

    if (argc > 2) {
        row_num = atoi(argv[1]);
        col_num = atoi(argv[2]);
    }
    if (argc > 3) {
        flag_out = atoi(argv[3]);
    }
    MPI_Init(&argc, &argv);
    MPI_Initialized(&flag);
    if (!flag) {
        std::cout << "Init MPI Error";
        return 0;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    sub_row_num = static_cast<int>(ceil(static_cast<double>(row_num)/
                                        (static_cast<double>(proc_num))));
    // quantity rows for one process

// memory alloc for all processes
    vector = new double[col_num];
    sub_matrix = new double[sub_row_num*col_num];
    sub_parallel_res = new double[sub_row_num];
//  end memory alloc for all processes

    if (proc_id == 0) {
        int tail = proc_num*sub_row_num-row_num;
        /* tail-> if the quantity rows is divided between the processes 
        is not entirely (for correct work Gather() )*/

        // init vector and matrix (memory and values)
        matrix = new double[row_num*col_num];
        serial_res = new double[row_num];
        parallel_res = new double[row_num+tail];
        for (int i = 0; i < col_num; i++) {
            vector[i] = -(std::rand()%100) + (std::rand()%200)/13.0;
            for (int j = 0; j < row_num; j++) {
                matrix[j*col_num+i] = ((std::rand()%300)/17.0)-
                                      (std::rand()%100)/3.0;
            }
        }
        // end init vector and matrix (memory and values)
        // serial multiplication
        serial_time = MPI_Wtime();
        for (int i = 0; i < row_num; i++) {
            serial_res[i] = scal_mult(vector, matrix+col_num*i, col_num);
            // (*) look description
        }
        serial_time = MPI_Wtime() - serial_time;
        std::cout << "serial time: " << serial_time << '\n';
        // serial multiplication end
        parallel_time = MPI_Wtime();
    }

    // (#)(start) code executed by all processes
    MPI_Bcast(vector, col_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // MPI_Bcast send vector to all procceses
    MPI_Scatter(matrix, sub_row_num*col_num, MPI_DOUBLE,
       sub_matrix, sub_row_num*col_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    /* send parts of matrix to all processes (each procces gets matrix
    with dim [sub_row_num x col_num])*/
    // begin calculate sub result
    for (int i = 0; i < sub_row_num; i++) {
            sub_parallel_res[i] = scal_mult(vector, (sub_matrix+col_num*i),
                                            col_num);  // (*) look description
        }
    // end calculate sub result
    MPI_Gather(sub_parallel_res, sub_row_num, MPI_DOUBLE, parallel_res,
               sub_row_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    /* collection result data in proc#0 (each procces send 
    part of result_vector, part had dim [sub_row_num x 1])*/

    // (#)(finish) code executed by all processes

    if (proc_id == 0) {
        parallel_time = MPI_Wtime() - parallel_time;
        if (flag_out) {
            std::cout << "Serial result vector" << '\n';
            for (int i = 0; i < row_num; i++) {
                std::cout << serial_res[i] << " ";
            }
            std::cout << '\n' << "Parallel result vector" << '\n';
            for (int i = 0; i < row_num; i++) {
                std::cout << parallel_res[i] << " ";
            }
        }

        std::cout << '\n' << "Parralel time: " << parallel_time << '\n';
        if (check(serial_res, parallel_res, row_num))
            std::cout << "Error: vectors are not equal" << '\n';
    }
    if (matrix       != nullptr ) { delete[]matrix;       }
    if (serial_res   != nullptr ) { delete[]serial_res;   }
    if (parallel_res != nullptr ) { delete[]parallel_res; }
    delete[]sub_parallel_res;
    delete[]vector;
    delete[]sub_matrix;
    MPI_Finalize();
    return 0;
}
