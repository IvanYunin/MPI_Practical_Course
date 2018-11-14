#include <iostream>
#include <mpi.h>
#include <string>
#include <ctime>
#include <cstdlib>
#include <assert.h>

using namespace std;
double max_func(double* first, int size){
	double max = first[0];
	for(int i = 0; i<size;i++){
	if (first[i]>max)
		max = first[i];
	}
	return max;
}

int main (int argc, char* argv[]){
	int tail_flag=0; //for сheck the existence  vector "tail"
	int size=42, proc_num,proc_id, block_size,flag;
	double *vector,*temp_vector;
	double serial_time, parallel_time;
	double serial_max, parallel_max, temp_max;
	MPI_Status status;
	srand((int)time(0));
	
	if (argc>1){
		size = atoi(argv[1]);
	}

	MPI_Init(&argc, &argv);
	MPI_Initialized(&flag);
	if (!flag){
		cout << "Error";
        return 0;
	}
	MPI_Comm_size(MPI_COMM_WORLD,&proc_num); 
	MPI_Comm_rank(MPI_COMM_WORLD,&proc_id);
	block_size = size/proc_num;
	if (proc_id == 0){
		//init vector
		vector = new double[size];
		for(int i = 0; i<size; i++)
			vector[i]=-(rand()%23) + (rand()%200)/13.0;
		//serial max search
		serial_time = MPI_Wtime();
		serial_max = max_func(vector,size); 
		serial_time = MPI_Wtime() - serial_time;
		cout<<"serial  max: "<<serial_max<<endl;
		cout<<"serial time: "<<serial_time<<endl;
		//serial max search end
		parallel_time = MPI_Wtime();
	}
	
	temp_vector = new double[block_size];
	MPI_Scatter(vector, block_size, MPI_DOUBLE,
       temp_vector, block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	temp_max = max_func(temp_vector,block_size);
	MPI_Reduce(&temp_max, &parallel_max, 1, MPI_DOUBLE,
                  MPI_MAX, 0, MPI_COMM_WORLD);
	if (proc_id == 0){
		parallel_time = MPI_Wtime() - parallel_time;
		cout<<"Parralel max: "<<parallel_max<<endl;
		cout<<"Parralel time: "<<parallel_time<<endl;
		delete[]vector;
	}
	delete[]temp_vector;
	MPI_Finalize();
	return 0;
}
