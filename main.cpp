#include <mpi.h>
#include <cstdio>
#include <fstream>
#include <cstring>

#include "wtime.h"
#include "psgl.h"

using namespace std;

// global members
char *data_graph_file_name, *query_graph_file_name, *preset_file_name;
double time1;

// Function declarations
void help();

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	int worker_count = 1;
	int worker_rank  = 0;

	MPI_Comm_size(MPI_COMM_WORLD, &worker_count);
	MPI_Comm_rank(MPI_COMM_WORLD, &worker_rank);

	if(argc != 4){
		printf("Wrong number of parameters\n");
		help();
		exit(1);
	}

	// Read query and data file names
	data_graph_file_name  	= argv[1];
	query_graph_file_name 	= argv[2];
	preset_file_name  		= argv[3];

	if(worker_rank == 0) printf("Data graph: %s and Query Graph %s\n", data_graph_file_name, query_graph_file_name);
	
	// Find Automorphs
	bool break_automorph = false;
	Psgl myAutomorphFinder(data_graph_file_name, query_graph_file_name, break_automorph, preset_file_name);
	myAutomorphFinder.generic_query_proc();

	if(worker_rank == 0){printf("**** subgraph Matching completed *****\n");}
	time1 = wtime();
	
	MPI_Finalize();
	return 0;
}

void help()
{
	std::cout << "Usage:" << "1. application.bin\n"<< "2. dataGraph\n" << "3. queryGraph\n" << "4. preset file\n";
}

