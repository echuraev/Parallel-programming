#include <stdio.h>
#include <string.h>
#include <mpi.h>

#define SIZE 100

int main(int argc, char ** argv)
{
	int myrank;
	char * msg = malloc(sizeof(char) * SIZE);
	MPI_Init(&argc, &argv);
	
	MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

	if (myrank == 0)
	{
		strcpy(msg, "Hello MPI");
		printf("sent: %s\n", msg);
		MPI_Send(msg, strlen(msg), MPI_CHAR, 1, 0, MPI_COMM_WORLD);
	}
	else if (myrank == 1)
	{
		MPI_Recv(msg, SIZE, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("recv: %s\n", msg);
	}
	
	MPI_Finalize();
}
