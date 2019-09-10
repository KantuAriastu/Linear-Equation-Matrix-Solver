#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h> 

typedef struct Matrix {
	int row;
	int col;
	int* data;

}Matrix;

typedef struct Vector {
	int row;
	int* data;

}Vector;

int VectorNew(Vector *m, int row)
{
	m->row = row;
	m->data = (int*)calloc(m->row, sizeof(int));

	if (m->data)
		return 1;
	else
		return 0;
}


int MatrixNew(Matrix *m, int row, int col)
{
	m->row = row;
	m->col = col;
	m->data = (int*)calloc(m->row*m->col , sizeof(int));

	if (m->data)
		return 1;
	else
		return 0;
}


void PrintMatrix(Matrix *m)
{
	int i, j;
	for (i = 0; i < m->row; i++)
	{
		for (j = 0; j < m->col; j++)
		{
			printf("%d	", m->data[i*m->col + j]);
		}
		printf("\n");
	}
}

void PrintVector(Vector *m)
{
	int i;
	for (i = 0; i < m->row; i++)
	{
		printf("%d\n", m->data[i]);
	}
}


#define NRM 10000		       /* number of rows in matrix */
#define NCM 10000              /* number of columns in matrix */
#define NRV 10000              /* number ofrow in vector */
#define FROM_MASTER 1          /* setting a message type */
#define FROM_WORKER 2          /* setting a message type */

int main(int argc, char *argv[])
{
	clock_t start, end;
	start = clock();
	double cpu_time_used;
	int	numtasks,              /* number of tasks in partition */
		taskid,                /* a task identifier */
		averow,				   /* used to determine rows sent to each worker */
		i,j,rc;				   /* misc */
	Vector Vect;			   /* global vector used in every process */

	
	/* Initializing MPI environment */
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	if (numtasks < 2) {
		printf("Need at least two MPI tasks. Quitting...\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
		exit(1);
	}
	
	/* Initializing global vector */
	if (VectorNew(&Vect, NRV) != 1)
	{
		printf("creating global vect fail\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
		exit(1);
	}
	for (i = 0; i < Vect.row; i++)
	{
		Vect.data[i] = 1;
	}
	
	/* dividing work per worker */
	averow = NRM / 4;

	if (taskid == 0)
	{
		/* Create local main Marix and Vector */
		Matrix Mat;
		Vector Res;
		if (MatrixNew(&Mat, NRM, NCM) != 1)
		{
			printf("creating local main matrix on worker 0 fail\n");
			MPI_Abort(MPI_COMM_WORLD, rc);
			exit(1);
		}
		if (VectorNew(&Res, NRM) != 1)
		{
			printf("creating local main vect on worker 0 fail\n");
			MPI_Abort(MPI_COMM_WORLD, rc);
			exit(1);
		}

		/* Initializing local main vector and matrix */ 
		for (i = 0; i < Mat.row; i++)
		{
			for (j = 0; j <Mat.col; j++)
			{
				Mat.data[i*Mat.col + j] = 1;
			}
		}
		
	
		/* sending other row data to other worker */ 
		MPI_Send(Mat.data + averow*Mat.col * 1, averow*Mat.col, MPI_INT, 1, 0, MPI_COMM_WORLD);	 // to worker 1
		MPI_Send(Mat.data + averow*Mat.col * 2, averow*Mat.col, MPI_INT, 2, 0, MPI_COMM_WORLD);  // to worker 2
		MPI_Send(Mat.data + averow*Mat.col * 3, averow*Mat.col, MPI_INT, 3, 0, MPI_COMM_WORLD);  // to worker 3

		/* starting work on this worker */
		printf("worker %d working on row %d -> %d\n", taskid, taskid*averow, taskid*averow + averow - 1);
		for (i = 0; i < averow; i++)
		{
			for (j = 0; j < Mat.col; j++)
			{
				Res.data[i] += Mat.data[i*Mat.col + j] * Vect.data[j];
			}
		}
		printf("worker %d finished calculating, receiving data from other workers\n", taskid);

		/* Receiving data from other worker */
		MPI_Recv(Res.data+1*averow, averow, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);	// from worker 1
		MPI_Recv(Res.data+2*averow, averow, MPI_INT, 2, 0, MPI_COMM_WORLD, &status);	// from worker 2
		MPI_Recv(Res.data+3*averow, averow, MPI_INT, 3, 0, MPI_COMM_WORLD, &status);	// from worker 3
		printf("finished receiving, matrix multiplication done\n\n");

		/* printing Result */
		//PrintVector(&Res);
		free(Mat.data);
		free(Res.data);
		Mat.data = NULL;
		Res.data = NULL;

	}

	if (taskid == 1)
	{
		/* local Matrix to store data from worker 0 */
		Matrix local_Mat1;
		if (MatrixNew(&local_Mat1, averow, NCM) != 1)
		{
			printf("creating local matrix on worker 1 fail\n");
			MPI_Abort(MPI_COMM_WORLD, rc);
			exit(1);
		}
		MPI_Recv(local_Mat1.data, averow*local_Mat1.col, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	
		/* local Vector to store result to send to worker 0  */
		Vector local_Vect1;
		if (VectorNew(&local_Vect1, averow) != 1)
		{
			printf("creating local vect on worker 1 fail\n");
			MPI_Abort(MPI_COMM_WORLD, rc);
			exit(1);
		}

		/* calculating local multiplication */
		printf("worker %d working on row %d -> %d\n", taskid, taskid*averow, taskid*averow + averow - 1);
		for (i = 0; i < averow ; i++)
		{
			for (j = 0; j < local_Mat1.col; j++)
			{
				local_Vect1.data[i] += local_Mat1.data[i*local_Mat1.col + j] * Vect.data[j];
			}
		}
		printf("worker %d finished calculating, send to central 0\n\n", taskid);
		
		/* send result to worker 0 */
		MPI_Send(local_Vect1.data, averow, MPI_INT, 0, 0, MPI_COMM_WORLD);
		free(local_Mat1.data);
		free(local_Vect1.data);
		local_Mat1.data = NULL;
		local_Vect1.data = NULL;
	}

	if (taskid == 2)
	{
		/* local Matrix to store data from worker 0 */
		Matrix local_Mat2;
		if (MatrixNew(&local_Mat2, averow, NCM) != 1)
		{
			printf("creating local matrix on worker 2 fail\n");
			MPI_Abort(MPI_COMM_WORLD, rc);
			exit(1);
		}

		MPI_Recv(local_Mat2.data, averow*local_Mat2.col, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	
		/* local Vector to store result to send to worker 0  */
		Vector local_Vect2;
		if (VectorNew(&local_Vect2, averow) != 1)
		{
			printf("creating local vector on worker 2 fail\n");
			MPI_Abort(MPI_COMM_WORLD, rc);
			exit(1);
		}

		/* calculating local multiplication */
		printf("worker %d working on row %d -> %d\n", taskid, taskid*averow, taskid*averow + averow - 1);
		for (i = 0; i < averow; i++)
		{
			for (j = 0; j < local_Mat2.col; j++)
			{
				local_Vect2.data[i] += local_Mat2.data[i*local_Mat2.col + j] * Vect.data[j];
			}
		}
		printf("worker %d finished calculating, send to central 0\n\n", taskid);
		
		/* send result to worker 0 */
		MPI_Send(local_Vect2.data, averow, MPI_INT, 0, 0, MPI_COMM_WORLD);
		free(local_Mat2.data);
		free(local_Vect2.data);
		local_Mat2.data = NULL;
		local_Vect2.data = NULL;
	}

	if (taskid == 3)
	{
		/* local Matrix to store data from worker 0 */
		Matrix local_Mat3;
		if (MatrixNew(&local_Mat3, averow, NCM) != 1)
		{
			printf("creating local matrix on worker 3 fail\n");
			MPI_Abort(MPI_COMM_WORLD, rc);
			exit(1);
		}
		MPI_Recv(local_Mat3.data, averow*local_Mat3.col, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		
		/* local Vector to store result to send to worker 0  */
		Vector local_Vect3;
		if (VectorNew(&local_Vect3, averow) != 1)
		{
			printf("creating local vector on worker 3 fail\n\n");
			MPI_Abort(MPI_COMM_WORLD, rc);
			exit(1);
		}

		/* calculating local multiplication */
		printf("worker %d working on row %d -> %d\n", taskid, taskid*averow, taskid*averow + averow - 1);
		for (i = 0; i < averow; i++)
		{
			for (j = 0; j < local_Mat3.col; j++)
			{
				local_Vect3.data[i] += local_Mat3.data[i*local_Mat3.col + j] * Vect.data[j];
			}
		}
		printf("worker %d finished calculating, send to central 0\n\n", taskid);

		/* send result to worker 0 */
		MPI_Send(local_Vect3.data, averow, MPI_INT, 0, 0, MPI_COMM_WORLD);
		free(local_Mat3.data);
		free(local_Vect3.data);
		local_Mat3.data = NULL;
		local_Vect3.data = NULL;
	}


	free(Vect.data);
	Vect.data = NULL;
	
	MPI_Finalize();

	end = clock();
	cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
	printf("%f seconds\n", cpu_time_used);
	
}