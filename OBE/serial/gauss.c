/* Linear equation solution using gauss elimination method, serial version	  */
/* Ngakan Putu Ariastu														  */
/*	get time = ~265ms										                  */
/*  																	      */
/*  tested on Intel i7-3517u 4 core											  */
/******************************************************************************/

/* Source Code:                                                               */
/* Include all library we need                                                */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <Windows.h>

/* type def struct for picture data											  */
typedef struct EquationM {
	int row;
	int col;
	float* data;
}EquationM ;

typedef struct Vector {
	int row;
	float* data;
}Vector;


/* Declare all constant                                                       */
//#define N_EQ 5 // maximum equation
volatile DWORD dwStart;
#define M_cell(m,i,j) (m->data[i*m->col +j])


/* Declare all function prototype                                             */
int EquationNew(EquationM *m, int x, int y);
int getEchelon (EquationM *m);
int getResult (EquationM *m, Vector *r);
void PrintEquation(EquationM *m);
int VectorNew(Vector *m, int row);
void PrintVector(Vector *m)
{
	int i;
	for (i = 0; i < m->row; i++)
	{
		printf("%f\n", m->data[i]);
	}
}


/* Begin the Main Function                                                    */
int main()
{
	EquationM eq;
	Vector res;
	int n_eq, n_var;
	int i,j;
	dwStart = GetTickCount();
	printf("input number of variable: ");
	scanf ("%d",&n_var);
	printf("input number of equation: ");
	scanf ("%d",&n_eq);
	
	
	if(n_var != n_eq)
	{
		printf("cant solve equation, number of unknown variable and equation not match\n");
		return 0;
	} 
	
	if (EquationNew(&eq, n_eq, n_var+1) != 1)
		printf("creating equation matrix failed\n");
	
	if (VectorNew(&res,n_var) != 1)
		printf("creating result failed\n");
		
	printf("Enter the elements of augmented matrix row-wise:\n\n");	
	for(i = 0; i < eq.row  ; i++)
	{
		for(j = 0; j <  eq.col; j++)
		{
			printf("A[%d][%d] : ", i,j);
			scanf("%f",&eq.data[i*eq.col + j]);
		}
	}
	
	
	
	if (getEchelon(&eq) != 1)
		printf ("echelon form failed\n");
	if (getResult(&eq,&res) != 1)
		printf ("get final result failed\n");
	
	
	//PrintEquation(&eq);
	PrintVector(&res);
	printf("time taken, %d milliseconds\n", GetTickCount() - dwStart);

	return(0);
	/* End Main Function                                                          */
}

int EquationNew(EquationM *m, int x, int y)
{
	m->row = x;
	m->col = y;
	m->data = (float*)malloc(m->row*m->col*sizeof(float));

	if (m->data)
		return 1;
	else
		return 0;
}

int VectorNew(Vector *m, int row)
{
	m->row = row;
	m->data = (float*)malloc(m->row*sizeof(float));

	if (m->data)
		return 1;
	else
		return 0;
}


int getEchelon (EquationM *m)
{
	int i,j,k;
	float c;
	for(j = 0; j < m->row ; j++) /* loop for the generation of upper triangular matrix*/
    {
        for(i = 0 ; i < m->row; i++)
        {
            if(i>j)
            {
				if (M_cell(m,j,j) == 0 )
					return 0;
				
                c = M_cell(m,i,j) / M_cell(m,j,j);
                for(k=0; k < m->col; k++)
                {
                    M_cell(m,i,k) = M_cell(m,i,k)-c*M_cell(m,j,k);
                }
            }
        }
    }
	return 1;
}

int getResult (EquationM *m, Vector *r)
{
	int i,j;
	float sum = 0.0;
	
	r->data[r->row -1]  = m->data[(m->row-1)*m->col + m->col-1] / m->data[(m->row-1)*m->col + m->col-2];
	
	for(i=m->row-2; i>=0; i--)
    {
        sum=0;
        for(j=i+1; j<=m->row-1; j++)
        {
            sum = sum + m->data[i * m->col + j] *r->data[j] ;
        }
        r->data[i] = (m->data[i*m->col + m->col-1]-sum) /  m->data[i*m->col+i];
    }
	return 1;
}



void PrintEquation(EquationM *m)
{
	int i, j;
	for (i = 0; i < m->row; i++)
	{
		for (j = 0; j < m->col; j++)
		{
			printf("%f	", m->data[i*m->col + j]);
		}
		printf("\n");
	}
}










