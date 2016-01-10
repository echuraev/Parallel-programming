#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <mpi.h>

#define MULT = 0
#define SUM = 1

void parse_params(int argc, char ** argv, double *** a, int * colA, int * rowA, double *** b, int * colB, int * rowB, int * oper);
void sum(double *** c, double ** a, int colA, int rowA, double ** b, int colB, int rowB);
void mult(double *** c, double ** a, int colA, int rowA, double ** b, int colB, int rowB);

int main(int argc, char ** argv)
{
	int oper;
	int colA, rowA, colB, rowB;
	double ** a;
	double ** b;
	double ** c;
	parse_params(argc, argv, &a, &colA, &rowA, &b, &colB, &rowB, &oper);
	printf("a\n");
	int i, j;
	for (i = 0; i < rowA; ++i)
	{
		for (j = 0; j < colA; ++j)
		{
			printf("%lf\t", a[i][j]);
		}
		printf("\n");
	}
	printf("\nb\n");
	for (i = 0; i < rowB; ++i)
	{
		for (j = 0; j < colB; ++j)
		{
			printf("%lf\t", b[i][j]);
		}
		printf("\n");
	}
	return 0;
}

void parse_params(int argc, char ** argv, double *** a, int * colA, int * rowA, double *** b, int * colB, int * rowB, int * oper)
{
	if (argc < 4)
	{
		printf("Use: %s file1 file2 oper\n", argv[0]);
		printf("\tfile1\tFile with first matrix\n");
		printf("\tfile2\tFile with second matrix\n");
		printf("\toper\tOne symbol: * of +\n");
		exit(1);
	}
	FILE* aFile = fopen(argv[1], "rt");
	if (aFile == NULL)
	{
		perror("Error! Cannot open file!\n");
		exit(2);
	}
	char * line = NULL;
	int len = 0;
	*rowA = 0;
	*colA = 0;
	int i, j;
	int read;
	while ((read = getline(&line, &len, aFile)) != -1) 
	{
		if (*rowA == 0)
		{
			for (i = 0; i < read; ++i)
			{
				if (line[i] == ' ') 
				{
					*colA++;
				}
			}
		}
		*rowA++;
	}
	fseek(aFile, 0, SEEK_SET);
	if (*rowA == 0 || *colA == 0)
	{
		perror("Error! Invalid matrix size!\n");
		exit(3);
	}
	*a = malloc(sizeof(double) * *rowA);
	for (i = 0; i < *rowA; ++i)
	{
		*a[i] = malloc(sizeof(double) * *colA);
		for (j = 0; j < *colA; ++j)
		{
			if (fscanf(aFile, "%lf", *a[i][j]) == EOF && (i < *rowA - 1 || j < *colA - 1))
			{
				perror("Error! Cannot read matrix!\n");
				exit(4);
			}
		}
	}
	FILE* bFile = fopen(argv[2], "rt");
	if (bFile == NULL)
	{
		perror("Error! Cannot open file!\n");
		exit(2);
	}
	*rowB = 0;
	*colB = 0;
	while ((read = getline(&line, &len, aFile)) != -1) 
	{
		if (*rowB == 0)
		{
			for (i = 0; i < read; ++i)
			{
				if (line[i] == ' ') 
				{
					*colB++;
				}
			}
		}
		*rowB++;
	}
	fseek(bFile, 0, SEEK_SET);
	if (*rowB == 0 || *colB == 0)
	{
		perror("Error! Invalid matrix size!\n");
		exit(3);
	}
	*b = malloc(sizeof(double) * *rowB);
	for (i = 0; i < *rowB; ++i)
	{
		*b[i] = malloc(sizeof(double) * *colB);
		for (j = 0; j < *colB; ++j)
		{
			if (fscanf(bFile, "%lf", *b[i][j]) == EOF && (i < *rowB - 1 || j < *colB - 1))
			{
				perror("Error! Cannot read matrix!\n");
				exit(4);
			}
		}
	}
	if (argv[3][0] == '*')
	{
		*oper = MULT;
	}
	else if (argv[3][0] == '+')
	{
		*oper = SUM;
	}
	else
	{
		perror("Error! Unknown operation!");
		exit(5);
	}
}

