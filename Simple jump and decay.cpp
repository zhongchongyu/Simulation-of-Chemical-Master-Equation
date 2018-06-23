#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<fstream>
using namespace std;

#define col 100
#define row 299
#define MAX_OUTPUT_LENGTH  2400000

void Gillespie(double stoich_matrix[row][col], double tspan[2], double x_start[col], double p, double lamda)
{
	double** result = (double**)malloc(sizeof(double*)* MAX_OUTPUT_LENGTH);	//第一列存放时间，其他列存放反应物数量
	for (int i = 0; i < MAX_OUTPUT_LENGTH; i++)
	{
		result[i] = (double*)malloc(sizeof(double)*(col + 1));
	}
	result[0][0] = tspan[0];	//时间赋初值
	for (int i = 1; i <= col; i++)	//反应物赋初值
	{
		result[0][i] = x_start[i - 1];
	}

	int rxn_count = 0; double a[row];

	while (result[rxn_count][0] < tspan[1])
	{
		// 设置初始a向量
		a[0] = lamda * 15;
		for (int i = 1; i < col; i++)
		{
			a[i] = lamda * result[rxn_count][i];

		}
		for (int i = col; i < 2 * col - 1; i++)
		{
			a[i] = lamda * result[rxn_count][2 * col - i];
		}
		for (int i = 2 * col - 1; i < row; i++)
		{
			a[i] = p * result[rxn_count][i - (2 * col - 2)];
		}

		double a_sum = 0;
		for (int i = 0; i < row; i++)
		{
			a_sum += a[i];
		}

		double r[2];
		r[0] = (rand() + 1) / (RAND_MAX + 2.0); r[1] = (rand() + 1) / (RAND_MAX + 2.0);
		double tau = -log(r[0]) / a_sum;

		// Sample identity of earliest reaction channel to fire (mu)
		int mu = 0; double s = a[0]; double r0 = r[1] * a_sum;

		while (s < r0)
		{
			mu += 1;
			s += a[mu];
		}

		if (rxn_count + 3 > MAX_OUTPUT_LENGTH)
		{
			printf("超出最大迭代次数");
			break;
		}

		result[rxn_count + 1][0] = result[rxn_count][0] + tau;
		for (int i = 1; i <= col; i++)
		{
			result[rxn_count + 1][i] = result[rxn_count][i] + stoich_matrix[mu][i - 1];
		}
		rxn_count += 1;
	}

	if (result[rxn_count][0] > tspan[1])
	{
		result[rxn_count][0] = tspan[1];
		for (int i = 1; i <= col; i++)
		{
			result[rxn_count][i] = result[rxn_count - 1][i];
		}
	}

	ofstream outfile;
	outfile.open("resule.txt");
	for (int i = rxn_count - 1000; i < rxn_count; i++)
	{
		for (int j = 0; j < col + 1; j++)
		{
			outfile << result[i][j] << " ";
		}
		outfile << endl;
	}
	outfile.close();

}




int main(void)
{
	double lamda = 0.9;  //lamda是一个跟dx有关的常数  这里lamda = p‘ * 0.01
	double p = 0.01;
	double tspan[2]; tspan[0] = 0; tspan[1] = 10000;
	double x_start[col];
	for (int i = 0; i < col; i++)
	{
		x_start[i] = 0;
	}


	double stoich_matrix[3 * col - 1][col];
	//这里初始化矩阵
	for (int i = 0; i < 3 * col - 1; i++)
	{
		for (int j = 0; j < col; j++)
		{
			stoich_matrix[i][j] = 0;
		}
	}
	//输入矩阵
	stoich_matrix[0][0] = 1;
	for (int i = 1; i < col; i++)
	{
		stoich_matrix[i][i] = 1;
		stoich_matrix[i][i - 1] = -1;
	}

	for (int i = col; i < 2 * col - 1; i++)
	{
		stoich_matrix[i][2 * col - 1 - i] = -1;
		stoich_matrix[i][2 * col - 2 - i] = 1;
	}

	for (int i = 2 * col - 1; i < row; i++)
	{

		stoich_matrix[i][i - (2 * col - 1)] = -1;

	}

	Gillespie(stoich_matrix, tspan, x_start, p, lamda);

	return 0;
}