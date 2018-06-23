#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<fstream>
using namespace std;

#define k 4		//有k种反应发生
#define num_species 2		//有num_species中反应物
#define MAX_OUTPUT_LENGTH 10000	//最大迭代次数




class parm
{
public:
	double kr, kp, gr, gp;
};

void Gillespie(double stoich_matrix[4][2], double tspan[2], double x0[2], parm p)
{
	double** result = (double**)malloc(sizeof(double*)*MAX_OUTPUT_LENGTH);	//第一列存放时间，其他列存放反应物数量
	for (int i = 0; i < MAX_OUTPUT_LENGTH; i++)
	{
		result[i] = (double*)malloc(sizeof(double)*(num_species + 1));
	}
	result[0][0] = tspan[0];	//时间赋初值
	for (int i = 1; i <= num_species; i++)	//反应物赋初值
	{
		result[0][i] = x0[i - 1];
	}

	int rxn_count = 0; double a[4];

	while (result[rxn_count][0] < tspan[1])
	{
		double mRNA = result[rxn_count][1];
		double protein = result[rxn_count][2];
		a[0] = p.kr; a[1] = p.kp * mRNA; a[2] = p.gr * mRNA; a[3] = p.gp * protein;
		double a_sum = a[0] + a[1] + a[2] + a[3];
		double r[2];
		r[0] = (rand()+1) / (RAND_MAX + 2.0); r[1] = (rand()+1) / (RAND_MAX + 2.0);
		double tau = -log(r[0]) / a_sum;

		// Sample identity of earliest reaction channel to fire (mu)
		int mu = 0; double s = a[0]; double r0 = r[1] * a_sum;

		while (s < r0)
		{
			mu += 1;
			s += a[mu];
		}

		if (rxn_count + 2 > MAX_OUTPUT_LENGTH)
		{
			printf("超出最大迭代次数");
			break;
		}

		result[rxn_count + 1][0] = result[rxn_count][0] + tau;
		for (int i = 1; i <= num_species; i++)
		{
			result[rxn_count + 1][i] = result[rxn_count][i] + stoich_matrix[mu][i - 1];
		}
		rxn_count += 1;
	}

	if (result[MAX_OUTPUT_LENGTH - 1][0] > tspan[1])
	{
		result[MAX_OUTPUT_LENGTH - 1][0] = tspan[1];
		for (int i = 1; i <= num_species; i++)
		{
			result[-1][i] = result[rxn_count - 1][i];
		}
	}

	ofstream outfile;
	outfile.open("resule.txt");
	for (int i = 0; i < rxn_count; i++)
	{

		outfile << result[i][0] << " " << result[i][1] << " " << result[i][2] << endl;

	}
	outfile.close();

}

int main(void)
{
	parm p;
	p.kr = 0.1; p.kp = 0.1; p.gr = 0.1; p.gp = 0.002;
	double tspan[2]; tspan[0] = 0;	tspan[1] = 10000;
	double x0[2];	 x0[0] = 0;	x0[1] = 0;

	double stoich_matrix[4][2];
	stoich_matrix[0][0] = 1; stoich_matrix[0][1] = 0; stoich_matrix[1][0] = 0; stoich_matrix[1][1] = 1;
	stoich_matrix[2][0] = -1; stoich_matrix[2][1] = 0; stoich_matrix[3][0] = 0; stoich_matrix[3][1] = -1;


	Gillespie(stoich_matrix, tspan, x0, p);



	return 0;
}