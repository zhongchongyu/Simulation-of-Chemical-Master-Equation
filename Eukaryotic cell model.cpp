#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<fstream>
using namespace std;

#define k 10		//有k种反应发生
#define num_species 6		//有num_species中反应物
#define MAX_OUTPUT_LENGTH 10000	//最大迭代次数




class parm
{
public:
	double OPEN_dna, CLOSE_dna, ZL, OUT_rna, FY, IN_pro, D_rna_in, D_pro_in, D_rna_out, D_pro_out;
};

void Gillespie(double stoich_matrix[10][6], double tspan[2], double x0[6], parm p)
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

	int rxn_count = 0; double a[k];

	while (result[rxn_count][0] < tspan[1])
	{
		double DNA_open = result[rxn_count][1];
		double DNA_close = result[rxn_count][2];
		double RNA_in = result[rxn_count][3];
		double pro_in = result[rxn_count][4];
		double RNA_out = result[rxn_count][5];
		double pro_out = result[rxn_count][6];

		a[0] = p.OPEN_dna * DNA_close;
		a[1] = p.CLOSE_dna * DNA_open;
		a[2] = p.ZL * DNA_open;
		a[3] = p.OUT_rna * RNA_in;
		a[4] = p.FY * RNA_out;
		a[5] = p.IN_pro * pro_out;
		a[6] = p.D_rna_in * RNA_in;
		a[7] = p.D_pro_in * pro_in;
		a[8] = p.D_rna_out * RNA_out;
		a[9] = p.D_pro_out * pro_out;

		double a_sum = a[0] + a[1] + a[2] + a[3] + a[4] + a[5] + a[6] + a[7] + a[8] + a[9];
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
	outfile.open("result.txt");
	for (int i = 0; i < rxn_count; i++)
	{
		for (int j = 0; j <= num_species; j++)
		{
			outfile << result[i][j] << " ";
		}
		outfile << endl;

		// outfile << result[i][0] << " " << result[i][1] << " " << result[i][2] << endl;
	}
	outfile.close();

}

int main(void)
{
	parm p;
	p.OPEN_dna = 0.1; p.CLOSE_dna = 0.8; p.ZL = 0.3; p.OUT_rna = 0.7; p.FY = 0.5; p.IN_pro = 0.005;
	p.D_rna_in = 0.02; p.D_pro_in = 0.02; p.D_rna_out = 0.02; p.D_pro_out = 0.02;

	double tspan[2]; tspan[0] = 0;	tspan[1] = 10000;
	double x0[6];	 x0[0] = 10;	x0[1] = 0; x0[2] = 0; x0[3] = 0; x0[4] = 0; x0[5] = 0;

	double stoich_matrix[10][6];
	stoich_matrix[0][0] = 1; stoich_matrix[0][1] = -1; stoich_matrix[0][2] = 0; stoich_matrix[0][3] = 0; stoich_matrix[0][4] = 0; stoich_matrix[0][5] = 0;		//DNA开
	stoich_matrix[1][0] = -1; stoich_matrix[1][1] = 1; stoich_matrix[1][2] = 0; stoich_matrix[1][3] = 0; stoich_matrix[1][4] = 0; stoich_matrix[1][5] = 0; 		//DNA关
	stoich_matrix[2][0] = 0; stoich_matrix[2][1] = 0; stoich_matrix[2][2] = 1; stoich_matrix[2][3] = 0; stoich_matrix[2][4] = 0; stoich_matrix[2][5] = 0; 		//转录
	stoich_matrix[3][0] = 0; stoich_matrix[3][1] = 0; stoich_matrix[3][2] = -1; stoich_matrix[3][3] = 0; stoich_matrix[3][4] = 1; stoich_matrix[3][5] = 0;		//RNA跳出
	stoich_matrix[4][0] = 0; stoich_matrix[4][1] = 0; stoich_matrix[4][2] = 0; stoich_matrix[4][3] = 0; stoich_matrix[4][4] = 0; stoich_matrix[4][5] = 1;		//翻译
	stoich_matrix[5][0] = 0; stoich_matrix[5][1] = 0; stoich_matrix[5][2] = 0; stoich_matrix[5][3] = 1; stoich_matrix[5][4] = 0; stoich_matrix[5][5] = -1;		//pro进入
	stoich_matrix[6][0] = 0; stoich_matrix[6][1] = 0; stoich_matrix[6][2] = -1; stoich_matrix[6][3] = 0; stoich_matrix[6][4] = 0; stoich_matrix[6][5] = 0;		//核RNA衰变
	stoich_matrix[7][0] = 0; stoich_matrix[7][1] = 0; stoich_matrix[7][2] = 0; stoich_matrix[7][3] = -1; stoich_matrix[7][4] = 0; stoich_matrix[7][5] = 0;		//核pro衰变
	stoich_matrix[8][0] = 0; stoich_matrix[8][1] = 0; stoich_matrix[8][2] = 0; stoich_matrix[8][3] = 0; stoich_matrix[8][4] = -1; stoich_matrix[8][5] = 0;		//质RNA衰变
	stoich_matrix[9][0] = 0; stoich_matrix[9][1] = 0; stoich_matrix[9][2] = 0; stoich_matrix[9][3] = 0; stoich_matrix[9][4] = 0; stoich_matrix[9][5] = -1;		//质pro衰变	

	Gillespie(stoich_matrix, tspan, x0, p);



	return 0;
}