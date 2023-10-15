#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <complex.h>
#include <math.h>


int main( void ){ 

	//変数の宣言：計算式中で使用する文字はここで宣言しておく
	int 	i, j;				// 整数
	double	temp[78], newtemp[78], dt, dr, a, Cp, h_k, h_th, h_u235, h_u238, c_k, c_th, c_u235, c_u238, t_half_k, t_half_th, t_half_u235, t_half_u238, q_k, q_th, q_u235, q_u238, q_all ;
	FILE *fp;				// 計算結果をファイルに出力するときに必要な特殊な変数
	dt = 3.1536e+11;   //10000年
	dr = 10000.0;        //単位m
	a = 2.0e-6;
	Cp = 1.0e3; //定圧比熱[J/kg/K]
	h_k = 29.17e-6; //単位質量当たりの発熱量[W/kg]
	h_th = 26.38e-6;
	h_u235 = 568.7e-6;
	h_u238 = 94.65e-6;
	c_k = 1104e-9; //濃集率[]
	c_th = 53.8e-9; 
	c_u235 = 8.2e-9;
	c_u238 = 26.2e-9;
	t_half_k = 1277; //半減期[Myr]
	t_half_th = 14030;
	t_half_u235 = 703.81;
	t_half_u238 = 4468;

	// ファイル名を日時で生成する
	time_t t = time(NULL);
	struct tm *now = localtime(&t);
	char filename_i[64], filename_j[64];
	strftime(filename_i, sizeof(filename_i), "ganymede_%m%d%H%M_i.csv", now);
	strftime(filename_j, sizeof(filename_j), "ganymede_%m%d%H%M_j.csv", now);

	//初期条件
	for(i=0; i<78; i++)		
		temp[i] = 300;	// その他の部分(temp[1]からtemp[99])の初期温度


	FILE *file_j = fopen(filename_j, "w");
	// 熱伝導方程式の計算（for文を使った繰り返し）
	for(j = 0; j < 460000; j++){	//時間ステップの繰り返し
		q_k = c_k * (h_k*dt) * exp(-log(2.0) * j / (t_half_k*100.0));
        q_th = c_th * (h_th*dt) * exp(-log(2.0)*j/(t_half_th*100.0));
        q_u235 = c_u235 * (h_u235*dt) * exp(-log(2.0)*j/(t_half_u235*100.0));
        q_u238 = c_u238 * (h_u238*dt) * exp(-log(2.0)*j/(t_half_u238*100.0));
        q_all = q_k + q_th + q_u235 + q_u238;

		for(i=1; i<77; i++)
			newtemp[i] = dt/dr/dr*a * ((1/(i+100)+1.0)*temp[i+1] - 2.0*temp[i] + (1.0-1/(i+100))*temp[i-1]) + temp[i] + q_all/Cp; 	// 差分方程式

		for(i=1; i<77; i++)
			temp[i] = newtemp[i];				// 温度を次の時間ステップのものに更新する
		
		//境界条件
		temp[0] = temp[1];					// 熱浴側の境界条件
		temp[77] = 278.6 ;				// 熱浴とは反対側の境界条件（太陽からのエネルギーと宇宙への放射が平衡）

		//温度の時間変化を出力
		printf("%d,%5.1f\n", j, temp[0]);
		fprintf(file_j, "%d,%5.1f\n", j, temp[0]); 
	}				// 時間ステップ繰り返し終了
	fclose(file_j);


	//温度の動径分布を出力
	FILE *file_i = fopen(filename_i, "w");
	for (int i=0; i<77; i++) {
		fprintf(file_i, "%d,%5.1f\n", i, temp[i]);
	}
	fclose(file_i);

}