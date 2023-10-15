#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <complex.h>
#include <math.h>

/*memo
・ファイル出力先の変更
・半減期の単位をyrに合わせる
・岩石の定数の整理
・r,tのグリッド数を自動で決定するコードを追加する
*/


int main( void ){ 
	//ファイル出力用
	FILE *fp;
	time_t t = time(NULL);
	struct tm *now = localtime(&t);
	char filename_i[64], filename_j[64], filename_g[64];
	strftime(filename_i, sizeof(filename_i), "ganymede_%m%d%H%M_i.csv", now);
	strftime(filename_j, sizeof(filename_j), "ganymede_%m%d%H%M_j.csv", now);
	strftime(filename_g, sizeof(filename_g), "ganymede_%m%d%H%M_g.csv", now);
	//ループ変数
	int i, j;
	//計算の幅
	double dt = 3.1536e+11;   //s (=10000yr)
	double dr = 10000.0;      //m
	//計算用の変数
	double	temp[78];
	double newtemp[78];
	double r[78];
	double mass[78];
	double g[78];
	double f_ad;
	double f_cmb;
	//物理定数
	double G = 6.6743e-11;     //m^3 kg^-1 s^-2
	//半径,密度
	double r_core = 1.0e6;     //m
	double dens_core = 6.0e3;  //kg m^-3
	double dens_rock = 3.3e3;  //kg m^-3
	//その他定数
	double a = 2.0e-6;         //
	double a_c = 8.0e2;        //J kg^-1 K^-1
	double alpha_c = 8.0e-5;   //K^-1
	double k_c = 5.0;          //
	double Cp_c = 1.0e3;       //J kg^-1 K^-1 (定圧比熱)
	//単位質量当たりの発熱量[W/kg]
	double h_k = 29.17e-6;
	double h_th = 26.38e-6;
	double h_u235 = 568.7e-6;
	double h_u238 = 94.65e-6;
	//濃集率[]
	double c_k = 1104e-9; 
	double c_th = 53.8e-9; 
	double c_u235 = 8.2e-9;
	double c_u238 = 26.2e-9;
	//半減期[Myr]
	double t_half_k = 1277; 
	double t_half_th = 14030;
	double t_half_u235 = 703.81;
	double t_half_u238 = 4468;
	//合計発熱量
	double q_k;
	double q_th;
	double q_u235;
	double q_u238;
	double q_all;


	//動径ステップの繰り返し
	////FILE *file_g = fopen(filename_g, "w");
	for(i=0; i<78; i++)	{
		//初期条件の決定
		temp[i] = 1100;
		//重力加速度の計算
		r[i] = r_core + i * dr;
		mass[i] = 4/3*M_PI * ( (dens_core-dens_rock)*pow(r_core,3.0) + dens_rock*pow(r[i],3.0) );
		g[i] = G * mass[i] / pow(r[i],2.0);
	}
	////fclose(file_g);


	//時間ステップの繰り返し
	FILE *file_j = fopen(filename_j, "w");
	for(j = 0; j < 460001; j++){
		//発熱量の計算
		q_k = c_k * (h_k*dt) * exp(-log(2.0) * j / (t_half_k*100.0));
        q_th = c_th * (h_th*dt) * exp(-log(2.0)*j/(t_half_th*100.0));
        q_u235 = c_u235 * (h_u235*dt) * exp(-log(2.0)*j/(t_half_u235*100.0));
        q_u238 = c_u238 * (h_u238*dt) * exp(-log(2.0)*j/(t_half_u238*100.0));
        q_all = q_k + q_th + q_u235 + q_u238;
		//境界条件を適用
		temp[0] = 3 * f_cmb * dt / (Cp_c * dens_core * r_core) + temp[0];				//内側の境界条件(∵ f_cmb(j-1) * S_c * dt = Cp_c * dens_c * V_c * (temp[0](j-1)-temp[0](j)) )
		temp[77] = 278.6 ;			                                                //外側の境界条件
		//温度の計算
		for(i=1; i<77; i++)
			newtemp[i] = dt/dr/dr*a * ((1.0/(i+100)+1.0)*temp[i+1] - 2.0*temp[i] + (1.0-1.0/(i+100))*temp[i-1]) + temp[i] + q_all/Cp_c;
		// printf("%d,%f+%f=%f\n", j, 3 * f_cmb * dt / (Cp_c * dens_core *r_core), temp[0], newtemp[0]);
		printf("%d,%f,%f,%f,%f\n", j, temp[0], newtemp[0], f_ad, f_cmb);
		//温度を更新
		for(i=1; i<77; i++)
			temp[i] = newtemp[i];
		//断熱温度勾配の計算
		f_ad = k_c * alpha_c * temp[0] * g[0] / Cp_c;
		//熱流量の計算
		f_cmb = k_c * (temp[1]-temp[0]) / dr;
		//時間変化を出力

		// printf("%d,%5f,%f\n", j, temp[0] );
		fprintf(file_j, "%d,%f,%f,%f\n", j, temp[0], f_ad, f_cmb); 
	}
	fclose(file_j);



	//動径分布を出力
	FILE *file_i = fopen(filename_i, "w");
	for (int i=0; i<78; i++) {
		fprintf(file_i, "%d,%5.1f\n", i, temp[i]);
	}
	fclose(file_i);

}