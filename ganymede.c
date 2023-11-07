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
	double f_tot;
	double f_conv;
	double f_core_tot;
	double f_core_conv;
	double base_ml;
	double base_mac;
	double v_ml;
	double v_mac;
	double R_ml;
	double R_mac;
	double lambda;
	//物理定数
	double G = 6.6743e-11;        //m^3 kg^-1 s^-2
	double mu = M_PI*4.0e-7;      //N A^-2
	//ガニメデ全体の定数
	double omega = 1.618e-6;	  //自転角速度[s^-1]
	double B = 7.2e-7;			  //赤道表面での磁場 [T = Wb m^-2 = kg A^-1 s^-2]
	//核の定数
	double r_core = 1.0e6;     	  //核の半径[m]
	double k_core = 5.0;          //核の熱伝導率
	double Cp_core = 8.0e2;       //核の定圧比熱[J kg^-1 K^-1]
	double rho_core = 6.0e3;  	  //核の密度[kg m^-3]
	double alpha_core = 8.0e-5;   //核の熱膨張率[K^-1]
	double sigma_core = 5.0e5;    //核の電気伝導率[S/m = kg^-1 m^-3 s^3 A^2]
	double a_core = 8.0e2;        // k/(ρ*Cp)
	//岩石層の定数
	double k_rock = 6.6;		  //岩石の熱伝導率（仮）
	double Cp_rock = 1.0e3;       //岩石の定圧比熱[J kg^-1 K^-1]
	double rho_rock = 3.3e3;      //岩石の密度[kg m^-3]
	double a_rock = 2.0e-6;       // k/(ρ*Cp)
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
		mass[i] = 4/3*M_PI * ( (rho_core-rho_rock)*pow(r_core,3.0) + rho_rock*pow(r[i],3.0) );
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
		newtemp[0] = 3 * f_cmb * dt / (Cp_rock * rho_core * r_core) + temp[0];				//内側の境界条件(∵ f_cmb(j-1) * S_c * dt = Cp_rock * rho_c * V_c * (temp[0](j-1)-temp[0](j)) )
		temp[77] = 278.6 ;			                                                	//外側の境界条件
		//温度の計算
		for(i=1; i<77; i++)
			newtemp[i] = dt/dr/dr*a_rock * ((1.0/(i+100)+1.0)*temp[i+1] - 2.0*temp[i] + (1.0-1.0/(i+100))*temp[i-1]) + temp[i] + q_all/Cp_rock;
			// newtemp[i] =  temp[i] + dt/(dr*dr) * k_rock/(Cp_rock*rho_rock) * ( (1.0+1.0/(i+100))*temp[i+1] - 2.0*temp[i] + (1.0-1.0/(i+100))*temp[i-1] ) + q_all/Cp_rock;
		//j+1での断熱温度勾配の計算
		f_ad = -k_core * alpha_core * newtemp[0] * g[0] / Cp_core;
		//j+1でのCMBでの熱流量の計算
		f_cmb = k_core * (newtemp[1]-newtemp[0]) / dr;
		//j+1でのコア内の半径rでの総熱流量
		f_core_tot = - rho_core * Cp_core * r_core * (newtemp[0] - temp[0])/dt /3;
		f_core_conv = f_core_tot - f_ad;
		//j+1での平均流速を計算
		base_ml = (r_core * f_core_conv * alpha_core * g[0]) / (rho_core * Cp_core);
		base_mac = (f_core_conv * alpha_core * g[0]) / (rho_core * Cp_core * omega);
		v_ml = 0.3 * pow( (r_core * f_core_conv * alpha_core * g[0]) / (rho_core * Cp_core), 1.0/3.0);
		v_mac = pow( (f_core_conv * alpha_core * g[0]) / (rho_core * Cp_core * omega), 1.0/2.0);
		// 磁気レイノルズ数
		R_ml = v_ml * r_core * mu * sigma_core;
		R_mac = v_mac * r_core * mu * sigma_core;
		// エルザッサー数
		lambda = sigma_core * pow(B,2) / ( 2 * rho_core * omega );
		//温度を更新
		for(i=0; i<77; i++)
			temp[i] = newtemp[i];
		//時間変化を出力
		printf("%d,%f,%f,%f,%f,%f,%f\n", j, temp[0], f_core_tot, f_core_conv, R_ml, R_mac, lambda);
		fprintf(file_j, "%d,%f,%f,%f,%f,%f,%f\n", j, temp[0], f_ad, f_cmb, R_ml, R_mac, lambda); 
	}
	fclose(file_j);


	//動径分布を出力
	FILE *file_i = fopen(filename_i, "w");
	for (int i=0; i<78; i++) {
		fprintf(file_i, "%d,%5.1f\n", i, temp[i]);
	}
	fclose(file_i);

}