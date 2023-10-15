#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <complex.h>
#include <math.h>


int main( void ){ 				// おまじない（必要） プログラムの本体



///////////変数の宣言：計算式中で使用する文字はここで宣言しておく/////////
	int 	i, j;				// 整数
	double	temp[78], newtemp[78], dt, dr, a, Cp, h_al, h_fe, h_mn, h_k, c_al, c_fe, c_mn, c_k, t_half_al, t_half_fe, t_half_mn, t_half_k, q_al, q_fe, q_mn, q_k, q_all ;	// 実数（temp[0]からtemp[99]の合計100個の配列）
	FILE *fp;				// 計算結果をファイルに出力するときに必要な特殊な変数
	dt = 3.1536e+11;   //10000年
	dr = 10000;        //単位m
	a = 2.0e-6;
	Cp = 1.0e3;
	h_al = 0.146;
	h_fe = 0.070;
	h_mn = 0.027;
	h_k = 29.17e-6;
	c_al = 87.5e-9;
	c_fe = 39.7e-9;
	c_mn = 17.7e-9;
	c_k = 1104;
	t_half_al = 0.72e6; 
	t_half_fe = 1.5e6;
	t_half_mn = 3.7e6;
	t_half_k = 1277e6;

///////////// 初期条件の設定 ////////////////
	for(i=0; i<78; i++)		// for:繰り返し　i=1:iが1からスタート　i<100:iが100より小さい間　i++:iを1つずつ追加
		temp[i] = 500;	// その他の部分(temp[1]からtemp[99])の初期温度

/////////////////////////////////////////////



/////// 熱伝導方程式の計算（for文を使った繰り返し）//////

	for(j = 0; j < 460000; j++){	//時間ステップの繰り返し

        // q_al = c_al * h_al * exp(-log(2)/t_half_al*j*100);
        // q_fe = c_fe * h_fe * exp(-log(2)/t_half_fe*j*100);
        // q_mn = c_mn * h_mn * exp(-log(2)/t_half_mn*j*100);
        q_k = c_k * h_k * exp(-log(2)/t_half_k*j*100);
        // q_all = q_al + q_fe + q_mn + q_k;

		for(i=1; i<77; i++)
			newtemp[i] = dt/dr/dr*a * ((1/(i+100)+1)*temp[i+1] - 2*temp[i] + (1-1/(i+100))*temp[i-1]) + temp[i] + q_k/Cp; 	// 差分方程式

		for(i=1; i<77; i++)
			temp[i] = newtemp[i];				// 温度を次の時間ステップのものに更新する
		
		//境界条件
		temp[0] = temp[1];					// 熱浴側の境界条件
		temp[77] = 278.6 ;				// 熱浴とは反対側の境界条件（太陽からのエネルギーと宇宙への放射が平衡）

		printf("%d,%f\n", j, temp[0]);

	}				// 時間ステップ繰り返し終了

///////////////////////////////////////////////////////////



	//////////////////出力(分布)///////////////////////////////
	// ファイル名を日時で生成する
	time_t t = time(NULL);
	struct tm *now = localtime(&t);
	char filename[64];
	strftime(filename, sizeof(filename), "ganymede_%m%d%H%M.csv", now);

	// ファイルを書き込みモードで新規に作成し、開く
	FILE *file = fopen(filename, "w");

	if (file == NULL) {
		printf("ファイルを作成または開くことができませんでした。\n");
		return 1;
	}

	for (int i=0; i<77; i++) {
		fprintf(file, "%d,%5.1f\n", i, temp[i]);
	}

	// ファイルを閉じる
	fclose(file);
	////////////////////////////////////////////////////////////



}				// main( void )の終わり（プログラムの終わり）



//VSCode閉じてコンパイルしないと計算結果がバグる