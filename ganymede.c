#include <stdio.h>		// おまじない（必要)
#include <stdlib.h>		// おまじない（必要）
#include <time.h>

// 「//」以降はコメント　プログラムには反映されない


int main( void ){ 				// おまじない（必要） プログラムの本体


///////////変数の宣言：計算式中で使用する文字はここで宣言しておく/////////

	int 	i, j;				// 整数
	double	temp[128], newtemp[128], dt, dr, a;	// 実数（temp[0]からtemp[99]の合計100個の配列）
	FILE *fp;				// 計算結果をファイルに出力するときに必要な特殊な変数
	dt = 3.1536e+11;   //10000年
	dr = 50000;        //50km*127=6350km
	a = 2.0e-6;

////////// 計算結果を書き込むファイル(例 "result.txt")を作成//////////////

	if ((fp=fopen("result.txt","w")) == NULL){
		printf("Cannot open the file\n");
		exit(1);
	}


///////////// 初期条件の設定 ////////////////

	for(i = 0; i<127; i++)		// for:繰り返し　i=1:iが1からスタート　i<100:iが100より小さい間　i++:iを1つずつ追加して 
		temp[i] = 2000;	// その他の部分(temp[1]からtemp[99])の初期温度

/////////////////////////////////////////////



/////// 熱伝導方程式の計算（for文を使った繰り返し）//////

	for(j = 0; j < 460000; j++){	//時間ステップの繰り返し j<○○ の○○の値を増やせば計算時間も増やせる

		for(i = 1; i<127; i++)
			newtemp[i] = dt/dr/dr*a * ( (1/i + 1)*temp[i+1] - 2*temp[i] + (1 - 1/i)*temp[i-1])+temp[i]; 	// 差分方程式

		for(i = 1; i<127; i++)
			temp[i] = newtemp[i];				// 温度を次の時間ステップのものに更新する
		
		//境界条件
		temp[0] = temp[1];					// 熱浴側の境界条件
		temp[127] = 278.6 ;				// 熱浴とは反対側の境界条件（太陽からのエネルギーと宇宙への放射が平衡）

		//出力部分1
		// printf("%d %5.1f\n", j, temp[20]);
		// fprintf(fp, "%d %5.1f\n", j, temp[20]);

	}				// 時間ステップ繰り返し終了
	
	//出力部分2
	for(i = 0; i < 127; i++){
		printf("%d %5.1f\n", i, temp[i]);			//画面にjとtemp[10]を出力 %d→整数の表示　%f→小数点つき数の表示（5.1は全部で5桁，小数点以下1桁の意味）　\n→改行
		fprintf(fp, "%d %5.1f\n", i, temp[i]);			//result.txtに出力
	}

	//csvファイル作成
	// ファイル名を日時で生成する
    time_t t = time(NULL);
    struct tm* now = localtime(&t);
    char filename[64];
    strftime(filename, sizeof(filename), "ganymede_%Y%m%d_%H%M%S.csv", now);

    // ファイルを書き込みモードで新規に作成し、開く
    FILE* file = fopen(filename, "w");

    if (file == NULL) {
        printf("ファイルを作成または開くことができませんでした。\n");
        return 1;
    }

    for (int i = 0; i < 127; i++) {
        fprintf(file, "%d,%5.1f\n", i, temp[i]);
    }

    // ファイルを閉じる
    fclose(file);

////////////////////////////////////////////////////////

	fclose(fp);		// 開いたresult.txtを閉じる。
}				// main( void )の終わり（プログラムの終わり）



//VSCode閉じてコンパイルしないと計算結果がバグる