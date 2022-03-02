#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define n 500//天体個数
#define seed 149//乱数生成用seed
#define nt 30//時間ステップ数上限
#define arraySize 100//格子数
#define G 1//万有引力定数
#define m 1//質量
#define H 3.5//ハッブル定数
#define ni 50//ガウス・ザイデル法循環上限
#define dt 0.1//時間ステップ
#define dx 1.0//x方向微小量
#define dy 1.0//y方向微小量

//乱数の発生用関数
double creat_ran(double a, double b){
	double x = a + ((b-a) * rand() / (RAND_MAX + 1.0));
	return x;
}

//データファイル出力用関数
void creat_data_file(char *fname, double x[], double y[]){
    FILE *fo;
	if ((fo=fopen(fname, "w")) == NULL){
		printf("File[%s] dose not open !!\n", fname);
		exit(0);
	}
	for(int i=0; i<n; i++) fprintf(fo, "%f,%f\n", x[i], y[i]);
	fclose(fo);
}


void main(void){
    //座標
    double x[n];
	double y[n];
    
    //速度
    double vx[n];
    double vy[n];

    double phi[arraySize][arraySize];//ポテンシャル
    double ro[arraySize][arraySize];//質量密度
    double Fx[arraySize][arraySize];//横方向重力場
    double Fy[arraySize][arraySize];//縦方向重力場

	srand(seed);//乱数seed設定

    //初期座標および初期速度の設定
	for(int j=0; j<n; j++){
		double a = creat_ran(-5.0, 5.0);
		double b = creat_ran(-sqrt(25.0 - pow(a, 2.0)), sqrt(25.0 - pow(a, 2.0)));
		x[j] = a;
		y[j] = b;
        vx[j] = H * x[j];
        vy[j] = H * y[j];
	}

    creat_data_file("data0.csv", x, y);//初期座標データを取得

    //演算ループ開始
    for(int p=1; p<=nt; p++){
        //初期値0で初期化する
        for(int i=0; i<arraySize; i++){
            for(int j=0; j<arraySize; j++){
                phi[i][j] = 0.0;
                ro[i][j] = 0.0;
                Fx[i][j] = 0.0;
                Fy[i][j] = 0.0;
            }
        }
        //NGP法で質量密度roの更新
        for(int j=0; j<n; j++){
            if(fmod(x[j],1.0) >= 0.5){
                if(fmod(y[j], 1.0) >= 0.5)ro[(int)ceil(x[j])][(int)ceil(y[j])] = ro[(int)ceil(x[j])][(int)ceil(y[j])] + m;
                else ro[(int)ceil(x[j])][(int)floor(y[j])] = ro[(int)ceil(x[j])][(int)floor(y[j])] + m;
            }else{
                if(fmod(y[j], 1.0) >= 0.5)ro[(int)floor(x[j])][(int)ceil(y[j])] = ro[(int)floor(x[j])][(int)ceil(y[j])] + m;
                else ro[(int)floor(x[j])][(int)floor(y[j])] = ro[(int)floor(x[j])][(int)floor(y[j])] + m;
            }
        }
        //ガウス・ザイデル法でポテンシャルphiの更新
        for(int k=0; k<ni; k++){
            for(int i=1; i<arraySize-1; i++){
                for(int t=1; t<arraySize-1; t++){
                    phi[i][t] = (phi[i+1][t] + phi[i-1][t] + phi[i][t-1] + phi[i][t+1] - dx*dx*G*ro[i][t])/4;//N+1回目の数値をガウス・ザイデル法で計算      
                }
            }
        }
        //重力場Fx, Fyの更新
        for(int i=1; i<arraySize-1; i++){
            for(int t=1; t<arraySize-1; t++){
                Fx[i][t] = -(phi[i+1][t]-phi[i][t])/dx;
                Fy[i][t] = -(phi[i][t+1]-phi[i][t])/dy;
            }
        }

        for(int i=0; i<n; i++){
            double Fpx = 0.0;//粒子x方向受力
            double Fpy = 0.0;//粒子y方向受力
            //粒子受力の更新    
            if(fmod(x[i], 1.0) >= 0.5){
                if(fmod(y[i], 1.0) >= 0.5){
                    Fpx = m*Fx[(int)ceil(x[i])][(int)ceil(y[i])];
                    Fpy = m*Fy[(int)ceil(x[i])][(int)ceil(y[i])];
                }else{
                    Fpx = m*Fx[(int)ceil(x[i])][(int)floor(y[i])];
                    Fpy = m*Fy[(int)ceil(x[i])][(int)floor(y[i])];
                }
            }else{
                if(fmod(y[i], 1.0) >= 0.5){
                    Fpx = m*Fx[(int)floor(x[i])][(int)ceil(y[i])];
                    Fpy = m*Fy[(int)floor(x[i])][(int)ceil(y[i])];
                }else{
                    Fpx = m*Fx[(int)floor(x[i])][(int)floor(y[i])];
                    Fpy = m*Fy[(int)floor(x[i])][(int)floor(y[i])];
                }
            }
        //粒子速度の更新
            vx[i] = vx[i] + (Fpx/m)*dt;
            vy[i] = vy[i] + (Fpy/m)*dt;
        //粒子位置の更新
            x[i] = x[i] + vx[i]*dt;
            y[i] = y[i] + vy[i]*dt;
        }
        //出力するデータ
        if(p==10)creat_data_file("data10.csv", x, y);
        if(p==20)creat_data_file("data20.csv", x, y);
        if(p==30)creat_data_file("data30.csv", x, y);
    }    
}
