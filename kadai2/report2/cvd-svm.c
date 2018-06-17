#include <iostream>
#include <sstream>
#include <string>
#include "QuadProg++.hh"
#include <stdio.h>
#include <string.h>
#include <math.h>

#define SCALE 1

int getmatrix(char *fname, int *dim, double vector[MATRIX_DIM][MATRIX_DIM]);
double calc(int i, int j, int dim, double train[MATRIX_DIM][MATRIX_DIM], int num, double SIGMA);
int maxalpha(double alpha[], int size);

int main(int argc, char *argv[])
{
  double G[MATRIX_DIM][MATRIX_DIM], g0[MATRIX_DIM],
		CE[MATRIX_DIM][MATRIX_DIM], ce0[MATRIX_DIM],
		CI[MATRIX_DIM][MATRIX_DIM], ci0[MATRIX_DIM],
		alpha[MATRIX_DIM];
	int n, m, p;
  int dim; //ベクトル次元数
  int split = atoi(argv[3]); //分割数
  FILE *gp, *fp;
  double Accuracy = 0;
  double SIGMA = 0;
  double bestSIGMA = 0;
  double bestAccuracy = 0;
  double bests = 0;

  if(atoi(argv[1]) != 0 && atoi(argv[1]) != 1 && atoi(argv[1]) != 2){
    printf("第一引数に番号を指定\n");
    printf("0 : カーネルトリック無し\n");
    printf("1 : 多項式カーネル\n");
    printf("2 : Gaussカーネル\n");
    exit(1);
  }

  if(argv[2] == NULL){
    printf("第二引数にファイルを指定\n");
    exit(1);
  }

  if(argv[3] == NULL){
    printf("第三引数に分割数を入力\n");
    exit(1);
  }


  double vector[MATRIX_DIM][MATRIX_DIM]; //特徴ベクトル

  n = getmatrix(argv[2], &dim, vector); //行列サイズを取得、同時に特徴ベクトルを配列に整理する

  for(int i = 0; i < n; i++){  //スケーリング
    for(int a = 0; a < dim; a++){
      vector[i][a] = vector[i][a] / SCALE;
    }
  }

  if((n % split) != 0){
    printf("データ数%dが分割数%dで割り切れません\n", n, split);
    exit(1);
  }


  int testnum = n / split; //評価データ数
  int tn = n - testnum; //訓練データ数


fp = fopen("cvresult.dat", "w"); //結果記録用

//SIGMA探索ループ開始////////////////////////////////////////////////////////////



  for(double s = 9; s < 13; s = s + 0.1){  //最適なパラメータSIGMAを探索
    if(atoi(argv[1]) == 0 || atoi(argv[1]) == 1){ //ガウスカーネル以外は再帰しない
      s = 100;
    }

    SIGMA = pow(2, s);


//交差検定ループ開始//////////////////////////////////////////////////////////////



  for(int cvd = 0; cvd < split; cvd++){  //分割数を元に交差検定を行う
    printf("テスト%d回目\t", cvd + 1);

    double train[MATRIX_DIM][MATRIX_DIM];  //訓練データ
    double test[MATRIX_DIM][MATRIX_DIM]; //評価データ

    for(int i = 0; i < n; i++){  //取得データを分割用の配列に入れる
      for(int a = 0; a < dim + 1; a++){
        train[i][a] = vector[i][a];
      }
    }


    for(int i = 0; i < testnum; i++){ //評価データを代入
      for(int a = 0; a < dim + 1; a++){
        test[i][a] = train[(cvd * testnum) + i][a];
      }
    }

    for(int i = 0; (cvd * testnum) + i < tn; i++){ //訓練データを代入
      for(int a = 0; a < dim + 1; a++){
        train[(cvd * testnum) + i][a] = train[(cvd * testnum) + i + testnum][a];
      }
    }


      for(int i = 0; i < tn; i++){  //行列Gの要素を計算して代入
        for(int j = 0; j < tn; j++){
          G[i][j] = train[i][dim] * train[j][dim] * calc(i, j, dim, train, atoi(argv[1]), SIGMA);
          if(i == j){
            if(atoi(argv[1]) == 1){
              G[i][j] += 1.0e-7;
            }else{
              G[i][j] += 1.0e-9;
            }
          }
        }
      }

      for(int i=0; i<tn; i++){
        g0[i] = -1.0;
      }

m = tn;
      for(int i=0; i<tn; i++){
        CE[i][0] = train[i][dim];
      }

      ce0[0] = 0.0;

p = 1;
      for(int i=0; i<tn; i++){
        for(int j=0; j<tn; j++){
          if(i == j){
            CI[i][j] = 1.0;
          }else{
            CI[i][j] = 0.0;
          }
        }
      }

      for(int i=0; i<tn; i++){
        ci0[i] = 0.0;
      }


solve_quadprog(G, g0, tn, CE, ce0, p, CI, ci0, m, alpha); //２次計画問題を解く


double w[MATRIX_DIM] = {0};
for(int a = 0; a < dim; a++){ //重み計算
  for(int i = 0; i < tn; i++){
    w[a] += alpha[i] * train[i][dim] * train[i][a];
  }
}

double theta = 0.0; //閾値
int k =  maxalpha(alpha, tn);//alphaが最大になる次数

if(atoi(argv[1]) == 0){ // カーネルトリックなしのtheta
  for(int a = 0; a < dim; a++){
    theta += w[a] * train[k][a];
  }
  theta = theta - train[k][dim];
}

if(atoi(argv[1]) == 1){ // 多項式カーネルのtheta
  double inner = 0;
  for(int i = 0; i < tn; i++){
    for(int a = 0; a < dim; a++){
      inner += train[i][a] * train[k][a]; //内積を計算
    }
    theta += alpha[i] * train[i][dim] * pow(1.0 + inner, dim + 1);
    inner = 0;
  }
}

if(atoi(argv[1]) == 2){ // ガウスカーネルのtheta
  double norm = 0.0;
  for(int i = 0; i < tn; i++){
    for(int a = 0; a < dim; a++){
      norm += pow(train[i][a] - train[k][a], 2);
    }
    theta += alpha[i] * train[i][dim] * exp(-norm / (2 * pow(SIGMA, 2)));
    norm = 0;
  }
}



//カーネルごとに識別開始///////////////////////////////////////////////////////////



if(atoi(argv[1]) == 0){ //カーネルトリックなしの識別器
  double a, b; // y = ax + b
  double inner = 0;
  double correct = 0; //正しく識別できた数をカウント

  for(int i = 0; i < testnum; i++){ //評価データが正しく識別されているか調べる
    for(int a = 0; a < dim; a++){
      inner += w[a] * test[i][a]; //内積を計算
    }

    if((inner - theta) > 0){
      if(test[i][dim] == 1){
        correct += 1;
      }
    }

    if((inner - theta) < 0){
      if(test[i][dim] == -1){
        correct += 1;
      }
    }
inner = 0;
}
printf("正答数 : %d / %d\t\tAccuracy : %lf%%\n", int(correct), testnum, (correct / testnum) * 100);
Accuracy += correct / testnum;



} //カーネルトリックなしの識別終わり



if(atoi(argv[1]) == 1 || atoi(argv[1]) == 2){ //多項式カーネルorガウスカーネル
  double x[MATRIX_DIM] = {0}; //xベクトル
  double kernel = 0;
  double inner = 0;
  double norm = 0.0;
  double f = 0;
  double check; //識別関数
  double correct = 0;



  for(int i = 0; i < testnum; i++){  //評価データが正しく識別されているか調べる
    for(int k = 0; k < tn; k++){

      if(atoi(argv[1]) == 1){ //多項式カーネルの場合
        for(int a = 0; a < dim; a++){
          inner += train[k][a] * test[i][a];
        }
        kernel = pow(1.0 + inner, dim + 1);
      }

        if(atoi(argv[1]) == 2){ //ガウスカーネルの場合
          for(int a = 0; a < dim; a++){
            norm += pow(train[k][a] - test[i][a], 2);
          }
          kernel = exp(-norm / (2 * pow(SIGMA, 2)));
        }

    f += alpha[k] * train[k][dim] * kernel;
    norm = 0;
    inner = 0;
    kernel = 0;
    }
    check = f - theta;

    if(check > 0){
      if(test[i][dim] == 1){
        correct += 1;
      }
    }

    if(check < 0){
      if(test[i][dim] == -1){
        correct += 1;
      }
    }
    f = 0;
  }
  printf("正答数 : %d / %d\t\tAccuracy : %lf%%\n", int(correct), testnum, (correct / testnum) * 100);
  Accuracy += correct / testnum;



} //多項式カーネルorガウスカーネルの識別終わり


} //交差検定を繰り返す



if(atoi(argv[1]) == 2){
  printf("σ : %lf(s = %lf)\t\tAccuracy : %lf%%\n", SIGMA, s, (Accuracy * 100) / split);
}else{
  printf("Accuracy : %lf%%\n", (Accuracy * 100) / split);
}

if(Accuracy > bestAccuracy){  //Accuracyが最大になる場所を探す
  bestAccuracy = Accuracy;
  bestSIGMA = SIGMA;
  bests = s;
}


fprintf(fp, "%lf\t%lf\n", s, (Accuracy * 100) / split); //結果をファイルに記録


Accuracy = 0;
} //SIGMA探索を繰り返す

if(atoi(argv[1]) == 2){
  printf("bestσ : %lf(s = %lf)\t\tbestAccuracy : %lf%%\n", bestSIGMA, bests, (bestAccuracy * 100) / split);
}

fclose(fp);
  return 1;
}





int getmatrix(char *fname, int *dim, double vector[MATRIX_DIM][MATRIX_DIM]) //Gの行列サイズを取得、同時に特徴ベクトルを配列に整理する
{
  FILE *fp;
  char s[100];
  int count = 0;
  int a = 0, b = 0;
  char *p;

  fp = fopen(fname, "r" );
  if(fp == NULL){
    printf("%sファイルが開けません\n", fname);
    exit(1);
  }

  while(fgets(s, 100, fp) != NULL){
    vector[a][b] = atof(strtok(s, " "));
    b += 1; //ベクトルの何番目の要素かを表す

    do{
        p = strtok(NULL, " ");
        if(p){
        vector[a][b] = atof(p);
        b += 1;
      }
    }while(p);

    a += 1; //次のベクトルへ
    *dim = b - 1; //ベクトル次元数を保存
    b = 0;
    count += 1; //行列サイズ
  }

  fclose(fp);
  return count;
}





double calc(int i, int j, int dim, double train[MATRIX_DIM][MATRIX_DIM], int num, double SIGMA) //カーネルを計算
{
  double inner = 0;

  for(int a = 0; a < dim; a++){
    inner += train[i][a] * train[j][a]; //内積を計算
  }

if(num == 0){
  return inner; // (x_i, x_j)
}else if(num == 1){
  return pow(1.0 + inner, dim + 1); //多項式カーネル
}else if(num == 2){ //ガウスカーネル
  double norm = 0.0;

  for(int a = 0; a < dim; a++){
    norm += pow(train[i][a] - train[j][a], 2);
  }
  return exp(-norm / (2 * pow(SIGMA, 2)));
}else{
  printf("カーネル計算エラー\n");
  return 1.0;
}
}





int maxalpha(double alpha[], int size) //alphaが最大になる次数を求める
{
  int i;
  int k = 0;

for(i = 0; i < size; i++){
  if(alpha[i] > alpha[k]){
    k = i; //alphaが最大になる次数を格納
  }
}

  return k;
}
