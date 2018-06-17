#include <iostream>
#include <sstream>
#include <string>
#include "QuadProg++.hh"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


#define EPSILON 0.1
#define VECTOR_DIM 10
#define TIME_DIM 6


int getmatrix(char *fname, double vector[MATRIX_DIM][VECTOR_DIM], int timestamp[MATRIX_DIM][TIME_DIM], int dim, double day2[MATRIX_DIM][VECTOR_DIM], int *prednum);
double calc(int i, int j, int dim, double vector[MATRIX_DIM][VECTOR_DIM], int num, double SIGMA);


int main(int argc, char *argv[])
{
  double G[MATRIX_DIM][MATRIX_DIM], g0[MATRIX_DIM],
		CE[MATRIX_DIM][MATRIX_DIM], ce0[MATRIX_DIM],
		CI[MATRIX_DIM][MATRIX_DIM], ci0[MATRIX_DIM],
		alpha[MATRIX_DIM];
  double SIGMA;
  double C;
	int n, m, p;
  int dim = atoi(argv[3]); //ベクトル次元数
  int datanum = 0;
  int split = 8; //分割数
  double errorsum = 0; //平均二乗誤差
  double minerror = 99999999; //平均二乗誤差の平均の内最も小さいもの
  double bestSIGMA;
  double bests;
  double bestC;
  double bestc;
  FILE *gp, *fp;
  double average = 0; //1日目入札額の平均
  double day1max = 0; //1日目入札額の最大
  int adv = 0; //広告掲示回数
  
  srand((unsigned)time(NULL));



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
    printf("第三引数に引数の数を入力\n");
    exit(1);
  }



  double vector[MATRIX_DIM][VECTOR_DIM]; //特徴ベクトル1日目
  double day2[MATRIX_DIM][VECTOR_DIM]; //特徴ベクトル2日目
  int timestamp[MATRIX_DIM][TIME_DIM]; //時系列データ
  int prednum = 0; //特徴ベクトル2日目の数


  datanum = getmatrix(argv[2], vector, timestamp, dim, day2, &prednum); //時系列データを取得
  n = datanum;


for(int i = 0; i < n; i++){ //引数の数に合わせて特徴ベクトルを作成
  for(int j = 1; dim-j >= 0; j++){
    if(i-j >= 0){
      vector[i][dim-j] = vector[i-j][dim];
    }else{
      vector[i][dim-j] = vector[0][dim];
    }
  }
}


for(int i = 0; i < n; i++){ //引数の数に合わせて特徴ベクトルを作成
  for(int j = 1; dim-j >= 0; j++){
    if(i-j >= 0){
      day2[i][dim-j] = day2[i-j][dim];
    }else{
      day2[i][dim-j] = day2[0][dim];
    }
  }
}


fp = fopen("aaa.dat", "w"); //デバッグ用
for(int i = 0; i < prednum; i++){
  fprintf(fp, "%lf %lf %lf %lf\n", day2[i][0], day2[i][1], day2[i][2], day2[i][3]);
}
fclose(fp);


for(int i = 0; i < n; i++){ //1日目入札データの平均、最大を調べる
  average += vector[i][dim];
  if(vector[i][dim] > day1max){
    day1max = vector[i][dim];
  }
}
average = average / n;
double budget = average * (prednum / 3); //予算
printf("1日目落札額平均 : %lf\t最大落札額 : %lf\n", average, day1max);
printf("予算 : %lf\n", budget);


  if((n % split) != 0){
    printf("データ数%dが分割数%dで割り切れません\n", n, split);
    exit(1);
  }


  int testnum = n / split; //評価データ数
  int tn = n - testnum; //訓練データ数


//C探索ループ開始開始//////////////////////////////////////////////////////////////
printf("実行中");
fflush(stdout);


  for(double c = -3; c < 5; c++){
    C = pow(2, c);


//SIGMA探索ループ開始/////////////////////////////////////////////////////////////


  for(double s = -3; s < 5; s++){
    if(atoi(argv[1]) == 0 || atoi(argv[1]) == 1){ //ガウスカーネル以外は再帰しない
      s = 100;
    }
    SIGMA = pow(2, s);


//交差検定ループ開始///////////////////////////////////////////////////////////////


for(int cvd = 0; cvd < split; cvd++){  //分割数を元に交差検定を行う
  //printf("テスト%d回目\n", cvd + 1);

  double train[MATRIX_DIM][VECTOR_DIM];  //訓練データ
  double test[MATRIX_DIM][VECTOR_DIM]; //評価データ

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






      for(int i = 0; i < 2*tn; i++){
        for(int j = 0; j < 2*tn; j++){
          if((i+1 <= tn) && (j+1 <= tn)){
            G[i][j] = calc(i, j, dim, train, atoi(argv[1]), SIGMA); //行列Gの要素を計算して代入
          }else if((i+1 <= tn) && (j+1 > tn)){
            G[i][j] = -1 * calc(i, j-tn, dim, train, atoi(argv[1]), SIGMA);
          }else if((i+1 > tn) && (j+1 > tn)){
            G[i][j] = calc(i-tn, j-tn, dim, train, atoi(argv[1]), SIGMA);
          }else if((i+1 > tn) && (j+1 <= tn)){
            G[i][j] = -1 * calc(i-tn, j, dim, train, atoi(argv[1]), SIGMA);
          }else{
            printf("行列計算エラー\n");
            exit(1);
          }

          if(i == j){
            if(atoi(argv[1]) == 1){
              G[i][j] += 1.0e-7;
            }else{
              G[i][j] += 1.0e-9;
            }
          }
        }
      }


      for(int i = 0; i < 2*tn; i++){
        if(i+1 <= tn){
          g0[i] = EPSILON - train[i][dim];
        }else if(i+1 > tn){
          g0[i] = EPSILON + train[i-tn][dim];
        }else{
          printf("行列計算エラー\n");
          exit(1);
        }
      }


m = tn;
      for(int i = 0; i < 2*tn; i++){
        if(i+1 <= tn){
          CE[i][0] = 1;
        }else if(i+1 > tn){
          CE[i][0] = -1;
        }else{
          printf("行列計算エラー\n");
          exit(1);
        }
      }


      ce0[0] = 0.0;


p = 1;



      for(int i = 0; i < 2*tn; i++){
        for(int j = 0; j < 4*tn; j++){
          if(j+1 <= 2*tn){
            if(i == j){
              CI[i][j] = 1;
            }else{
              CI[i][j] = 0;
            }
          }else if(j+1 > 2*tn){
            if(i == (j-2*tn)){
              CI[i][j] = -1;
            }else{
              CI[i][j] = 0;
            }
          }else{
            printf("行列計算エラー\n");
            exit(1);
          }
        }
      }



      for(int i = 0; i < 4*tn; i++){
        if(i+1 <= 2*tn){
          ci0[i] = 0.0;
        }else if(i+1 > 2*tn){
          ci0[i] = C;
        }else{
          printf("行列計算エラー\n");
          exit(1);
        }
      }


int solven = 2*tn;
int solvem = 4*tn;


solve_quadprog(G, g0, solven, CE, ce0, p, CI, ci0, solvem, alpha); //２次計画問題を解く


for(int i = 0; i < 2*tn; i++){ //0以下のalphaを0とする
  if(alpha[i] < 0.1){
    alpha[i] = 0;
  }
}




double w[MATRIX_DIM] = {0};
for(int a = 0; a < dim; a++){ //重み計算
  for(int i = 0; i < tn; i++){
    w[a] += (alpha[i] - alpha[i+tn]) * train[i][a];
  }
}



double theta = 0.0; //閾値
double thetasum = 0;
double thetanum = 0;


for(int i = 0; i < tn; i++){ // thetaを計算
  theta = 0;

  if(alpha[i] > 0 && alpha[i] < C){
    for(int k = 0; k < tn; k++){
      theta += (alpha[k] - alpha[k+tn]) * calc(i, k, dim, train, atoi(argv[1]), SIGMA);
    }
    theta += -train[i][dim] + EPSILON;
    thetasum += theta;
    thetanum += 1;
  }


  if(alpha[i+tn] > 0 && alpha[i+tn] < C){
    for(int k = 0; k < tn; k++){
      theta += (alpha[k] - alpha[k+tn]) * calc(i, k, dim, train, atoi(argv[1]), SIGMA);
    }
    theta += -train[i][dim] - EPSILON;
    thetasum += theta;
    thetanum += 1;
  }
}

theta = thetasum / thetanum;


/*
printf("閾値θ: %lf\n", theta);

printf("重みw: ");
for(int a = 0; a < dim; a++){
  printf("%lf ", w[a]);
}
printf("\n");
*/


double error = 0; //平均二乗誤差


if(atoi(argv[1]) == 0){ //カーネルトリックなしの場合
  for(int i = 0; i < testnum; i++){
    double predval = 0;


    for(int a = 0; a < dim; a++){
      predval += w[a] * test[i][a] - theta;
    }

    error += pow(predval - test[i][dim], 2);
  }

  //printf("平均二乗誤差 : %lf\n", error / testnum);
}



if(atoi(argv[1]) == 1){ //多項式カーネルの場合
  for(int i = 0; i < testnum; i++){
    double predval = 0;


    for(int k = 0; k < tn; k++){
      double inner = 0;
      double kernel = 0;

      for(int a = 0; a < dim; a++){
        inner += vector[k][a] * vector[i][a];
      }
      kernel = pow(1.0 + inner, dim + 1);
      predval += (alpha[k] - alpha[k+tn]) * kernel;
    }

    predval = predval - theta;
    error += pow(predval - test[i][dim], 2);
  }

  //printf("平均二乗誤差 : %lf\n", error / testnum);
}



if(atoi(argv[1]) == 2){ //ガウスカーネルの場合
  for(int i = 0; i < testnum; i++){
    double predval = 0;


    for(int k = 0; k < tn; k++){
      double norm = 0;
      double kernel = 0;

      for(int a = 0; a < dim; a++){
        norm += pow(train[k][a] - test[i][a], 2);
      }
      kernel = exp(-norm / (2 * pow(SIGMA, 2)));
      predval += (alpha[k] - alpha[k+tn]) * kernel;
    }

    predval = predval - theta;
    error += pow(predval - test[i][dim], 2);
  }

  //printf("平均二乗誤差 : %lf\n", error / testnum);
}


errorsum += error / testnum; //平均二乗誤差の合計
} //交差検定を繰り返す


/*
if(atoi(argv[1]) == 2){
  printf("\nσ : %lf(s = %lf)\t\tC : %lf(c = %lf)\t\t平均二乗誤差の平均 : %lf\n", SIGMA, s, C, c, errorsum / split);
}else{
  printf("\nC : %lf(c = %lf)\t\t平均二乗誤差の平均 : %lf\n", C, c, errorsum / split);
}
*/


if((errorsum / split) < minerror){  //minerrorが最小になる場所を探す
  minerror = errorsum / split;
  bestSIGMA = SIGMA;
  bests = s;
  bestC = C;
  bestc = c;
}


errorsum = 0;
printf(".");
fflush(stdout);
} //SIGMA探索を繰り返す
} //C探索を繰り返す

/*
printf("\nパラメータ探索終了\n");
if(atoi(argv[1]) == 2){
  printf("bestσ : %lf(s = %lf)\t\tbestC : %lf(c = %lf)\t\t平均二乗誤差の平均の最小 : %lf\n", bestSIGMA, bests, bestC, bestc, minerror);
}else{
  printf("bestC : %lf(c = %lf)\t\t平均二乗誤差の平均の最小 : %lf\n", bestC, bestc, minerror);
}
*/
printf("\n");


//ここから回帰を開始///////////////////////////////////////////////////////////////



//printf("\n回帰プログラム開始\n");
C = bestc;
SIGMA = bestSIGMA;


for(int i = 0; i < 2*n; i++){
  for(int j = 0; j < 2*n; j++){
    if((i+1 <= n) && (j+1 <= n)){
      G[i][j] = calc(i, j, dim, vector, atoi(argv[1]), SIGMA); //行列Gの要素を計算して代入
    }else if((i+1 <= n) && (j+1 > n)){
      G[i][j] = -1 * calc(i, j-n, dim, vector, atoi(argv[1]), SIGMA);
    }else if((i+1 > n) && (j+1 > n)){
      G[i][j] = calc(i-n, j-n, dim, vector, atoi(argv[1]), SIGMA);
    }else if((i+1 > n) && (j+1 <= n)){
      G[i][j] = -1 * calc(i-n, j, dim, vector, atoi(argv[1]), SIGMA);
    }else{
      printf("行列計算エラー\n");
      exit(1);
    }


    if(i == j){
      if(atoi(argv[1]) == 1){
        G[i][j] += 1.0e-7;
      }else{
        G[i][j] += 1.0e-9;
      }
    }
  }
}


for(int i = 0; i < 2*n; i++){
  if(i+1 <= n){
    g0[i] = EPSILON - vector[i][dim];
  }else if(i+1 > n){
    g0[i] = EPSILON + vector[i-n][dim];
  }else{
    printf("行列計算エラー\n");
    exit(1);
  }
}


m = n;
for(int i = 0; i < 2*n; i++){
  if(i+1 <= n){
    CE[i][0] = 1;
  }else if(i+1 > n){
    CE[i][0] = -1;
  }else{
    printf("行列計算エラー\n");
    exit(1);
  }
}


ce0[0] = 0.0;


p = 1;



for(int i = 0; i < 2*n; i++){
  for(int j = 0; j < 4*n; j++){
    if(j+1 <= 2*n){
      if(i == j){
        CI[i][j] = 1;
      }else{
        CI[i][j] = 0;
      }
    }else if(j+1 > 2*n){
      if(i == (j-2*n)){
        CI[i][j] = -1;
      }else{
        CI[i][j] = 0;
      }
    }else{
      printf("行列計算エラー\n");
      exit(1);
    }
  }
}



for(int i = 0; i < 4*n; i++){
  if(i+1 <= 2*n){
    ci0[i] = 0.0;
  }else if(i+1 > 2*n){
    ci0[i] = C;
  }else{
    printf("行列計算エラー\n");
    exit(1);
  }
}


int solven = 2*n;
int solvem = 4*n;


solve_quadprog(G, g0, solven, CE, ce0, p, CI, ci0, solvem, alpha); //２次計画問題を解く


for(int i = 0; i < 2*n; i++){ //0以下のalphaを0とする
if(alpha[i] < 0.1){
alpha[i] = 0;
}
}




double w[MATRIX_DIM] = {0};
for(int a = 0; a < dim; a++){ //重み計算
for(int i = 0; i < n; i++){
w[a] += (alpha[i] - alpha[i+n]) * vector[i][a];
}
}



double theta = 0.0; //閾値
double thetasum = 0;
double thetanum = 0;


for(int i = 0; i < n; i++){ // thetaを計算
theta = 0;

if(alpha[i] > 0 && alpha[i] < C){
for(int k = 0; k < n; k++){
theta += (alpha[k] - alpha[k+n]) * calc(i, k, dim, vector, atoi(argv[1]), SIGMA);
}
theta += -vector[i][dim] + EPSILON;
thetasum += theta;
thetanum += 1;
}


if(alpha[i+n] > 0 && alpha[i+n] < C){
for(int k = 0; k < n; k++){
theta += (alpha[k] - alpha[k+n]) * calc(i, k, dim, vector, atoi(argv[1]), SIGMA);
}
theta += -vector[i][dim] - EPSILON;
thetasum += theta;
thetanum += 1;
}
}

theta = thetasum / thetanum;

/*
printf("閾値θ: %lf\n", theta);

printf("重みw: ");
for(int a = 0; a < dim; a++){
printf("%lf ", w[a]);
}
printf("\n");
*/


fp = fopen("trueval.dat", "w"); //真値ファイル

for(int i = 0; i < prednum; i++){
  fprintf(fp, "%lf\n", day2[i][dim]);
}

gp = popen("gnuplot -persist", "w");
fprintf(gp, "plot \"trueval.dat\"\n");

fclose(fp);


fp = fopen("predval.dat", "w"); //予測値ファイル
printf("\nオークション開始\n");


if(atoi(argv[1]) == 0){ //カーネルトリックなしの場合
  for(int i = 0; i < prednum; i++){
    double predval = 0;


    for(int a = 0; a < dim; a++){
      predval += w[a] * day2[i][a] - theta;
    }

  fprintf(fp, "%lf\n", predval);
  }
}


if(atoi(argv[1]) == 1){ //多項式カーネルの場合
  for(int i = 0; i < prednum; i++){
    double predval = 0;


    for(int k = 0; k < n; k++){
      double inner = 0;
      double kernel = 0;

      for(int a = 0; a < dim; a++){
        inner += vector[k][a] * day2[i][a];
      }
      kernel = pow(1.0 + inner, dim + 1);
      predval += (alpha[k] - alpha[k+n]) * kernel;
    }

    predval = predval - theta;
    fprintf(fp, "%lf\n", predval);
  }
}


if(atoi(argv[1]) == 2){ //ガウスカーネルの場合
  for(int i = 0; i < prednum; i++){
    double predval = 0;


    for(int k = 0; k < n; k++){
      double norm = 0;
      double kernel = 0;

      for(int a = 0; a < dim; a++){
        norm += pow(vector[k][a] - day2[i][a], 2);
      }
      kernel = exp(-norm / (2 * pow(SIGMA, 2)));
      predval += (alpha[k] - alpha[k+n]) * kernel;
    }

    predval = predval - theta;
    fprintf(fp, "%lf\n", predval);


    printf("予測値 : %lf  ", predval);

    if(budget > (prednum - i)*average){ //残予算が、残オークション回数*1日目落札平均額を上回る場合
      average = day1max; //1日目の最大落札額を入札の目安に変更する
    }

    if(predval <= average && (predval+0.3) > day2[i][dim]){ //入札&落札
      if(budget < (predval+0.3)){ //残予算チェック
        printf("予算が足りません");
      }else if((rand()%6+1) <= 3){
        budget = budget - (predval+0.3);
        adv++;
        printf("落札額 : %lf  残予算 : %lf     ", (predval+0.3), budget);
      }
    }
    printf("\n");
  }
}

printf("広告掲示回数 : %d\n", adv);

fprintf(gp, "replot \"predval.dat\"\n");
fclose(gp);


return 1;
}





int getmatrix(char *fname, double vector[MATRIX_DIM][VECTOR_DIM], int timestamp[MATRIX_DIM][TIME_DIM], int dim, double day2[MATRIX_DIM][VECTOR_DIM], int *prednum) //時系列データを取得
{
  FILE *fp;
  char buf[4][256];
  int i = 0;
  int vec = 0;
  double maxdata = 0;
  double getdata[MATRIX_DIM][VECTOR_DIM]; //一時的にデータを格納



  fp = fopen(fname, "r");
  if(fp == NULL){
    printf("%sファイルが開けません\n", fname);
    exit(1);
  }


  fscanf(fp, "%[^,],%[^,],%[^,],%s", buf[0], buf[1], buf[2], buf[3]); //最初の行を読み飛ばす


while(1){ //1日目のデータを得る

  while(1){ //同じ時間の入札データを得る
    fscanf(fp, "%d/%d/%d %d:%d:%d,%lf,%lf,%lf", &timestamp[i][0], &timestamp[i][1], &timestamp[i][2], &timestamp[i][3], &timestamp[i][4], &timestamp[i][5], &getdata[i][0], &getdata[i][1], &getdata[i][2]);

    if(i != 0){
      if(timestamp[i-1][4] != timestamp[i][4] || timestamp[i-1][3] != timestamp[i][3]){
        break;
      }
    }

    i++;
  }

  for(int j = 0; j < i; j++){
    if(maxdata < getdata[j][2]){
      maxdata = getdata[j][2];
    }
  }

  vector[vec][dim] = maxdata;
  vec++;
  maxdata = 0;


  if(timestamp[i-1][1] != timestamp[i][1]){ //もし読み込んだ次の時間が翌日であれば、そこで1日目終了
    for(int a = 0; a <= 5; a++){ //読み込んでしまった次の時間のデータを次のループへ持ち越す
      timestamp[0][a] = timestamp[i][a];
    }

    for(int a = 0; a <= 2; a++){ //読み込んでしまった次の時間のデータを次のループへ持ち越す
      getdata[0][a] = getdata[i][a];
    }

    break;
  }


  for(int a = 0; a <= 5; a++){ //読み込んでしまった次の時間のデータを次のループへ持ち越す
    timestamp[0][a] = timestamp[i][a];
  }

  for(int a = 0; a <= 2; a++){ //読み込んでしまった次の時間のデータを次のループへ持ち越す
    getdata[0][a] = getdata[i][a];
  }

  i = 1;
}

i = 1;

while(1){ //2日目のデータを得る

  while(1){ //同じ時間の入札データを得る
    fscanf(fp, "%d/%d/%d %d:%d:%d,%lf,%lf,%lf", &timestamp[i][0], &timestamp[i][1], &timestamp[i][2], &timestamp[i][3], &timestamp[i][4], &timestamp[i][5], &getdata[i][0], &getdata[i][1], &getdata[i][2]);

    if(i != 0){
      if(timestamp[i-1][4] != timestamp[i][4] || timestamp[i-1][3] != timestamp[i][3]){
        break;
      }
    }

    i++;
  }

  for(int j = 0; j < i; j++){
    if(maxdata < getdata[j][2]){
      maxdata = getdata[j][2];
    }
  }

  day2[*prednum][dim] = maxdata;
  *prednum = *prednum + 1;
  maxdata = 0;


  if(timestamp[i-1][1] != timestamp[i][1]){ //もし読み込んだ次の時間が翌日であれば、そこで2日目終了
    break;
  }


  for(int a = 0; a <= 5; a++){ //読み込んでしまった次の時間のデータを次のループへ持ち越す
    timestamp[0][a] = timestamp[i][a];
  }

  for(int a = 0; a <= 2; a++){ //読み込んでしまった次の時間のデータを次のループへ持ち越す
    getdata[0][a] = getdata[i][a];
  }

  i = 1;
}


return vec;
}





double calc(int i, int j, int dim, double vector[MATRIX_DIM][VECTOR_DIM], int num, double SIGMA) //カーネルを計算
{
  double inner = 0;

  for(int a = 0; a < dim; a++){
    inner += vector[i][a] * vector[j][a]; //内積を計算
  }

if(num == 0){
  return inner; // (x_i, x_j)
}else if(num == 1){
  return pow(1.0 + inner, dim + 1); //多項式カーネル
}else if(num == 2){ //ガウスカーネル
  double norm = 0.0;

  for(int a = 0; a < dim; a++){
    norm += pow(vector[i][a] - vector[j][a], 2);
  }
  return exp(-norm / (2 * pow(SIGMA, 2)));
}else{
  printf("カーネル計算エラー\n");
  return 1.0;
}
}
