#include <iostream>
#include <sstream>
#include <string>
#include "QuadProg++.hh"
#include <stdio.h>
#include <string.h>
#include <math.h>

#define SIGMA 10

int getmatrix(char *fname, int *dim, double vector[MATRIX_DIM][MATRIX_DIM]);
double calc(int i, int j, int dim, double vector[MATRIX_DIM][MATRIX_DIM], int num);
int maxalpha(double alpha[], int size);

int main(int argc, char *argv[])
{
  double G[MATRIX_DIM][MATRIX_DIM], g0[MATRIX_DIM],
		CE[MATRIX_DIM][MATRIX_DIM], ce0[MATRIX_DIM],
		CI[MATRIX_DIM][MATRIX_DIM], ci0[MATRIX_DIM],
		alpha[MATRIX_DIM];
	int n, m, p;
  int dim; //ベクトル次元数
  FILE *gp, *fp;

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

  double vector[MATRIX_DIM][MATRIX_DIM]; //特徴ベクトル

  n = getmatrix(argv[2], &dim, vector); //行列サイズを取得、同時に特徴ベクトルを配列に整理する


      for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
          G[i][j] = vector[i][dim] * vector[j][dim] * calc(i, j, dim, vector, atoi(argv[1])); //行列Gの要素を計算して代入
          if(i == j){
            if(atoi(argv[1]) == 1){
              G[i][j] += 1.0e-7;
            }else{
              G[i][j] += 1.0e-9;
            }
        }
      }
    }

      for(int i=0; i<n; i++){
        g0[i] = -1.0;
      }

m = n;
      for(int i=0; i<n; i++){
        CE[i][0] = vector[i][dim];
      }

      ce0[0] = 0.0;

p = 1;
      for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
          if(i == j){
            CI[i][j] = 1.0;
          }else{
            CI[i][j] = 0.0;
          }
        }
      }

      for(int i=0; i<n; i++){
        ci0[i] = 0.0;
      }

solve_quadprog(G, g0, n, CE, ce0, p, CI, ci0, m, alpha); //２次計画問題を解く

double w[MATRIX_DIM] = {0};
for(int a = 0; a < dim; a++){ //重み計算
  for(int i = 0; i < n; i++){
    w[a] += alpha[i] * vector[i][dim] * vector[i][a];
  }
}

double theta = 0.0; //閾値
int k =  maxalpha(alpha, n);//alphaが最大になる次数

if(atoi(argv[1]) == 0){ // カーネルトリックなしのtheta
  for(int a = 0; a < dim; a++){
    theta += w[a] * vector[k][a];
  }
  theta = theta - vector[k][dim];
}

if(atoi(argv[1]) == 1){ // 多項式カーネルのtheta
  double inner = 0;
  for(int i = 0; i < n; i++){
    for(int a = 0; a < dim; a++){
      inner += vector[i][a] * vector[k][a]; //内積を計算
    }
    theta += alpha[i] * vector[i][dim] * pow(1.0 + inner, dim + 1);
    inner = 0;
  }
}

if(atoi(argv[1]) == 2){ // ガウスカーネルのtheta
  double norm = 0.0;
  for(int i = 0; i < n; i++){
    for(int a = 0; a < dim; a++){
      norm += pow(vector[i][a] - vector[k][a], 2);
    }
    theta += alpha[i] * vector[i][dim] * exp(-norm / (2 * pow(SIGMA, 2)));
    norm = 0;
  }
}

printf("閾値θ: %lf\n", theta);

printf("重みw: ");
for(int a = 0; a < dim; a++){
  printf("%lf ", w[a]);
}
printf("\n");


fp = fopen("plusclass.dat", "w"); //クラス１のファイルを作成
for(int i = 0; i < n; i++){
  if(vector[i][dim] == 1){
    for(int a = 0; a < dim; a++){
      fprintf(fp, "%lf\t", vector[i][a]);
    }
    fprintf(fp, "\n");
  }
}
fclose(fp);

fp = fopen("minusclass.dat", "w"); //クラス-１のファイルを作成
for(int i = 0; i < n; i++){
  if(vector[i][dim] == -1){
    for(int a = 0; a < dim; a++){
      fprintf(fp, "%lf\t", vector[i][a]);
    }
    fprintf(fp, "\n");
  }
}
fclose(fp);


gp = popen("gnuplot -persist", "w");
fprintf(gp, "plot \"plusclass.dat\"\n"); //特徴ベクトルプロット
fprintf(gp, "replot \"minusclass.dat\"\n");


if(atoi(argv[1]) == 0){ //カーネルトリックなしの識別器
  double a, b; // y = ax + b

  a = -w[0] / w[1];
  b = theta / w[1];
  fprintf(gp, "f(x) = %lf * x + %lf\n", a, b);
  fprintf(gp, "replot f(x)\n");
}

fp = fopen("boundary.dat", "w");

if(atoi(argv[1]) == 1 || atoi(argv[1]) == 2){ //多項式カーネルorガウスカーネル
  double x[MATRIX_DIM] = {0}; //xベクトル
  double kernel = 0;
  double inner = 0;
  double norm = 0.0;
  double f = 0;
  double check;

  while(x[0] <= 50){
    while(x[1] <= 50){
      for(int i = 0; i < n; i++){

        if(atoi(argv[1]) == 1){ //多項式カーネルの場合
          for(int a = 0; a < dim; a++){
            inner += vector[i][a] * x[a];
          }
          kernel = pow(1.0 + inner, dim + 1);
        }

          if(atoi(argv[1]) == 2){ //ガウスカーネルの場合
            for(int a = 0; a < dim; a++){
              norm += pow(vector[i][a] - x[a], 2);
            }
            kernel = exp(-norm / (2 * pow(SIGMA, 2)));
          }

      f += alpha[i] * vector[i][dim] * kernel;
      norm = 0;
      inner = 0;
      kernel = 0;
      }
      check = f - theta;
      if(check > -0.05 && check < 0.05){
        fprintf(fp, "%lf\t%lf\n", x[0], x[1]);
      }
      x[1] += 0.1;
      f = 0;
    }
    x[0] += 0.1;
    x[1] = 0;
  }
  fprintf(gp, "replot \"boundary.dat\"\n");
  fclose(fp);
}

pclose(gp);

for(int i = 0; i < n; i++){
  printf("alpha[%d]: ", i+1);
  printf("%lf\n", alpha[i]);
}

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

double calc(int i, int j, int dim, double vector[MATRIX_DIM][MATRIX_DIM], int num) //カーネルを計算
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

int maxalpha(double alpha[], int size)
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
