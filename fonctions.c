#include <math.h>
#include <stdio.h>

#include "FonctionDemo4.h"

void rotation(int length, int width, float angle, float ** mat_src, float ** mat_dst) {
    // correction rotation dans la mauvaise direction
    angle = -angle;
    
    // Le code
    for (int i = 0; i < length; i++) {
        for (int j = 0; j < width; j++) {

            int half_width = width/2;
            int half_length = length/2;

            // Les coordonnées centrés
            float x_centered = j - half_width;
            float y_centered = i - half_length;


            // La matrice de rotation est la suivante:


            // | cos(θ) -sin(θ) |
            // | sin(θ)  cos(θ) |
            
            // Puisqu'elle est orthonormale, son inverse est sa transposée

            // Les coordonnées où on devrait prendre les pixels
            float x_src = (x_centered*cos(angle) + y_centered*sin(angle)) + half_width;
            float y_src = (x_centered*-sin(angle) + y_centered*cos(angle)) + half_length;

            int up = (int)floor(y_src);
            int down = (int)ceil(y_src);
            int left = (int)floor(x_src);
            int right = (int)ceil(x_src);

            if (left>=0 && right < width && up >= 0 && down < length) {


                // Les quatre pixels réels qui encadrent notre pixel destination



                // à quel distance somme nous du point entier en haut à gauche?
                float x_off = x_src - (float)left;
                float y_off = y_src - (float)up;

                // interpolation linéaire pour le côté supérieur et inférieur du pixel

                /*           up          */
                /*        +------+       */
                /*   left |  le  | right */
                /*        | pixel|       */
                /*        +------+       */
                /*          down         */
                float down_avg = mat_src[down][left] + x_off * (mat_src[down][right] - mat_src[down][left]);
                float up_avg = mat_src[up][left] + x_off * (mat_src[up][right] - mat_src[up][left]);


                // interpolation entre les deux interpolations. le up, down est
                // inversé, car la composante verticale pointe vers le bas.
                float vert_avg = up_avg + y_off * (down_avg - up_avg);

                mat_dst[i][j]= (int)vert_avg;
            } else {
                mat_dst[i][j]= 0;
            }
        }
    }
}

int max(int a, int b) {
    if (a > b) {
        return a;
    }else{
        return b;
    }
}

int min(int a, int b) {
    if (a < b) {
        return a;
    }else{
        return b;
    }
}

float interpol_lin(float** img, int h, int w, float x_src, float y_src)
{
    int up = (int)floor(y_src);
    up = max(0, up);
    int down = (int)ceil(y_src);
    down = min(h - 1, down);
    int left = (int)floor(x_src);
    left = max(0, left);
    int right = (int)ceil(x_src);
    right = min(w - 1, right);

    float col_up =
        (x_src - left) * img[right][up] +
        (right - x_src) * img[left][up];
    float col_down = 
        (x_src - left) * img[right][down] +
        (right - x_src) * img[left][down];

    float col_vert = (y_src - up) * col_down +
        (down - y_src) * col_up;

    return col_vert;
}

void fmatrix_move(int length, int width, float ** src, float ** dst) {
    for (int i = 0;i < length;i++) {
        for (int j = 0;j < width;j++) {
            dst[i][j] = src[i][j];
        }
    }
}

void fmatrix_module(int length, int width, float ** r, float ** im, float ** dst) {
    for (int i = 0;i < length;i++) {
        for (int j = 0;j < width;j++) {
            float real = r[i][j];
            float ima = im[i][j];
            if (isnan(real) || isnan(ima)) {
                // printf("erreur dans module\n");
            }
            dst[i][j] = sqrtf(real * real + ima * ima);
        }
    }
}

void fmatrix_zero(int length, int width, float ** mat) {
    for (int i = 0;i < length;i++) {
        for (int j = 0;j < width;j++) {
            mat[i][j] = 0.0;
        }
    }
}

//Fonction pour centrer l'image, tirée verbatim de FonctionDemo2.c
void CenterImg(float** mat,int lgth,int wdth)
{
    int i,j;
    int ci,cj;
    float** mattmp;

    /*Initialisation*/
    ci=(int)(lgth/2);
    cj=(int)(wdth/2);

    /*Allocation memoire*/
    mattmp=fmatrix_allocate_2d(lgth,wdth);

    /*Recadrage*/
    for(i=0;i<ci;i++) for(j=0;j<cj;j++)
                          mattmp[ci+i][cj+j]=mat[i][j];

    for(i=ci;i<lgth;i++) for(j=cj;j<wdth;j++)
                             mattmp[i-ci][j-cj]=mat[i][j];

    for(i=0;i<ci;i++) for(j=cj;j<wdth;j++)
                          mattmp[ci+i][j-cj]=mat[i][j];

    for(i=ci;i<lgth;i++) for(j=0;j<cj;j++)
                             mattmp[i-ci][cj+j]=mat[i][j];

    /*Transfert*/
    for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
                            mat[i][j]=mattmp[i][j];

    /*desallocation memoire*/
    free_fmatrix_2d(mattmp);
}

float norm(float a, float b) {
    return sqrt(a * a + b * b);
}

void threshold(float**matR, float**matI, int length, int width, float norm_min) {
    for (int i = 0;i < length;i++) {
        for (int j = 0;j < width;j++) {
            if (norm(matR[i][j], matI[i][j]) < norm_min) {
                matR[i][j] = 0.0;
                matI[i][j] = 0.0;
            }
        }
    } 
}

// fonction pour multiplier 2 matrices complexes, tirée verbatim de FonctionDemo3.c
void MultMatrix(float** matRout,float** matIout,float** mat1Rin,float** mat1Iin,
float** mat2Rin,float** mat2Iin,int lgth,int wdth)
{
 int i,j;

 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
   { matRout[i][j]=mat1Rin[i][j]*mat2Rin[i][j]-mat1Iin[i][j]*mat2Iin[i][j];
     matIout[i][j]=mat1Rin[i][j]*mat2Iin[i][j]+mat2Rin[i][j]*mat1Iin[i][j]; }
}

void Recal(float** mat,int lgth,int wdth)
{
 int i,j;
 float max,min;

 /*Initialisation*/
 min=mat[0][0];

 /*Recherche du min*/
  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    if (mat[i][j]<min) min=mat[i][j];

 /*plus min*/
   for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    mat[i][j]-=min;

   max=mat[0][0];
 /*Recherche du max*/
  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
    if (mat[i][j]>max) max=mat[i][j];

 /*Recalibre la matrice*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
   mat[i][j]*=(GREY_LEVEL/max);      
}

// Recal l'image coupant les valeurs trop hauts et trop basses, en les remplaçant par 0 ou 255 
void Recal2(float** mat, int length, int width){
    int i,j;
    for (i = 0; i < length; i++) {
    	for (j = 0; j < width; j++) {
	        if (mat[i][j]<0) {
	            mat[i][j]=0;
            } else if (mat[i][j]>255) {
	            mat[i][j]=255;
            }
        }
    }
}
