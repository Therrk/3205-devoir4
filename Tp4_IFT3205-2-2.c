//------------------------------------------------------
// Prog    : Tp4_IFT3205                          
// Auteur  : Élie Leblanc, Justin Veilleux
// Date    :                                  
// version :                                             
// langage : C                                          
// labo    : DIRO                                       
//------------------------------------------------------

//------------------------------------------------
// FICHIERS INCLUS -------------------------------
//------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "FonctionDemo4.h"
#include "fonctions.c"

//------------------------------------------------
// DEFINITIONS -----------------------------------   
//------------------------------------------------
#define NAME_VISUALISER "display "
#define NAME_IMG_IN  "cameraman"
#define NAME_IMG_OUT1 "cameraman_restored"
#define NAME_IMG_OUT2 "cameraman_degraded"

//------------------------------------------------
// PROTOTYPE DE FONCTIONS  -----------------------   
//------------------------------------------------

//------------------------------------------------
// VARIABLES GLOBALES  ---------------------------   
//------------------------------------------------

//>Taille Image
int length=256;
int width=256;
int size_image=256;

//>Parametre de degradation
int size=9;
float variance_noise=0.5;
int psf=1;

//>Nb d'Iterations
int nbiter=20;

//>ImprovmentSNR
float isnr;

//------------------------------------------------
//------------------------------------------------
// FONCTIONS  ------------------------------------   
//------------------------------------------------
//------------------------------------------------


//---------------------------------------------------------
//---------------------------------------------------------
// PROGRAMME PRINCIPAL   ----------------------------------                     
//---------------------------------------------------------
//---------------------------------------------------------
int main(int argc,char **argv)
 {
  int i,j,k;
  char BufSystVisuImg[100];
  float total = 0;
  float denum, num, isnr;

  float** mat_img=LoadImagePgm(NAME_IMG_IN,&length,&width);
  
  //Allocation memoire matrice
  float** mat=fmatrix_allocate_2d(length,width);
  float** mat_rest=fmatrix_allocate_2d(length,width);
  float** mat_rest_prec=fmatrix_allocate_2d(length,width); 
  float** mat_rest_best=fmatrix_allocate_2d(length,width); 
  float** mat_psf=fmatrix_allocate_2d(length,width);

  float** mat_I=fmatrix_allocate_2d(length,width);
  float** mat_mod=fmatrix_allocate_2d(length,width);
  float** mat_h_R=fmatrix_allocate_2d(length,width);
  float** mat_h_I=fmatrix_allocate_2d(length,width);
  float** mat_rest_I=fmatrix_allocate_2d(length,width);
  float** mat_h_conj_R=fmatrix_allocate_2d(length,width);
  float** mat_h_conj_I=fmatrix_allocate_2d(length,width);
  float** mat_h_mod=fmatrix_allocate_2d(length,width);
  float** mat_tmp10=fmatrix_allocate_2d(length,width);
  float** mat_tmp11=fmatrix_allocate_2d(length,width);
  float** mat_tmp12=fmatrix_allocate_2d(length,width);

    // zeros
  fmatrix_zero(length, width, mat_h_R);
  fmatrix_zero(length, width, mat_h_I);
  fmatrix_zero(length, width, mat_I);
  fmatrix_zero(length, width, mat_mod);
  fmatrix_zero(length, width, mat_rest);
  fmatrix_zero(length, width, mat_rest_I);
  fmatrix_zero(length, width, mat_h_conj_R);
  fmatrix_zero(length, width, mat_h_conj_I);
  fmatrix_zero(length, width, mat_h_mod);
 
  //=========================================================
  //== PROG =================================================
  //=========================================================

  //>Lecture Image

  //--------------------------------------------------------
  //>Degradation
  //--------------------------------------------------------
  // Cette fonction ajoute un flou cr�� par une psf uniforme 
  // (fonction porte) de taille sizexsize puis ajoute sur
  // cette image rendue floue, un bruit Gaussien de variance
  // variance_noise
  //
  // Entr�e : mat_img :: image non d�grad�e
  // Sortie : mat     :: image d�grad�e
  //--------------------------------------------------------
  degradation(mat_img,mat,length,width,psf,size,variance_noise);

  //============
  //WIENER
  //============
  // Création h(x,y)
  for (i = 0; i < length; i++) {
	  for (j = 0; j < width; j++) {
	    if (i<=((length/2)+4)&&i>=((length/2)-4)&&j<=((length/2)+4)&&j>=((length/2)-4)) {
       mat_h_R[i][j]=1;
       total+=1;
      }
    }
  }
  
  // centrage correct
  CenterImg(mat_h_R, length, width);
  
  // normalisation
  for (i = 0; i < length; i++) {
	  for (j = 0; j < width; j++) {
      mat_h_R[i][j]*=((length*width)/total);
    }
  }

  FFTDD(mat_h_R, mat_h_I, length, width);
  FFTDD(mat, mat_I, length, width);

  // H conjugué
  for (i = 0; i < length; i++) {
  	for (j = 0; j < width; j++) {
	    mat_h_conj_R[i][j]=mat_h_R[i][j];
	    mat_h_conj_I[i][j]=-mat_h_I[i][j];
    }
  }
  
  // modules au carré
  fmatrix_module(length, width, mat_h_R, mat_h_I, mat_h_mod);
  fmatrix_module(length, width, mat, mat_I, mat_mod);
  for (i = 0; i < length; i++) {
    	for (j = 0; j < width; j++) {
  	    mat_h_mod[i][j]= CARRE(mat_h_mod[i][j]);
  	    mat_mod[i][j]= CARRE(mat_mod[i][j]);
      }
  }
  // restauration
  for (i = 0; i < length; i++) {
  	for (j = 0; j < width; j++) {
	    denum = mat_h_mod[i][j]+(variance_noise/(length*width))/mat_mod[i][j];
      mat_h_conj_R[i][j] = mat_h_conj_R[i][j]/denum;
      mat_h_conj_I[i][j] = mat_h_conj_I[i][j]/denum;
    }
  }
  MultMatrix(mat_rest, mat_rest_I, mat, mat_I, mat_h_conj_R, mat_h_conj_I, length, width);
  
  IFFTDD(mat_rest, mat_rest_I, length, width);
  Recal(mat_rest, length, width);
  IFFTDD(mat, mat_I, length, width);

  // wow, ça marche!
  // Calcul isnr
  num = 0;
  denum = 0;
  for (i = 0; i < length; i++) {
	  for (j = 0; j < width; j++) {
	    num+=CARRE(mat_img[i][j]-mat[i][j]);
	    denum+=CARRE(mat_img[i][j]-mat_rest[i][j]);
    }
  }
  isnr = 10*log10(num/denum);
  printf("isnr: %f\n",isnr);
  //---------------------------------------------
  // SAUVEGARDE et VISU
  // -------------------
  // Le resultat de la restoration > mat_rest
  // L'image d�grad�e              > mat
  // L'image non d�grad�e          > mat_img
  //----------------------------------------------
  SaveImagePgm(NAME_IMG_OUT1,mat_rest,length,width);
  SaveImagePgm(NAME_IMG_OUT2,mat,length,width);
  
  strcpy(BufSystVisuImg,NAME_VISUALISER);
  strcat(BufSystVisuImg,NAME_IMG_OUT2);
  strcat(BufSystVisuImg,".pgm&");
  printf("\n > %s",BufSystVisuImg);
  system(BufSystVisuImg);
  strcpy(BufSystVisuImg,NAME_VISUALISER);
  strcat(BufSystVisuImg,NAME_IMG_IN);
  strcat(BufSystVisuImg,".pgm&");
  printf("\n > %s",BufSystVisuImg);
  system(BufSystVisuImg);
  strcpy(BufSystVisuImg,NAME_VISUALISER);
  strcat(BufSystVisuImg,NAME_IMG_OUT1);
  strcat(BufSystVisuImg,".pgm&");
  printf("\n > %s",BufSystVisuImg);
  system(BufSystVisuImg);
  
  
  //Liberation memoire pour les matrices
  // if (mat)            free_fmatrix_2d(mat);
  // if (mat_img)        free_fmatrix_2d(mat_img);
  // if (mat_rest)       free_fmatrix_2d(mat_rest);
  // if (mat_rest_prec)  free_fmatrix_2d(mat_rest_prec);
  // if (mat_rest_best)  free_fmatrix_2d(mat_rest_best);
  // if (mat_psf)        free_fmatrix_2d(mat_psf);
  // if (mat_tmp0)  free_fmatrix_2d(mat_tmp0);
  // if (mat_tmp1)  free_fmatrix_2d(mat_tmp1);
  // if (mat_tmp2)  free_fmatrix_2d(mat_tmp2);
  // if (mat_tmp3)  free_fmatrix_2d(mat_tmp3);
  // if (mat_tmp4)  free_fmatrix_2d(mat_tmp4);
  // if (mat_tmp5)  free_fmatrix_2d(mat_tmp5);
  // if (mat_tmp6)  free_fmatrix_2d(mat_tmp6);
  // if (mat_tmp7)  free_fmatrix_2d(mat_tmp7);
  // if (mat_tmp8)  free_fmatrix_2d(mat_tmp8);
  // if (mat_tmp9)  free_fmatrix_2d(mat_tmp9);
  // if (mat_tmp10) free_fmatrix_2d(mat_tmp10);
  // if (mat_tmp11) free_fmatrix_2d(mat_tmp11);
  // if (mat_tmp12) free_fmatrix_2d(mat_tmp12);

  //retour sans probleme 
  printf("\n C'est fini ... \n\n");
  return 0; 	 
}


