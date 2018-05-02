//  10.04.2018 Weijie YE

# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>
# include <string.h>
# include <math.h>
# include <time.h>
# include <gmp.h>
#include "fonctions.h"


//on test les deux fonctoins avec les meme coefficients,a,k ;
int main(){



    unsigned long int deg;
    clock_t t1, t2; // les deux variable pour stocker le temps .
    float temps,temps2,temps3;
    printf("saisir le degre max au quel on veut atteindre /\nn=");

    scanf("%lu", &deg);


    mpz_t *tab = malloc((deg+1)*sizeof(mpz_t));

    //x=2puissance (le troisieme argument), x represente le nombre maximum qu on peux avoir;
    set_coefficients("coefficients_p",deg,48); //Ecriture de coefficients de type mpz_t aléatoires dans le fichier avec degre deg


    //Lecture des coefficients initiqlise par set_coefficients et de l ecrire dans le tableau tab;
    parse_file("coefficients_p",tab,deg);

    long long a, k;
    printf("Calcul de q(a/(2^k)), veuillez rentrer les valeurs de:\na:\n");
    scanf("%llu", &a);
    printf("\nk:\n");
    scanf("%llu", &k);
    mpz_t res_num1, res_denom1, res_num2, res_denom2,  res_num3, res_denom3;
     mpz_init(res_num1);//Initialistaion de numerateur pour la premier methode
     mpz_init(res_denom1);//Initialistaion de denominateur
     mpz_init(res_num2);//Pour la deuxieme methode
     mpz_init(res_denom2);
     mpz_init(res_num3);//Pour la methode d'horner
     mpz_init(res_denom3);


    effacer_fichier("E1_f1_num.txt");
    effacer_fichier("E1_f1_denom.txt");
    effacer_fichier("E1_temps.txt");
    effacer_fichier("E1_f2_num.txt");
    effacer_fichier("E1_f2_denom.txt");
    effacer_fichier("E1_Horner_num.txt");
    effacer_fichier("E1_Horner_denom.txt");



    unsigned long int i;

    for(i=0; i<=deg; i++){

        t1=clock();

        eval_poly_1(tab, a, k, i, &res_num1, &res_denom1);
        //Pour affichier le resultat a chaque incrementation de n
        //gmp_printf("\n\nResultat méthode 1: num = %Zd\t \n den = %Zd\n", res_num1, res_denom1);

        t2=clock();
        //printf("%lu\n",i);
        temps=(float)(t2-t1)/CLOCKS_PER_SEC;

        t1=clock();
        eval_poly_2(tab, a, k, i, &res_num2, &res_denom2);
        //gmp_printf("\n\nResultat méthode 2: num = %Zd\t \n den = %Zd\n", res_num2, res_denom2);
        t2=clock();

        temps2=(float)(t2-t1)/CLOCKS_PER_SEC;
        t1=clock();
        evaluateHorner(tab, a, k, i, &res_num3, &res_denom3);
        //gmp_printf("\n\nResultat méthode 3: num = %Zd\t \n den = %Zd\n", res_num3, res_denom3);
        t2=clock();

        temps3=(float)(t2-t1)/CLOCKS_PER_SEC;




        ecrit_result("E1_f1_num.txt", res_num1);//FICHIER qui contient les numerateurs pour M1
        ecrit_result("E1_f1_denom.txt", res_denom1);//FICHIER qui contient les dénominateurs M1

        ecrit_result("E1_f2_num.txt", res_num2);
        ecrit_result("E1_f2_denom.txt", res_denom2);

        ecrit_result("E1_Horner_num.txt", res_num3);
        ecrit_result("E1_Horner_denom.txt", res_denom3);

        //on multiplie le temps par 10000 pour insertion dans excel
        enregistrement_3("E1_temps.txt",(unsigned int)i,temps*10000,temps2*10000,temps3 *10000);




    }
    gmp_printf("\n\nResultat méthode 1: num = %Zd\t den = %Zd\n", res_num1, res_denom1);// resultat final
    gmp_printf("\n\nResultat méthode 2: num = %Zd\t den = %Zd\n", res_num2, res_denom2);
    gmp_printf("\n\nResultat méthode 3: num = %Zd\t den = %Zd\n", res_num3, res_denom3);

    printf(" \nfin\n");


    free(tab);
    return 0;
}
