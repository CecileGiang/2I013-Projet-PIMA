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
    clock_t t1, t2;
    float temps,temps2;
    unsigned int max;
    printf("max des coeff \nmax=");
    scanf("%u", &max);
    
    printf("saisir le degre max \nn=");
    scanf("%lu", &deg);
    

    mpz_t *tab = malloc((deg+1)*sizeof(mpz_t));
       long long a, k;
    printf("Calcul de q(a/(2^k)), veuillez rentrer les valeurs de:\na:\n");
    scanf("%llu", &a);
    printf("\nk:\n");
    scanf("%llu", &k);
    mpz_t res; mpz_init(res);
    mpz_t res1; mpz_init(res1);
    unsigned long int i;
    /* Evaluation de polynôme : les valeurs de a et k seront demandées.*/
    /*
     for(i=0; i<=deg; i++){
     
     gmp_printf("%Zd\n",tab[i]);
     
     }*/
    
    
    for(i=1; i<=max; i++){
        
        set_coefficients("coefficients_p",deg,i); //Ecriture de coefficients de type mpz_t aléatoires dans le fichier
        parse_file("coefficients_p",tab,deg); //Ecriture des coefficients dans le tableau tab

        t1=clock();
        //    eval_poly_1(tab, a, k, i, &res);
        gmp_printf("Resultat méthode 1: %Lf\n", eval_poly_1(tab, a, k, deg, &res));
        t2=clock();
        //printf("%lu\n",i);

        temps=(float)(t2-t1)/CLOCKS_PER_SEC;
        
        // enregistrement("temps_methode1.txt",(unsigned int)i,temps*10000);
        
        t1=clock();
        // eval_poly_2(tab, a, k, i, &res1);
        gmp_printf("Resultat méthode 2: %Lf\n", eval_poly_2(tab, a, k, deg, &res1));
        t2=clock();
        //printf("%lu\n",i);

        temps2=(float)(t2-t1)/CLOCKS_PER_SEC;
        
        ecrit_result("eval2_resultat_f1.txt", res);
        ecrit_result("eval2_resultat_f2.txt", res1);
        
        enregistrement("eval2_temps.txt",(unsigned int)i,temps*10000,temps2*10000);
        
        
        
    }
    printf("\nfin");
    
    
    free(tab);
    return 0;
}

