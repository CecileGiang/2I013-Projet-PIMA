# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>
# include <string.h>
# include <math.h>
# include <time.h>
# include <gmp.h>
#include "fonctions.h"

int main(){
    
    /* Ecriture/lecture de polynôme */
    
    unsigned long int deg;
    printf("Génération d'un fichier de degré n/\nn=");
    scanf("%lu", &deg);
    
    mpz_t *tab = malloc((deg+1)*sizeof(mpz_t));
    
    set_coefficients("coefficients_p",deg); //Ecriture de coefficients de type mpz_t aléatoires dans le fichier
    parse_file("coefficients_p",tab,deg); //Ecriture des coefficients dans le tableau tab
    
    
    /* Evaluation de polynôme : les valeurs de a et k seront demandées.*/
    unsigned long int i;
    for(i=0; i<=deg; i++) gmp_printf("%Zd\n", tab[i]);
    
    long long a, k;
    printf("Calcul de q(a/(2^k)), veuillez rentrer les valeurs de:\na:\n");
    scanf("%llu", &a);
    printf("\nk:\n");
    scanf("%llu", &k);
    
    mpz_t res; mpz_init(res);
    mpz_t res1; mpz_init(res1);
    
    /* Méthode 1 */
    
    
    gmp_printf("Resultat méthode 1: %Lf\n", eval_poly_1(tab, a, k, deg, &res));
    
    /* Méthode 2 */
    
    gmp_printf("Resultat méthode 2: %Lf\n", eval_poly_2(tab, a, k, deg, &res1));
    
    
    free(tab);
    return 0;
}

