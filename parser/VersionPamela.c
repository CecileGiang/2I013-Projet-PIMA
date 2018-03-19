# include <stdio.h>
# include <gmp.h>
# include <stdint.h>
# include <stdlib.h>

mpz_t* parser(char *fichier, long long int *degree){
    unsigned int i;
    FILE *f = fopen(fichier, "r");

    if (f == NULL) {
        printf("Erreur lors de l ouverture du fichier!\n");
        exit(1);
    }
    fscanf(f, "%lld", degree);
    mpz_t *tab = malloc(sizeof(mpz_t)*(*degree+1));

    for (i = 0; i < *degree; i++){
        fscanf(f, "%Zd", tab[i]);
    }
    
    fclose(f);
    return tab;
}
int main(){
    long long int degree;
    int i;

    mpz_t *poly;
    poly = parser("test1.txt", &degree);
    for (i = 0; i < degree; i++){
        printf("%Zd", poly[i]);
    } 
    return 0;
}
