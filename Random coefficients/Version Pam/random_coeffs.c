#include <gmp.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>


/* Cette fonction ecrit des coefficients entiers de type mpz sur un fichier donne en argument. 
Ceci se fait en implantantant un etat random a l aide de l algorithme de Mersenne Twister. 
Ce noyau donne plus tard des valeurs aleatoires (entre 0 et 2^(n-1), n etant le troisieme argument de la fonction)
a chaque coefficient.
Comme nous souhaitons aussi avoir des nombres negatifs, nous les generons aleatoirement, c est-a-dire, nous creons une variable signe,
qui determinera si un nombre est positif ou negatif avec une equiprobabilite */

void random_coeff(char *nom_fichier, unsigned long int deg) {
	
	FILE *f = fopen(nom_fichier, "w");

	if (f == NULL){
		printf("Erreur lors de l ouverture du fichier\n");
		exit(1);
	}

	unsigned long int i;
	int sign;
	mpz_t coeff;
	
	srand(time(NULL));
	
	gmp_randstate_t state; 
	gmp_randinit_mt(state); //initialise state pour l algorithme de Mersenne Twister 
	

	mpz_init(coeff);       //initialisation du coefficient a zero


	for (i = 0; i <= deg; i++){
		mpz_rrandomb(coeff, state, 120); //coeff prend une valeur random entre 0 et 2**(n-1) grace a l algorithme pris par state
		if (rand()%2 == 0){ //nous creons artificiellement des entiers negatifs aussi, avec une chance sur deux
			sign = 1;
		}
		else {
			sign = -1;
		}
		mpz_mul_si(coeff, coeff, sign); 
		mpz_out_str(f, 10, coeff); //nous imprimons le nombre dans le fichier
		fputs("\n", f); //retour a la ligne
	}
	gmp_randclear(state);
	mpz_clear(coeff);
	fclose(f);


//commentaires: nous devrions avoir des nombres aleatoires entre 0 et 2**(n-1), pourtant nous trouvons que des nombres avec le meme nombre de chiffres a chaque fois
}



int main(){
	unsigned long int deg;
	printf("Veuillez saisir le degre du polynome: ");
	scanf("%lu", &deg);

	random_coeff("coeffs.txt", deg);
	return 0;
}
