# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>
# include <string.h>
# include <math.h>
# include <time.h>
# include <gmp.h>

/* FONCTIONS DE LECTURE/ECRITURE DE POLYNOMES DANS UN FICHIER C */
	
void parse_file(char *nom_fichier, mpz_t *polynome, unsigned long int deg){

	FILE *f;
	f = fopen(nom_fichier, "r");

	unsigned long int i = 0;
	
	if(f==NULL){
		fprintf(stderr, "Impossible d'ouvrir le fichier\n");
		exit(1);
	}
	
	while(i<=deg){
		mpz_inp_str(polynome[i], f, 10); //lit le mpz_t stocké dans le fichier, le stocke dans polynôme[i] puis passe à la ligne
		i++;
	}
		fclose(f);
}
	
	
void set_coefficients(char *nom_fichier, unsigned long int deg){
	
	FILE *f;
	f = fopen(nom_fichier, "w+");
	
	if(f==NULL){
		fprintf(stderr, "Impossible d'ouvrir le fichier\n");
		exit(1);
	}
	
	mpz_t coeff;
	mpz_init(coeff); //coeff prendra les valeurs des mpz_t générés aléatoirement
	
	gmp_randstate_t seed;
	gmp_randinit_default(seed);
	gmp_randseed_ui(seed, time(NULL)); //Paramétrage du noyau pour la génération aléatoire (cf. NOTE)	
	
	unsigned long int i;
	for(i=0; i<=deg; i++){
		mpz_rrandomb(coeff, seed, 5);
		mpz_out_str (f, 10, coeff); //mpz_out_str (FILE *stream, int base, const mpz_t op) écrit dans le flux entré en paramètre le mpz_t op dans la base indiquée
		fputs("\n",f); //retour à la ligne
	}
	
	gmp_randclear(seed);
	mpz_clear(coeff);
	
	fclose(f);
}
	
/* NOTE: 
Génération aléatoire de mpz_t. Ces mpz_t seront ensuite stockés dans un fichier qui pourra être lu par la fonction parse_file. Le nombre de mpz_t créés sera de n avec n e degré entré par l'utilisateur. 
	
Pour initialiser le noyau (seed) permettant de générer aléatoirement des nombres, il nous faut déclarer un paramètre (dit "random state parameter"), et l'initialiser par la fonction gmp_randinit_defaut.
Le noyau initial sera ensuite modifié avec la fonction gmp_randseed_ui.
	
La taille du noyau détermine le nombre de séquences de nombres aléatoires qu'il est possible de générer.
	
Une fois le noyau correctement paramétré, l'appel à mpz_rrandomb (mpz_t rop, gmp_randstate_t state, mp_bitcnt_t n) permet de générer un mpz_t aléatoire, compris entre 0 et 2^n-1 inclus, et de le stocker dans rop.
*/
	
	
/* FONCTIONS D'EVALUATION DE POLYNOMES SELON LES 2 METHODES DECRITES DANS LE FICHIER .h */
	
	
/*PROBLÈME FONCTION 1: calcul de xi par décalage ne donne pas le résultat attendu */
	
	
/* METHODE 1: q(a/(2^k)) = c0*(a/(2^k))⁰ + c1*(a/(2^k))¹ +...+ cn*(a/(2^k))^n, où n est le degré du polynôme.*/
	
void eval_poly_1(mpz_t *coeff, int a, unsigned int k, unsigned long int deg, mpz_t *res){

		mpz_set(*res, coeff[0]); //On initialise le résultat à c0*(a/(2^k))⁰ = c0
	
		// Traitement du cas deg=0
		if(deg==0) exit(0);

		int ai = a;
	
		mpz_t tmp; 
		mpz_init(tmp);
	
		unsigned long int i;
		for(i=1; i<=deg; i++){
			mpz_mul_si(tmp, coeff[i], ai);
			ai = ai*a;
			mpz_div_2exp(tmp, tmp, k*i);
			mpz_add(*res, *res, tmp);
		}

		mpz_clear(tmp);
		
}

void eval_poly_1f(mpz_t *coeff, int a, unsigned int k, unsigned long int deg, mpf_t *res){
			
		mpf_set_z(*res, coeff[0]); //On initialise le résultat à c0*(a/(2^k))⁰ = c0
	
		// Traitement du cas deg=0
		if(deg==0) exit(0);

		int ai = a;
	
		mpf_t tmp; 
		mpf_init(tmp);
		
		mpf_t coef;
		mpf_init(coef);
	
		mpf_t ai_;
		mpf_init(ai_);
	
		unsigned long int i;
		for(i=1; i<=deg; i++){
			mpf_set_z(coef, coeff[i]);
			mpf_set_si(ai_, ai);
			mpf_mul(tmp, coef, ai_);
			ai = ai*a;
			mpf_div_2exp(tmp, tmp, k*i);
			mpf_add(*res, *res, tmp);
		}

		mpf_clear(tmp);
		mpf_clear(coef);
		
}
	
/* PROBLÈME FONCTION 2: division entière du res */	


/* METHODE 2: q(a/(2^k)) = (c0*a⁰*2^(k*n) + c1*a¹*2^(k*(n-1)) +...+ cn*a^n*2^(k*(n-n)))/2^(n*k), où n est le degré du polynôme. */
	
void eval_poly_2(mpz_t *coeff, int a, unsigned int k, unsigned long int deg, mpz_t *num, int *den){
	

	int ai = a; //ai est la valeur a^i pour le coefficient ci
	int ai_decale; //ai_decale prendra la valeur de ai*2^ki
	

	mpz_t tmp; 
	mpz_init(tmp);
		
	//Traitement cas deg=0
	if(deg==0){
		mpz_set(*num, coeff[0]);
		exit(0);
	}
		
	else{
		//Premier terme du dénominateur: coeff[0]*2^(k*deg)
		mpz_mul_2exp(*num, coeff[0], k*deg);
	
		unsigned long int i;
		for(i=1; i<=deg; i++){
			ai_decale = ai<<(k*(deg-i));
			mpz_mul_si(tmp, coeff[i], ai_decale);	
			mpz_add(*num, *num, tmp);
			ai = ai*a;	
		}

		*den = 1<<deg*k;				
	}

	mpz_clear(tmp);
}
	
int main(){
	
	/* Ecriture/lecture de polynôme */

	unsigned long int deg;
	printf("Génération d'un fichier de degré n/\nn=");
	scanf("%lu", &deg);
	
	mpz_t *poly = malloc((deg+1)*sizeof(mpz_t));
	
	set_coefficients("coefficients_p.c",deg); //Ecriture de coefficients de type mpz_t aléatoires dans le fichier
	parse_file("coefficients_p.c",poly,deg); //Ecriture des coefficients dans le tableau tab
	
		
	/* Evaluation de polynôme : les valeurs de a et k seront demandées.*/
	unsigned long int i;
	for(i=0; i<=deg; i++) gmp_printf("%Zd\n", poly[i]);	
	int a;
	unsigned int k;
	printf("Calcul de q(a/(2^k)), veuillez rentrer les valeurs de:\na:\n");
	scanf("%d", &a);
	printf("\nk:\n");
	scanf("%u", &k);


	/* Méthode 1 */		
	mpz_t res1; 
	mpz_init(res1);
	eval_poly_1(poly, a, k, deg, &res1);
	gmp_printf("Resultat méthode 1: %Zd\n", res1);


	/* Méthode 1 avec des floats */		
	//mpf_t res2; 
	//mpf_init(res2);
	//eval_poly_1f(poly, a, k, deg, &res2);
	//gmp_printf("Resultat méthode 1 avec des floats: %Zd\n", res2);

	/* Méthode 2 */
	int den;
	mpz_t num;
	mpz_init(num);
	
	eval_poly_2(poly, a, k, deg, &num, &den);
	gmp_printf("Resultat méthode 2: %Zd / %d\n", num, den);
	
	free(poly);
	return 0;
}
	
	
/* PROBLÈMES RENCONTRÉS LORS DE L'ÉVALUATION DE POLYNÔMES:
	PROBLÈME COMMUN AUX 2 FONCTIONS:
	La version flotant de la premiere methode possede un probleme - segfault apparement - 
*/
