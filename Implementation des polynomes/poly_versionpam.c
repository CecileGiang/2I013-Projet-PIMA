
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


/* FONCTIONS D EVALUATION DE POLYNOMES SELON LES 2 METHODES DECRITES DANS LE FICHIER .h */


/*PROBLEME FONCTION 1: calcul de xi par decalage ne donne pas le resultat attendu */


/* METHODE 1: q(a/(2^k)) = c0*(a/(2^k))⁰ + c1*(a/(2^k))¹ +...+ cn*(a/(2^k))^n, ou n est le degre du polynome*/

//Cette premiere fonction ne marche pas adequadement parce qu'elle n'est pas precise a cause de la facon dont la division est traitee, pourtant, elle nous sert pour avoir une approximation du resultat et pouvoir comparer

void eval_poly_1(mpz_t *coeff, int a, unsigned int k, unsigned long int deg, mpz_t *res){

		mpz_set(*res, coeff[0]); //On initialise le résultat à c0*(a/(2^k))⁰ = c0

		// Traitement du cas deg=0
		if(deg==0) exit(0);

		int ai = a;

		mpz_t tmp;
		mpz_init(tmp);

		unsigned long int i;
		for(i=1; i<=deg; i++){
			mpz_mul_si(tmp, coeff[i], ai); //on multiplie a**i par la valeur du i-eme coefficient
			ai = ai*a;			//on augmente la puissance de a par 1
			mpz_div_2exp(tmp, tmp, k*i);	//on divisse le numerateur par la i-eme puissance de 2**k
			mpz_add(*res, *res, tmp);  	//on fait la somme des deux termes
		}

		mpz_clear(tmp);

}


// Cette fonction ne marche pas: les flotants marchent differement sur gmp

void eval_poly_1f(mpz_t *coeff, int a, unsigned int k, unsigned long int deg, mpf_t *res){

		mpf_set_z(*res, coeff[0]); //On initialise le resultat a c0*(a/(2^k))⁰ = c0

		// Traitement du cas deg=0
		if(deg==0) return;

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

// Ceci est la fonction qui corresponde a la methode 1

void eval_poly_1bis(mpz_t *coeff, int a, unsigned int k, unsigned long int deg, mpz_t *num, mpz_t *den){
	
	//creation de variables mpz_t temporels pour stocker les additions et produits
	mpz_t tmp1;
	mpz_t tmp2;
	mpz_init(tmp1);
	mpz_init(tmp2);

	// creation de variables qui simulent les deux numerateurs pour faire une somme deux a deux
	mpz_t num1;
	mpz_t num2;
	mpz_init(num1);
	mpz_init(num2);

	unsigned long int i;

	long int ai = 1;
	long int ai_decale = a;

	mpz_t den1;
	mpz_init(den1);
	mpz_set_si(den1, 1);	//le premier terme a comme denominateur 1, qui est (2**k)**0

	mpz_t den2;
	mpz_init(den2);
	mpz_set_si(den2, 1);
	mpz_mul_2exp(den2,den2, k); //le deuxieme terme a comme denominateur (2**k)**1

	mpz_mul_si(num1, coeff[0], ai);
	mpz_mul_si(num2, coeff[1], ai_decale);

	if (deg == 0){
		mpz_set(*num, num1);
		mpz_set(*den, den1);
		return;
	}

	mpz_mul(tmp1, num1, den2);
	mpz_mul(tmp2, num2, den1);
	mpz_add(*num, tmp1, tmp2);
	mpz_mul(*den, den1, den2);

	for (i = 2; i <= deg; i++){
		mpz_set(num1, *num);
		mpz_set(den1, *den);
		ai_decale = ai_decale * a;
		mpz_mul_si(tmp1, coeff[i], ai_decale);
		mpz_set(num2, tmp1);
		mpz_mul_2exp(den2, den1, k);
		mpz_mul(tmp1, num1, den2);
		mpz_mul(tmp2, num2, den1);
		mpz_add(*num, tmp1, tmp2);
		mpz_mul(*den, den1, den2);
	}
	mpz_clear(num1);
	mpz_clear(num2);
	mpz_clear(tmp1);
	mpz_clear(tmp2);
	mpz_clear(den1);
	mpz_clear(den2);
}

/* PROBLÈME FONCTION 2: division entière du res */


/* METHODE 2: q(a/(2^k)) = (c0*a⁰*2^(k*n) + c1*a¹*2^(k*(n-1)) +...+ cn*a^n*2^(k*(n-n)))/2^(n*k), où n est le degré du polynôme. */

void eval_poly_2(mpz_t *coeff, int a, unsigned int k, unsigned long int deg, mpz_t *num, mpz_t *den){


	int ai = a; //ai est la valeur a^i pour le coefficient ci
	int ai_decale; //ai_decale prendra la valeur de ai*2^ki


	mpz_t tmp;
	mpz_init(tmp);
	mpz_set_si(*den, 1);
	//Traitement cas deg=0
	if(deg==0){
		mpz_set(*num, coeff[0]);
		return;
	}

	else{
		//Premier terme du dénominateur: coeff[0]*2^(k*deg)
		mpz_mul_2exp(*num, coeff[0], k*deg);

		unsigned long int i;
		for(i=1; i<=deg; i++){
			ai_decale = ai*pow(2,k*(deg-i));
			mpz_mul_si(tmp, coeff[i], ai_decale);
			mpz_add(*num, *num, tmp);
			ai = ai*a;
		}

		mpz_mul_2exp(*den, *den, k*deg);
	}

	mpz_clear(tmp);
}

//Methode de Horner

void eval_poly_horner(mpz_t *coeff, int a, unsigned int k, unsigned long int deg, mpz_t *num, mpz_t *den){
    mpz_t tmp1;
    mpz_init(tmp1);
    mpz_t tmp2;
    mpz_init(tmp2);

    mpz_set_si(*den, 1);
    long int d = deg;
    long int i;
    
    for (i = d; i >= 0; i--){
	 if(i!=d) mpz_mul_2exp(*den, *den, k);
      	mpz_mul_si(tmp1, *num, a);
	mpz_mul(tmp2, coeff[i], *den);
 	mpz_add(*num, tmp2 , tmp1);
	//gmp_printf("tmp1: %Zd\n tmp3: %Zd\n tmp2: %Zd\n num:%Zd\n", tmp1, tmp3, tmp2, *num);
  		
	
  }
      mpz_clear(tmp1);
      mpz_clear(tmp2);
}

void reduire_fraction(mpz_t *den, mpz_t *num){
	
	mpz_t num1;
	mpz_init(num1);
	mpz_set(num1, *num);

	mpz_t reste_num;
	mpz_init(reste_num);
	
	do{
		mpz_tdiv_r_2exp(reste_num, num1, 1); //reste_num stockera le reste de la division euclidienne de num1 par 2**1
	}
	while (mpz_cmp(reste_num, 0) == 0); // mpz_cmp(mpz_t op1, double op2) renvoie 0 lorsque op1 = op2
	mpz_mul_2exp(num1, num1, 1);
	mpz_add_ui(*num, num1, 1);

}

int main(){

	/* Ecriture/lecture de polynôme */

	unsigned long int deg;
	printf("Génération d'un fichier de degré n/\nn=");
	scanf("%lu", &deg);

	mpz_t *poly = malloc((deg+1)*sizeof(mpz_t));

	set_coefficients("coefficients_p.c",deg); //Ecriture de coefficients de type mpz_t aléatoires dans le fichier
	parse_file("coefficients_p.c",poly,deg); //Ecriture des coefficients dans le tableau tab


	/* Evaluation de polynôme : les valeurs de a et k seront demandees*/
	unsigned long int i;
	for(i=0; i<=deg; i++) {
		gmp_printf("%Zd\n", poly[i]);
	}
	int a;
	unsigned int k;
	printf("Calcul de q(a/(2^k)), veuillez rentrer les valeurs de:\na:\n");
	scanf("%d", &a);
	printf("\nk:\n");
	scanf("%u", &k);


	/* Methode 1 */
	mpz_t res1;
	mpz_init(res1);
	eval_poly_1(poly, a, k, deg, &res1);
	gmp_printf("Resultat méthode 1: %Zd\n", res1);

	mpz_clear(res1);

	/* Methode 1 avec des floats */
	//mpf_t res2;
	//mpf_init(res2);
	//eval_poly_1f(poly, a, k, deg, &res2);
	//gmp_printf("Resultat methode 1 avec des floats: %Zd\n", res2);

	/* Methode 1-bis*/
	mpz_t den1;
	mpz_t num1;
	mpz_init(num1);
	mpz_init(den1);

	eval_poly_1bis(poly, a, k, deg, &num1, &den1);
	gmp_printf("Resultat méthode 1-bis: %Zd / %d\n", num1, den1);

	mpz_clear(num1);
	mpz_clear(den1);

	/* Methode 2 */
	mpz_t den2;
	mpz_t num2;
	mpz_init(num2);
	mpz_init(den2);

	eval_poly_2(poly, a, k, deg, &num2, &den2);
	gmp_printf("Resultat méthode 2: %Zd / %Zd\n", num2, den2);

	mpz_clear(num2);
	mpz_clear(den2);

	//Methode de Horner
	mpz_t den3;
	mpz_t num3;
	mpz_init(num3);
	mpz_init(den3);

	eval_poly_horner(poly, a, k, deg, &num3, &den3);
	gmp_printf("Resultat methode 3 (Horner): %Zd / %Zd\n", num3, den3);
	mpz_clear(num3);
	mpz_clear(den3);

	free(poly);
	return 0;
}
	
/* PROBLÈMES RENCONTRÉS LORS DE L'ÉVALUATION DE POLYNÔMES:
	La version flotant de la premiere methode possede un probleme - segfault apparement - 
*/ 
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


/* FONCTIONS D EVALUATION DE POLYNOMES SELON LES 2 METHODES DECRITES DANS LE FICHIER .h */


/*PROBLEME FONCTION 1: calcul de xi par decalage ne donne pas le resultat attendu */


/* METHODE 1: q(a/(2^k)) = c0*(a/(2^k))⁰ + c1*(a/(2^k))¹ +...+ cn*(a/(2^k))^n, ou n est le degre du polynome*/

//Cette premiere fonction ne marche pas adequadement parce qu'elle n'est pas precise a cause de la facon dont la division est traitee, pourtant, elle nous sert pour avoir une approximation du resultat et pouvoir comparer

void eval_poly_1(mpz_t *coeff, int a, unsigned int k, unsigned long int deg, mpz_t *res){

		mpz_set(*res, coeff[0]); //On initialise le résultat à c0*(a/(2^k))⁰ = c0

		// Traitement du cas deg=0
		if(deg==0) exit(0);

		int ai = a;

		mpz_t tmp;
		mpz_init(tmp);

		unsigned long int i;
		for(i=1; i<=deg; i++){
			mpz_mul_si(tmp, coeff[i], ai); //on multiplie a**i par la valeur du i-eme coefficient
			ai = ai*a;			//on augmente la puissance de a par 1
			mpz_div_2exp(tmp, tmp, k*i);	//on divisse le numerateur par la i-eme puissance de 2**k
			mpz_add(*res, *res, tmp);  	//on fait la somme des deux termes
		}

		mpz_clear(tmp);

}


// Cette fonction ne marche pas: les flotants marchent differement sur gmp

void eval_poly_1f(mpz_t *coeff, int a, unsigned int k, unsigned long int deg, mpf_t *res){

		mpf_set_z(*res, coeff[0]); //On initialise le resultat a c0*(a/(2^k))⁰ = c0

		// Traitement du cas deg=0
		if(deg==0) return;

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

// Ceci est la fonction qui corresponde a la methode 1

void eval_poly_1bis(mpz_t *coeff, int a, unsigned int k, unsigned long int deg, mpz_t *num, mpz_t *den){
	
	//creation de variables mpz_t temporels pour stocker les additions et produits
	mpz_t tmp1;
	mpz_t tmp2;
	mpz_init(tmp1);
	mpz_init(tmp2);

	// creation de variables qui simulent les deux numerateurs pour faire une somme deux a deux
	mpz_t num1;
	mpz_t num2;
	mpz_init(num1);
	mpz_init(num2);

	unsigned long int i;

	long int ai = 1;
	long int ai_decale = a;

	mpz_t den1;
	mpz_init(den1);
	mpz_set_si(den1, 1);	//le premier terme a comme denominateur 1, qui est (2**k)**0

	mpz_t den2;
	mpz_init(den2);
	mpz_set_si(den2, 1);
	mpz_mul_2exp(den2,den2, k); //le deuxieme terme a comme denominateur (2**k)**1

	mpz_mul_si(num1, coeff[0], ai);
	mpz_mul_si(num2, coeff[1], ai_decale);

	if (deg == 0){
		mpz_set(*num, num1);
		mpz_set(*den, den1);
		return;
	}

	mpz_mul(tmp1, num1, den2);
	mpz_mul(tmp2, num2, den1);
	mpz_add(*num, tmp1, tmp2);
	mpz_mul(*den, den1, den2);

	for (i = 2; i <= deg; i++){
		mpz_set(num1, *num);
		mpz_set(den1, *den);
		ai_decale = ai_decale * a;
		mpz_mul_si(num2, coeff[i], ai_decale);
		mpz_mul_2exp(den2, den1, k);
		mpz_mul(tmp1, num1, den2);
		mpz_mul(tmp2, num2, den1);
		mpz_add(*num, tmp1, tmp2);
		mpz_mul(*den, den1, den2);
	}
	mpz_clear(num1);
	mpz_clear(num2);
	mpz_clear(tmp1);
	mpz_clear(tmp2);
	mpz_clear(den1);
	mpz_clear(den2);
}

/* PROBLÈME FONCTION 2: division entière du res */


/* METHODE 2: q(a/(2^k)) = (c0*a⁰*2^(k*n) + c1*a¹*2^(k*(n-1)) +...+ cn*a^n*2^(k*(n-n)))/2^(n*k), où n est le degré du polynôme. */

void eval_poly_2(mpz_t *coeff, int a, unsigned int k, unsigned long int deg, mpz_t *num, mpz_t *den){


	int ai = a; //ai est la valeur a^i pour le coefficient ci
	int ai_decale; //ai_decale prendra la valeur de ai*2^ki


	mpz_t tmp;
	mpz_init(tmp);
	mpz_set_si(*den, 1);
	//Traitement cas deg=0
	if(deg==0){
		mpz_set(*num, coeff[0]);
		return;
	}

	else{
		//Premier terme du dénominateur: coeff[0]*2^(k*deg)
		mpz_mul_2exp(*num, coeff[0], k*deg);

		unsigned long int i;
		for(i=1; i<=deg; i++){
			ai_decale = ai*pow(2,k*(deg-i));
			mpz_mul_si(tmp, coeff[i], ai_decale);
			mpz_add(*num, *num, tmp);
			ai = ai*a;
		}

		mpz_mul_2exp(*den, *den, k*deg);
	}

	mpz_clear(tmp);
}

//Methode de Horner

void eval_poly_horner(mpz_t *coeff, int a, unsigned int k, unsigned long int deg, mpz_t *num, mpz_t *den){
    mpz_t tmp1;
    mpz_init(tmp1);
    mpz_t tmp2;
    mpz_init(tmp2);

    mpz_set_si(*den, 1);
    long int d = deg;
    long int i;
    
    for (i = d; i >= 0; i--){
	 if(i!=d) mpz_mul_2exp(*den, *den, k);
      	mpz_mul_si(tmp1, *num, a);
	mpz_mul(tmp2, coeff[i], *den);
 	mpz_add(*num, tmp2 , tmp1);
	//gmp_printf("tmp1: %Zd\n tmp3: %Zd\n tmp2: %Zd\n num:%Zd\n", tmp1, tmp3, tmp2, *num);
  		
	
  }
      mpz_clear(tmp1);
      mpz_clear(tmp2);
}

void reduire_fraction(mpz_t *den, mpz_t *num){
	
	mpz_t num1;
	mpz_init(num1);
	mpz_set(num1, *num);

	mpz_t reste_num;
	mpz_init(reste_num);
	
	do{
		mpz_tdiv_r_2exp(reste_num, num1, 1); //reste_num stockera le reste de la division euclidienne de num1 par 2**1
		mpz_div_2exp(*den, *den, 1);
	}
	while (mpz_cmp(reste_num, 0) == 0); // mpz_cmp(mpz_t op1, double op2) renvoie 0 lorsque op1 = op2
	mpz_mul_2exp(num1, num1, 1);
	mpz_add_ui(*num, num1, 1);
	mpz_mul_2exp(*den, *den, 1); 

}

int main(){

	/* Ecriture/lecture de polynôme */

	unsigned long int deg;
	printf("Génération d'un fichier de degré n/\nn=");
	scanf("%lu", &deg);

	mpz_t *poly = malloc((deg+1)*sizeof(mpz_t));

	set_coefficients("coefficients_p.c",deg); //Ecriture de coefficients de type mpz_t aléatoires dans le fichier
	parse_file("coefficients_p.c",poly,deg); //Ecriture des coefficients dans le tableau tab


	/* Evaluation de polynôme : les valeurs de a et k seront demandees*/
	unsigned long int i;
	for(i=0; i<=deg; i++) {
		gmp_printf("%Zd\n", poly[i]);
	}
	int a;
	unsigned int k;
	printf("Calcul de q(a/(2^k)), veuillez rentrer les valeurs de:\na:\n");
	scanf("%d", &a);
	printf("\nk:\n");
	scanf("%u", &k);


	/* Methode 1 */
	mpz_t res1;
	mpz_init(res1);
	eval_poly_1(poly, a, k, deg, &res1);
	gmp_printf("Resultat méthode 1: %Zd\n", res1);

	mpz_clear(res1);

	/* Methode 1 avec des floats */
	//mpf_t res2;
	//mpf_init(res2);
	//eval_poly_1f(poly, a, k, deg, &res2);
	//gmp_printf("Resultat methode 1 avec des floats: %Zd\n", res2);

	/* Methode 1-bis*/
	mpz_t den1;
	mpz_t num1;
	mpz_init(num1);
	mpz_init(den1);

	eval_poly_1bis(poly, a, k, deg, &num1, &den1);
	gmp_printf("Resultat méthode 1-bis: %Zd / %d\n", num1, den1);

	mpz_clear(num1);
	mpz_clear(den1);

	/* Methode 2 */
	mpz_t den2;
	mpz_t num2;
	mpz_init(num2);
	mpz_init(den2);

	eval_poly_2(poly, a, k, deg, &num2, &den2);
	gmp_printf("Resultat méthode 2: %Zd / %Zd\n", num2, den2);

	mpz_clear(num2);
	mpz_clear(den2);

	//Methode de Horner
	mpz_t den3;
	mpz_t num3;
	mpz_init(num3);
	mpz_init(den3);

	eval_poly_horner(poly, a, k, deg, &num3, &den3);
	gmp_printf("Resultat methode 3 (Horner): %Zd / %Zd\n", num3, den3);
	mpz_clear(num3);
	mpz_clear(den3);

	free(poly);
	return 0;
}
	
/* PROBLÈMES RENCONTRÉS LORS DE L'ÉVALUATION DE POLYNÔMES:
	La version flotant de la premiere methode possede un probleme - segfault apparement - 
*/
