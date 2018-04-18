/* --------------------------------------- PROJET 2I013 : VERSION CÉCILE ------------------------------------------------ */

# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>
# include <string.h>
# include <math.h>
# include <time.h>
# include <gmp.h>

/* ------------------------- FONCTIONS DE LECTURE/ÉCRITURE DE POLYNÔMES DANS UN FICHIER C --------------------------------- */

//
///* Fonction de lecture des coefficients d'un polynôme depuis un fichier c (parseur) */
//
//void parse_file(char *nom_fichier, mpz_t *polynome, unsigned long int deg){
//
//	FILE *f;
//	f = fopen(nom_fichier, "r");
//
//	unsigned long int i = 0;
//
//	if(f==NULL){
//		fprintf(stderr, "Impossible d'ouvrir le fichier\n");
//		exit(1);
//	}
//
//	while(i<=deg){
//		mpz_inp_str(polynome[i], f, 10);
//		//lit le mpz_t stocké dans le fichier, le stocke dans polynôme[i] puis passe à la ligne
//		i++;
//	}
//
//	fclose(f);
//}
//
//
///* Fonction d'écriture de coefficients aléatoires d'un polynôme dans un fichier c */
//
//void set_coefficients(char *nom_fichier, unsigned long int deg){
//
//	FILE *f;
//	f = fopen(nom_fichier, "w+");
//
//	if(f==NULL){
//		fprintf(stderr, "Impossible d'ouvrir le fichier\n");
//		exit(1);
//	}
//
//	mpz_t coeff;
//	mpz_init(coeff); //coeff prendra les valeurs des mpz_t générés aléatoirement
//
//	gmp_randstate_t seed;
//	gmp_randinit_default(seed);
//	gmp_randseed_ui(seed, time(NULL)); //Paramétrage du noyau pour la génération aléatoire (cf. NOTE)
//
//	unsigned long int i;
//	for(i=0; i<=deg; i++){
//		mpz_rrandomb(coeff, seed, 5);
//		mpz_out_str (f, 10, coeff);
//    		//mpz_out_str(FILE *stream, int base, const mpz_t op) écrit dans le flux en paramètre le mpz_t op dans la base indiquée
//		fputs("\n",f); //retour à la ligne
//	}
//
//	gmp_randclear(seed);
//	mpz_clear(coeff);
//
//	fclose(f);
//}
//
/* NOTE:
Génération aléatoire de mpz_t. Ces mpz_t seront ensuite stockés dans un fichier qui pourra être lu par la fonction parse_file.
Le nombre de mpz_t créés sera de n avec n e degré entré par l'utilisateur.

Pour initialiser le noyau (seed) permettant de générer aléatoirement des nombres, il nous faut déclarer un paramètre
(dit "random state parameter"), et l'initialiser par la fonction gmp_randinit_defaut.
Le noyau initial sera ensuite modifié avec la fonction gmp_randseed_ui.

La taille du noyau détermine le nombre de séquences de nombres aléatoires qu'il est possible de générer.

Une fois le noyau correctement paramétré, l'appel à mpz_rrandomb (mpz_t rop, gmp_randstate_t state, mp_bitcnt_t n) permet de
générer un mpz_t aléatoire, compris entre 0 et 2^n-1 inclus, et de le stocker dans rop.
*/


/* -------------------------- FONCTIONS D'ÉVALUATION DE POLYNÔMES SELON LES 2 MÉTHODES DÉCRITES CI-DESSOUS -------------------*/


/* METHODE 1: q(a/(2^k)) = c0*(a/(2^k))⁰ + c1*(a/(2^k))¹ +...+ cn*(a/(2^k))^n, où n est le degré du polynôme.
   NOTE: On suppose res_num et res_den déjà initialisés avant l'appel de lafonction ci-dessous. */

void eval_poly_1(mpz_t *coeff, long long a, long long k, unsigned long int deg, mpz_t *res_num, mpz_t *res_den){

	mpz_set(*res_num, coeff[0]); //On initialise le numérateur à c0
	mpz_set_ui(*res_den, 1);

	mpz_t ai;
	mpz_init(ai); //ai prendra la valeur de a^i à chaque tour

	// Traitement du cas deg!=0
	if(deg!=0){

		unsigned long int i;
		mpz_set_ui(ai, a);


		for(i=1; i<=deg; i++){

			//*res_num = (*res_num)*2^k + ci*a^i, avec *res_num = c0 au premier tour
			mpz_mul_2exp(*res_num, *res_num, k);
			mpz_addmul(*res_num, coeff[i], ai);
			mpz_mul_ui(ai, ai, a);
		}

		//La valeur du dénominateur sera toujours 2^(k*deg)
		mpz_mul_2exp(*res_den, *res_den, k*deg);

		//RÉDUCTION DE LA FRACTION

		mpz_t gcd; mpz_init(gcd);
		mpz_gcd(gcd, *res_num, *res_den);

		//On divise *res_num et *res_den par leur gcd et on réduit ainsi la fraction
		mpz_divexact(*res_num, *res_num, gcd);
		mpz_divexact(*res_den, *res_den, gcd);

		/*NOTE: mpz_divexact(mpz_t a, const mpz_t , const mpz_t d)
			sets q to n/p. It produces correct results only when it is known in advance that n divides p. */
	}
	mpz_clear(ai);
}


/* METHODE 2: q(a/(2^k)) = (c0*a⁰*2^(k*n) + c1*a¹*2^(k*(n-1)) +...+ cn*a^n*2^(k*(n-n)))/2^(n*k), où n est le degré du polynôme.
   NOTE: On suppose res_num et res_den déjà initialisés avant l'appel de lafonction ci-dessous. */

void eval_poly_2(mpz_t *coeff, long long a, long long k, unsigned long int deg, mpz_t *res_num, mpz_t *res_den){

	long long ki = k; //ki est la valeur de la puissance de 2 au numérateur pour le coefficient ci
	mpz_t ai; mpz_init(ai); //ai est la valeur a^i pour le coefficient ci
	mpz_set_ui(ai, 1);
	mpz_t ai_decale; mpz_init(ai_decale); //ai_decale prendra la valeur de ai*2^ki

	mpz_t tmp;
	mpz_init(tmp); /* Il n'existe pas de fonction mpz_add_mul pour les long long signés.
			  On fera donc l'addition et la multiplication à part. */

	//Traitement cas deg=0
	if(deg==0){
		mpz_set(*res_num, coeff[0]);
		mpz_set_ui(*res_den, 1);
	}

	else{
		//Premier terme du dénominateur: coeff[0]*2^(k*deg)
		mpz_mul_2exp(*res_num, coeff[0], k*deg);

		unsigned long int i;
		for(i=1; i<=deg; i++){
			mpz_mul_ui(ai, ai, a);
			ki = k*(deg-i);
			mpz_mul_2exp(ai_decale, ai, ki);
			mpz_mul(tmp, coeff[i], ai_decale);
			mpz_add(*res_num, *res_num, tmp);
		}
			mpz_set_ui(*res_den, 1);
			mpz_mul_2exp(*res_den, *res_den, k*deg);

		//RÉDUCTION DE LA FRACTION

		mpz_t gcd; mpz_init(gcd);
		mpz_gcd(gcd, *res_num, *res_den);

		//On divise *res_num et *res_den par leur gcd et on réduit ainsi la fraction
		mpz_divexact(*res_num, *res_num, gcd);
		mpz_divexact(*res_den, *res_den, gcd);

		/*NOTE: mpz_divexact(mpz_t a, const mpz_t , const mpz_t d)
			sets q to n/p. It produces correct results only when it is known in advance that n divides p. */
	}
	mpz_clear(tmp);
}


/* --------------------------------------- IMPLÉMENTATION DE LA MÉTHODE DE HORNER ------------------------------------------- */

/* La méthode de Horner pour un polynôme P(X) = c0 + c1*X + c2*X² + ... + cn*X^n donne :
              P(X) = c0 + X(c1 + X(c2 + ... + X(cn-2 + X(cn-1 + cn*X))...)).
Appliquée dans notre cas, ie. pour un polynôme évalué en X = a/2^k, cela donne:
              P(X) = c0 + (c1 + ... (cn-2 + (cn-1 + (cn*a)/2^k)*a/2^k)*a/2^k) ... )*a/2^k

Soit, en mettant tout sur le même dénominateur:
	a) au numérateur: c0*2^(nk) + c1*a*2^(n-1)k + ... + cn-3*a*2^(3k) + cn-2*a*2^(2k) + cn-1*a*2^k + cn*a
	b) au dénominateur: 2^(nk)
*/

void evaluateHorner(mpz_t *coeff, long long a, long long k, unsigned long int deg, mpz_t *res_num, mpz_t *res_den){

	//printf("\n\n---EVALUATION DE HORNER---\n");

	mpz_set_ui(*res_num, 0);
	mpz_set_ui(*res_den, 1);

	unsigned long int i;

	mpz_t tmp; mpz_init(tmp);
  	/*Pour chaque calcul de num_i = coeff[i]*2^k + a*num_i (tour de boucle i), tmp prendra la valeur de coeff[i]*2^k*/

  for(i=deg; i>0; i--){
		mpz_mul_ui(*res_num, *res_num, a); //Calcul de a*(*res_num)
		//Calcul de coeff[i]*2^k
		mpz_mul_2exp(tmp, coeff[i], (deg-i)*k);
		//*res_num = coeff[i]*2^k + a*(*res_num)
		mpz_add(*res_num, *res_num, tmp);
	}

	//Traitement du coeff c0
	mpz_mul_ui(*res_num, *res_num, a); //Calcul de a*(*res_num)
	mpz_mul_2exp(tmp, coeff[0], deg*k);
	mpz_add(*res_num, *res_num, tmp);

	//La valeur du dénominateur sera toujours 2^(k*deg)
	mpz_mul_2exp(*res_den, *res_den, k*deg);

	//gmp_printf("Horner avant réduction: num = %Zd, den = %Zd\n", *res_num, *res_den);

	//RÉDUCTION DE LA FRACTION

	mpz_t gcd; mpz_init(gcd);
	mpz_gcd(gcd, *res_num, *res_den);

	//On divise *res_num et *res_den par leur gcd et on réduit ainsi la fraction
	mpz_divexact(*res_num, *res_num, gcd);
	mpz_divexact(*res_den, *res_den, gcd);

	/*NOTE: mpz_divexact(mpz_t a, const mpz_t , const mpz_t d)
		sets q to n/p. It produces correct results only when it is known in advance that n divides p. */

	mpz_clear(tmp);
}


/* ------------------------------------------------- TEST MAIN ------------------------------------------------------------ */
//
//
//int main(){
//
//	/* Ecriture/lecture de polynôme */
//
//	unsigned long int deg;
//	printf("Génération d'un fichier de degré n/\nn=");
//	scanf("%lu", &deg);
//
//	mpz_t *tab = malloc((deg+1)*sizeof(mpz_t));
//
//	set_coefficients("coefficients_p.c",deg); //Ecriture de coefficients de type mpz_t aléatoires dans le fichier
//	parse_file("coefficients_p.c",tab,deg); //Ecriture des coefficients dans le tableau tab
//
//
//	/* Evaluation de polynôme : les valeurs de a et k seront demandées.*/
//
//	unsigned long int i;
//	for(i=0; i<=deg; i++) gmp_printf("%Zd\n", tab[i]);
//	long long a, k;
//	printf("Calcul de q(a/(2^k)), veuillez rentrer les valeurs de:\na:\n");
//	scanf("%llu", &a);
//	printf("\nk:\n");
//	scanf("%llu", &k);
//
//	printf("\n\n\t----------------- EVALUATION DES POLYNÔMES -----------------\n");
//
//	/* Méthode 1 */
//
//	mpz_t res_num1, res_den1;
//	mpz_init(res_num1);
//	mpz_init(res_den1);
//
//	eval_poly_1(tab, a, k, deg, &res_num1, &res_den1);
//	gmp_printf("\n\nResultat méthode 1: num = %Zd\t den = %Zd\n", res_num1, res_den1);
//
//	/* Méthode 2 */
//
//	mpz_t res_num2, res_den2;
//	mpz_init(res_num2);
//	mpz_init(res_den2);
//
//	eval_poly_2(tab, a, k, deg, &res_num2, &res_den2);
//	gmp_printf("\n\nResultat méthode 2: num = %Zd\t den = %Zd\n", res_num2, res_den2);
//
//	/* TEST MÉTHODE DE HORNER */
//
//	mpz_t res_num3, res_den3;
//	mpz_init(res_num3);
//	mpz_init(res_den3);
//
//	evaluateHorner(tab, a, k, deg, &res_num3, &res_den3);
//	gmp_printf("\n\nResultat Horner: num = %Zd\t den = %Zd\n", res_num3, res_den3);
//
//
//	mpz_clears(res_num1, res_den1, res_num2, res_den2, res_num3, res_den3, NULL);
//
//	free(tab);
//	return 0;
//}
//
