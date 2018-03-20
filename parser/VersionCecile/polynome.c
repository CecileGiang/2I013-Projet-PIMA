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

void eval_poly_1(mpz_t *coeff, long long a, long long k, unsigned long int deg, mpz_t *res){

	printf("a = %lld, k = %lld, a>>k = %lf\n", a, k, ((long double)a)>>(long double)k);
	long double x = a>>k;
	gmp_printf("x = %lf\n", x);

	//A enlever
	printf("k = %lld\n", k);
	gmp_printf("2^k = %lu\n", pow(2,2));

	mpz_set(*res, coeff[0]); //On initialise le résultat à c0*(a/(2^k))⁰ = c0

	// Traitement du cas deg=0
	if(deg==0) exit(0);
	
	
	long long xi = 1; //xi prendra les valeurs de x^i = (a/(2^k))^i
	mpz_t tmp; mpz_init(tmp); /* Il n'existe pas de fonction mpz_add_mul pour les long long signés.
		   		   On fera donc l'addition et la multiplication à part. */
	//A enlever
	gmp_printf("0: %Zd\n",*res);

	unsigned long int i;
	printf("%lld\n", xi);
	for(i=1; i<=deg; i++){
		gmp_printf("%lu: res = %Zd\n", i, *res);
		xi = xi*x;
		printf("xi = %lld\n", xi);
		mpz_mul_si(tmp, coeff[i], xi);	//Bizzarement c'est la multiplication de c[i]*xi qui plante, pour de grandes valeurs
						//Pourtant tmp est bien de type mpz_t
		mpz_add(*res, *res, tmp);
		//A enlever
		gmp_printf("c[i] = %Zd\n", coeff[i]);
		gmp_printf("c[i]*xi = %Zd\n", tmp);
		gmp_printf("res = %Zd\n", *res);
	}

	mpz_clear(tmp);
}

/* PROBLÈME FONCTION 2: division entière du res */


/*METHODE 2: q(a/(2^k)) = (c0*a⁰*2^(k*n) + c1*a¹*2^(k*(n-1)) +...+ cn*a^n*2^(k*(n-n)))/2^(n*k), où n est le degré du polynôme. */

void eval_poly_2(mpz_t *coeff, long long a, long long k, unsigned long int deg, mpz_t *res){

	long long ki = k; //ki est la valeur de la puissance de 2 au numérateur pour le coefficient ci
	long long ai = 1; //ai est la valeur a^i pour le coefficient ci
	long long ai_decale; //ai_decale prendra la valeur de ai*2^ki

	mpz_t tmp; mpz_init(tmp); /* Il n'existe pas de fonction mpz_add_mul pour les long long signés.
		   		   On fera donc l'addition et la multiplication à part. */
	
	//Traitement cas deg=0
	if(deg==0){
		mpz_set(*res, coeff[0]);
		exit(0);
	}
	
	else{

	//Premier terme du dénominateur: coeff[0]*2^(k*deg)
		mpz_mul_2exp(*res, coeff[0], k*deg);

		unsigned long int i;
		for(i=1; i<=deg; i++){
			ai *= a;
			ki = k*(deg-i);
			ai_decale = ai<<ki;
			mpz_mul_si(tmp, coeff[i], ai_decale);						
			mpz_add(*res, *res, tmp);
		}

		mpz_tdiv_q_2exp(*res, *res, k*deg);

/*	unsigned long int res_ui; //valeur de *res en ui. A changer parce que la plage de valeurs n'est pas la même

	res_ui = mpz_get_ui(*res);
	res_ui >>= k*deg;

	mpz_set_ui(*res, res_ui);*/

	}

	mpz_clear(tmp);
}


int main(){

	/* Ecriture/lecture de polynôme */

	unsigned long int deg;
	printf("Génération d'un fichier de degré n/\nn=");
	scanf("%lu", &deg);

	mpz_t *tab = malloc((deg+1)*sizeof(mpz_t));

	set_coefficients("coefficients_p.c",deg); //Ecriture de coefficients de type mpz_t aléatoires dans le fichier
	parse_file("coefficients_p.c",tab,deg); //Ecriture des coefficients dans le tableau tab

	
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
	
	eval_poly_1(tab, a, k, deg, &res);
	gmp_printf("Resultat méthode 1: %Zd\n", res);

	/* Méthode 2 */

	eval_poly_2(tab, a, k, deg, &res1);
	gmp_printf("Resultat méthode 2: %Zd\n", res1);


	free(tab);
	return 0;
}


/* PROBLÈMES RENCONTRÉS LORS DE L'ÉVALUATION DE POLYNÔMES:

MÉTHODE 1:
* x désigne la valeur de a/2^k. Techniquement xi peut être flottant, mais les opérations avec des mpz_t s'effectuent obligatoirement avec des entiers. Cela rend le calcul faux.
* De plus il semble y avoir un problème avec les décalages pour xi trop petit.
* mais si l'on essaie de changer tous les mpz_t en mpf_t, on se retrouve avec un autre problème: les opérations mpz_t ne possèdent pas toujours d'équivalent en mpf_t.
* erreur d'approximation ? à vérifier

MÉTHODE 2:
* Division entière 

PROBLÈME COMMUN AUX 2 FONCTIONS:
* Les calculs ne fonctionnent pas pour de grandes valeurs de a et/ou de k (le résultat est pourtant bien donné en mpz_t): effet de bord.
*/

