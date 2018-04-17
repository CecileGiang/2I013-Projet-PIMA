# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>
# include <string.h>
# include <math.h>
# include <time.h>
# include <gmp.h>

/* _______________FONCTIONS DE LECTURE/ECRITURE DE POLYNOMES DANS UN FICHIER C___________ */

/* Cette fonction nous permet de lire les coefficients de type mpz_t ecrits dans un fichier sur la forme d une liste. 
Nous lisons donc chaque ligne, ensuite nous stockons le coefficient dans un tableau qui representera notre polynome et puis nous faisons un retour a la ligne. 
Nous faisons cela jusqu a avoir n+1 coefficients, n etant le degre du polynome et qui a ete saisi precedemment.
L indice du tableau polynome correspondra bien au degre du terme, ainsi le premier terme (d indice 0) sera le terme constant, le deuxieme terme (d indice 1) celui de degre 1, et ainsi de suite. */

void parse_file(char *nom_fichier, mpz_t *polynome, unsigned long int deg){

	FILE *f;
	f = fopen(nom_fichier, "r");

	unsigned long int i = 0;

	if(f==NULL){
		fprintf(stderr, "Impossible d'ouvrir le fichier\n");
		exit(1);
	}

	while(i<=deg){
		mpz_inp_str(polynome[i], f, 10); //lit le mpz_t stocké dans le fichier, le stocke dans polynome[i] puis passe a la ligne
		i++;
	}
		fclose(f);
}

/* __________GENERATION ALEATOIRE DE MPZ_T___________ */

/*

Cette fonction ecrit des coefficients entiers de type mpz sur un fichier donne en argument, qui pourrq etre lu plus tard par la fonction parse_file.
Le nombre de coefficients sera de n+1 avec n le degre donne par l utilisateur.
Ceci se fait en implantantant un etat random a l aide de l algorithme de Mersenne Twister. 
Ce noyau donne plus tard des valeurs aleatoires (entre 0 et 2^(n-1), n etant le troisieme argument de la fonction)
a chaque coefficient.
Comme nous souhaitons aussi avoir des nombres negatifs, nous les generons aleatoirement, c est-a-dire, nous creons une variable signe,
qui determinera si un nombre est positif ou negatif avec une equiprobabilite.

*/

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
		mpz_rrandomb(coeff, state, 120); //coeff prend une valeur random entre 0 et 2^(n-1) grace a l algorithme pris par state
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


/*      Commentaires: nous devrions avoir des nombres aleatoires entre 0 et 2^(n-1), 
	pourtant nous trouvons que des nombres avec le meme nombre de chiffres a chaque fois  */
}

/* ______FONCTION AUXILIAIRE: REDUCTION D UNE FRACTION_________ */

/* 

Nous utiliserons cette fonction pour avoir la fraction qui donne comme resultat l evaluation du polynome, sur sa forme reduite.
Nous profitons le fait que le denominateur est toujours une puissance de 2, donc, tant que le reste de la division du numerateur par 2 soit 0,
nous en deduisons que le numerateur est pair, est donc nous pouvons le diviser par 2, ainsi que le denominateur.
Comme pour faire ce teste nous aurons deja divise, nous avons une division de plus. 
Ce pourquoi a la fin nous ajoutons 1 au resultat et nous multiplions par 2

*/

void reduire_fraction(mpz_t *den, mpz_t *num){
	
	mpz_t num1;
	mpz_init(num1);
	mpz_set(num1, *num);

	mpz_t reste_num;
	mpz_init(reste_num);
	
	do{
		mpz_tdiv_r_2exp(reste_num, num1, 1); //reste_num stockera le reste de la division euclidienne de num1 par 2^1
		mpz_div_2exp(*den, *den, 1);
	}
	while (mpz_cmp(reste_num, 0) == 0); // mpz_cmp(mpz_t op1, double op2) renvoie 0 lorsque op1 = op2
	mpz_mul_2exp(num1, num1, 1);
	mpz_add_ui(*num, num1, 1);
	mpz_mul_2exp(*den, *den, 1); 

}


/* __________FONCTIONS D EVALUATION DE POLYNOMES__________ */

/* METHODE 1: 

Soit Q un polynome a coefficients entiers notes c0, c1, ... , cn et n son degre:

Q(a/(2^k)) = c0*(a/(2^k))⁰ + c1*(a/(2^k))¹ +...+ cn*(a/(2^k))^n */

/* 

Avec cette premiere fonction, nous obtenons un resultat qui n est pas precis a cause de la facon dont la division est traitee (c est-a-dire division entiere), 
pourtant elle nous sert pour avoir une approximation du resultat et ainsi avoir un point de comparaison.
Elle marche de la maniere suivant: nous calculons a chaque iteration le produit du coefficient par la valeur d evaluation, 
qui peut bien ne pas etre un entier. Pourtant, nous traitons cette valeur comme un entier et nous rajoutons ce produit a une variable qui representera le resultat. 

*/

void eval_poly_1(mpz_t *coeff, int a, unsigned int k, unsigned long int deg, mpz_t *res){

		mpz_set(*res, coeff[0]); //On initialise le résultat à c0*(a/(2^k))⁰ = c0

		// Traitement du cas deg=0
		if(deg==0) exit(0);

		int ai = a;

		mpz_t tmp;
		mpz_init(tmp);

		unsigned long int i;
		for(i=1; i<=deg; i++){
			mpz_mul_si(tmp, coeff[i], ai); //on multiplie a^i par la valeur du i-eme coefficient
			ai = ai*a;			//on augmente la puissance de a par 1
			mpz_div_2exp(tmp, tmp, k*i);	//on divisse le numerateur par la i-eme puissance de 2^k
			mpz_add(*res, *res, tmp);  	//on fait la somme des deux termes
		}

		mpz_clear(tmp);

}


/* METHODE 1: 

Soit Q un polynome a coefficients entiers notes c0, c1, ... , cn et n son degre:

Q(a/(2^k)) = c0*(a/(2^k))⁰ + c1*(a/(2^k))¹ +...+ cn*(a/(2^k))^n */

/* 

La valeur obtenu comme resultat en utilisant cette fonction est bien precise, elle marche de la facon suivante:
L idee est de faire la somme deux a deux des termes ci-dessus, sans utiliser le fait qu on sait que le denominateur commun sera celui du deuxieme terme.
D abord nous creons les variables temporaires qui stockerons les resultats des sommes et multiplications a chaque iteration, egalement des variables
qui prendront les valeurs des numerateurs et denominateurs.
Donc, nous avons le num1 et den1 qui seront le numerateur et denominateur du terme a gauche. Aussi, nous avons num2 et den2 qui represent le terme a droite.
A chaque iteration nous faisons le produit de num1 et den2, ainsi que celui de den1 et num2. Apres nous faisons la somme de ce produit, qui sera
utilise comme num1 pour la prochaine iteration, et nous donnerons a num2 la valeur qui correspond a son degre, c est-a-dire, le produit de son coefficient correspondant par a^i,
i etant le degre du terme.
Pour le traitement de denominateur, nous faisons le produit de deux denominateurs a chaque iteration, pourvu que nous devons sommer en utilisant un denominateur commun

*/

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
	mpz_set_si(den1, 1);	//le premier terme a comme denominateur 1, qui est (2^k)^0

	mpz_t den2;
	mpz_init(den2);
	mpz_set_si(den2, 1);
	mpz_mul_2exp(den2,den2, k); //le deuxieme terme a comme denominateur (2^k)^1

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
	gmp_printf("numerateur: %Zd\n", *num);
	gmp_printf("denominateur: %Zd\n", *den);

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
		gmp_printf("denominateur: %Zd\n", *den);
		gmp_printf("numerateur: %Zd\n", *num);
		gmp_printf("den1: %Zd\n", den1);
		gmp_printf("num1: %Zd\n", num1);
		gmp_printf("den2: %Zd\n", den2);
		gmp_printf("num2: %Zd\n", num2);
	}
	
	mpz_clear(num1);
	mpz_clear(num2);
	mpz_clear(tmp1);
	mpz_clear(tmp2);
	mpz_clear(den1);
	mpz_clear(den2);
}

/* METHODE 2: 

Soit Q un polynome a coefficients entiers dont la taille peut etre tres grande, notes c0, c1, ... , cn 
et n son degre:
Q(a/(2^k)) = (c0*a⁰*2^(k*n) + c1*a¹*2^(k*(n-1)) +...+ cn*a^n*2^(k*(n-n)))/2^(n*k) */

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

/*	 EVALUATION D APRES LA METHODE DE HORNER

La methode de Horner pour un polynome P(X) = c0 + c1*X + c2*X^2 + ... + cn*X^n  ou n est son degre et donc
le coefficient dominant n est pas nul, permet d evaluer un polynome quelconque de la maniere suivante:
              P(X) = c0 + X(c1 + X(c2 + ... + X(cn-1 + cn*X)))
	      
Appliquee dans notre cas, pour un polynome evalue en X = a/2^k, cela donne:
              Q(a/2^k) = c0 + (a/2^k)(c1 + (a/2^k) * (c2 + ... + (a/2^k) * (cn-1 + cn * (a/2^k))))
	      
Soit, en mettant tout sur le meme denominateur:
	a) au numerateur: c0*2^(nk) + c1*a*2^(n-1)k + ... + cn-3*a*2^(3k) + cn-2*a*2^(2k) + cn-1*a*2^k + cn*a
	b) au denominateur: 2^(nk)
*/

void eval_poly_horner(mpz_t *coeff, int a, unsigned int k, unsigned long int deg, mpz_t *num, mpz_t *den){
    mpz_t tmp1;
    mpz_init(tmp1);
    mpz_t tmp2;
    mpz_init(tmp2);

    mpz_set_si(*den, 1);
    mpz_set_si(*num, 0);
    long int d = deg;
    long int i;
    
    for (i = d; i >= 0; i--){
      	mpz_mul_si(tmp1, *num, a);
	mpz_set_si(*den, 1);
	mpz_mul_2exp(*den, *den, (deg-i)*k);
	mpz_mul(tmp2, coeff[i], *den);
 	mpz_add(*num, tmp2 , tmp1);
	//gmp_printf("tmp1: %Zd\n tmp2: %Zd\n num:%Zd\n den: %Zd\n ", tmp1, tmp2, *num, *den);
  		
	
  }
      mpz_clear(tmp1);
      mpz_clear(tmp2);
}

int main(){

	/* Ecriture/lecture de polynôme */

	unsigned long int deg;
	printf("Génération d'un fichier de degré n/\nn=");
	scanf("%lu", &deg);

	mpz_t *poly = malloc((deg+1)*sizeof(mpz_t));

	random_coeff("coefficients_p.c",deg); //Ecriture de coefficients de type mpz_t aléatoires dans le fichier
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


	/* Methode 1-bis*/
	mpz_t den1bis;
	mpz_t num1bis;
	mpz_init(num1bis);
	mpz_init(den1bis);

	eval_poly_1bis(poly, a, k, deg, &num1bis, &den1bis);
	gmp_printf("Resultat méthode 1-bis: %Zd / %Zd\n", num1bis, den1bis);

	mpz_clear(num1bis);
	mpz_clear(den1bis);

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
