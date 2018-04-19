# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>
# include <string.h>
# include <math.h>
# include <time.h>
# include <gmp.h>

/* _______________FONCTIONS DE LECTURE/ECRITURE DE POLYNOMES DANS UN FICHIER C___________ */

/* 

Cette fonction nous permet de lire les coefficients de type mpz_t ecrits dans un fichier sur la forme d une liste. 
Nous lisons donc chaque ligne, ensuite nous stockons le coefficient dans un tableau qui representera notre polynome et puis nous faisons un retour a la ligne. 
Nous faisons cela jusqu a avoir n+1 coefficients, n etant le degre du polynome et qui a ete saisi precedemment.
L indice du tableau polynome correspondra bien au degre du terme, ainsi le premier terme (d indice 0) sera le terme constant, le deuxieme terme (d indice 1) celui de degre 1, et ainsi de suite. 

*/

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
Ce noyau donne plus tard des valeurs aleatoires (entre 0 et 2^(n-1), n etant le troisieme argument de la fonction qui represente la taille binaire souhaitee)
a chaque coefficient.
Comme nous souhaitons aussi avoir des nombres negatifs, nous les generons aleatoirement, c est-a-dire, nous creons une variable signe,
qui determinera si un nombre est positif ou negatif avec une equiprobabilite.

*/

void random_coeff(char *nom_fichier, unsigned long int deg, int taille_bin) {
	
	FILE *f = fopen(nom_fichier, "w+");

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
		mpz_rrandomb(coeff, state, taille_bin); //coeff prend une valeur random entre 0 et 2^(n-1) grace a l algorithme pris par state
		if (rand()%2 == 0){ //nous creons artificiellement des entiers negatifs aussi, avec une chance sur deux
			sign = 1;
		}
		else {
			sign = -1;
		}
		mpz_mul_si(coeff, coeff, sign); 
		mpz_out_str(f, 10, coeff); //nous imprimons le nombre en base 10 dans le fichier
		fputs("\n", f); //retour a la ligne
	}
	gmp_randclear(state);
	mpz_clear(coeff);
	fclose(f);


/*      Commentaires: nous devrions avoir des nombres aleatoires entre 0 et 2^(n-1), 
	pourtant nous trouvons que des nombres avec le meme nombre de chiffres a chaque fois  */
}

/* ______FONCTION AUXILIAIRE: REDUCTION D UNE FRACTION_________ */

/* 

Nous utiliserons cette fonction pour avoir la fraction qui donne comme resultat l evaluation du polynome, sur sa forme reduite.
L idee est simple: nous calculons d abord le plus grand commun diviseur du numerateur et denominateur, puis nous divisons chacun par ce valeur et donc, nous aurons ainsi une primalite relative entre le numerateur et le denominateur.

*/

void reduire_fraction(mpz_t *den, mpz_t *num){

	mpz_t gcd;
	mpz_init(gcd);
	mpz_gcd(gcd, *num, *den);

	mpz_divexact(*num, *num, gcd);
	mpz_divexact(*den, *den, gcd);

	mpz_clear(gcd);
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
		if(deg==0) return;

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
L idee est de faire la somme deux a deux des termes ci-dessus, en utilisant le fait qu on sait que le denominateur commun sera celui du deuxieme terme.
Nous initialisons donc le numerateur et le denominateur comme ca se ferait pour un polynome de degre 0, donc le numerateur vaut le coefficient constant, stocke dans le tableau de coefficients sur la case 0, et le denominateur vaut 1.
Nous faisons le test du degre, pour ne plus modifier si c est le cas et le degre vaut bien 0.
Sinon, nous continuos a modifier la valeur du numerateur en la multipliant d abord par 2^k, qui nous permettra d avoir un denominateur commun entre le deux termes, et ainsi faire la somme.
A ce valeur, nous rajouterons le produit entre le coefficient du terme suivant et son a^i correspondant au degre du coefficient. Cette variable ai est modifie dernierement pour que son exposant corresponde au degre.

Apres etre sortis de la boucle, nous mettons a jour la valeur du denominateur en lui donnant 2^(k*deg), car ca sera la valeur necessaire pour mettre en denominateur commun tous les termes deja sommes et dont le resultat de la somme sera dans la varible num.

*/

void eval_poly_1bis(mpz_t *coeff, int a, unsigned int k, unsigned long int deg, mpz_t *num, mpz_t *den){
	
	mpz_set(*num, coeff[0]); //initialisation du numerateur a c0
	mpz_set_si(*den, 1);	//le premier terme a comme denominateur 1, qui est (2^k)^0
	
	if (deg == 0) return;

	unsigned long int i;

	mpz_t ai;
	mpz_init(ai);
	mpz_set_si(ai, a);


	for (i = 1; i <= deg; i++){
		
		//*num = (*num)(2^k) + (ci)(a^i), avec *num = c0 au premier tour
		mpz_mul_2exp(*num, *num, k);
		mpz_addmul(*num, coeff[i], ai);
		mpz_mul_si(ai, ai, a);

	}
	
	//La valeur du denominateur sera toujours 2^(k*deg)
	mpz_mul_2exp(*den, *den, k*deg);

	mpz_clear(ai);
	reduire_fraction(den, num);
}

/* METHODE 2: 
Soit Q un polynome a coefficients entiers dont la taille peut etre tres grande, notes c0, c1, ... , cn et n son degre:
Q(a/(2^k)) = (c0*a⁰*2^(k*n) + c1*a¹*2^(k*(n-1)) +...+ cn*a^n*2^(k*(n-n)))/2^(n*k) */

void eval_poly_2(mpz_t *coeff, int a, unsigned int k, unsigned long int deg, mpz_t *num, mpz_t *den){

	mpz_t ai;
	mpz_t ai_decale;//ai_decale prendra la valeur de ai*2^ki
	mpz_init(ai_decale);
	mpz_init(ai);
	mpz_set_si(ai, a); //ai est la valeur a^i pour le coefficient ci

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
			mpz_mul_2exp(ai_decale, ai, k*(deg-i));
			mpz_mul(tmp, coeff[i], ai_decale);
			mpz_add(*num, *num, tmp);
			mpz_mul_si(ai, ai, a);
		}

		mpz_mul_2exp(*den, *den, k*deg);
	}
	
	mpz_clear(tmp);
	mpz_clear(ai_decale);
	mpz_clear(ai);
	reduire_fraction(den, num);

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

	//Creation des variables temporaires pour stocker le resultats des operations intermediaires pour calculer le denominateur et le numerateur
    mpz_t tmp1;
    mpz_init(tmp1);
    mpz_t tmp2;
    mpz_init(tmp2);

    mpz_set_si(*den, 1);
    mpz_set_si(*num, 0);

	//Nous avons besoin des entiers signes pour faire la distinction entre un nombre positif et un nombre negatif lors de la condition "i >= 0" dans la boucle for
    long int d = deg;
    long int i;
    
    for (i = d; i >= 0; i--){
      	mpz_mul_si(tmp1, *num, a);
	mpz_set_si(*den, 1);
	mpz_mul_2exp(*den, *den, (d-i)*k);
	mpz_mul(tmp2, coeff[i], *den);
 	mpz_add(*num, tmp2 , tmp1);
		
	
  }
      mpz_clear(tmp1);
      mpz_clear(tmp2);
	reduire_fraction(den, num);
}


/* RESULTATS EXPERIMENTAUX EN FAISANT VARIER LE DEGRE */
/*

Cette fonction ecrira dans 6 fichiers de type texte (2 pour chaque methode: 1 pour le numerateur et 1 pour le denominateur) ou nous pourrons voir le denominateur et le numerateur en fonction du degre.
L affichage se fera en deux colonnes, la premiere qui represente le degre et la deuxieme qui donnera la valeur du numerateur ou denominateur.
Le degre maximum qu on veut atteindre sera saisi par l utilisateur, ainsi que la partition.
 
Nous pourrons utiliser ses fichiers ulterieurement pour la creation des graphes a l aide de gnuplot.

D autre part, ce fonction calcule aussi le temps d execution de chaque methode et l affiche directement sur le terminal. Les resultats sont donnes en millisecondes.

 
*/


void imprimer_resultats_deg(char *fichier1num, char *fichier1den, char *fichier2num, char *fichier2den, char *fichier_hornernum, char *fichier_hornerden, mpz_t *poly, int a, unsigned int k, unsigned long int deg_max, int partition, mpz_t *num1bis, mpz_t *den1bis, mpz_t *num, mpz_t *den, mpz_t *num3, mpz_t *den3){ 
	unsigned long int i;
	
			/*Resultats en utilisant la methode 1 */
	FILE *fd = fopen(fichier1num, "w+");
	
	for (i = 0; i <= deg_max; i+=partition){
		clock_t temps_initial = clock();
		eval_poly_1bis(poly, a, k, i, num1bis, den1bis);
		clock_t temps_final = clock();
		double temps_cpu = (double)(temps_final - temps_initial)/CLOCKS_PER_SEC*1000;
		printf("Temps d'execution pour la methode 1, au degre %lu: %1.18e\n", i, temps_cpu);
		fprintf(fd, "%ld\t", i);
		mpz_out_str(fd, 10, *num1bis); //nous imprimons le numerateur en base 10 dans le fichier
		fputs("\n", fd);
	}
	fclose(fd);

	fd = fopen(fichier1den, "w+");	
	for (i = 0; i <= deg_max; i+=partition){
		eval_poly_2(poly, a, k, i, num, den);
		fprintf(fd, "%ld\t", i);
		mpz_out_str(fd, 10, *den); //nous imprimons le nombre dans le fichier
		fputs("\n", fd);
	}
	fclose(fd);	

			/*Resultats en utilisant la methode 2*/
	fd = fopen(fichier2num, "w+");	
	for (i = 0; i <= deg_max; i+=partition){
		clock_t temps_initial = clock();
		eval_poly_2(poly, a, k, i, num, den);
		clock_t temps_final = clock();
		double temps_cpu = (double)(temps_final - temps_initial)/CLOCKS_PER_SEC*1000;
		printf("Temps d'execution pour la methode 2, au degre %lu: %1.18e\n", i, temps_cpu);
		fprintf(fd, "%ld\t", i);
		mpz_out_str(fd, 10, *num); //nous imprimons le nombre dans le fichier
		fputs("\n", fd);
	}
	fclose(fd);	

	fd = fopen(fichier2den, "w+");	
	for (i = 0; i <= deg_max; i+=partition){
		eval_poly_2(poly, a, k, i, num, den);
		fprintf(fd, "%ld\t", i);
		mpz_out_str(fd, 10, *den); //nous imprimons le nombre dans le fichier
		fputs("\n", fd);
	}
	fclose(fd);	


			/* Resultats en utilisant la methode de Horner */

	fd = fopen(fichier_hornernum, "w+");	
	for (i = 0; i <= deg_max; i+=partition){
		clock_t temps_initial = clock();
		eval_poly_2(poly, a, k, i, num, den);
		clock_t temps_final = clock();
		double temps_cpu = (double)(temps_final - temps_initial)/CLOCKS_PER_SEC*1000;
		printf("Temps d'execution pour la methode de Horner, au degre %lu: %1.18e\n", i, temps_cpu);
		fprintf(fd, "%ld\t", i);
		mpz_out_str(fd, 10, *num3); //nous imprimons le nombre dans le fichier
		fputs("\n", fd);
	}
	fclose(fd);
	
	fd = fopen(fichier_hornerden, "w+");	
	for (i = 0; i <= deg_max; i+=partition){
		eval_poly_2(poly, a, k, i, num, den);
		fprintf(fd, "%ld\t", i);
		mpz_out_str(fd, 10, *den3); //nous imprimons le nombre dans le fichier
		fputs("\n", fd);
	}
	fclose(fd);	
}

int main(){

	/* Ecriture/lecture de polynôme */

	unsigned long int deg;
	printf("Veuillez saisir le degre maximum \nn=");
	scanf("%lu", &deg);
	
	int partition;
	printf("Veuillez saisir le pas pour la partition\np=");
	scanf("%d", &partition);

	int taille_bin;
	printf("Veuillez entrer la taille binaire des coefficients\nt=");
	scanf("%d", &taille_bin);

	mpz_t *poly = malloc((deg+1)*sizeof(mpz_t));

	random_coeff("coefficients_p.c",deg, taille_bin); //Ecriture de coefficients de type mpz_t aléatoires dans le fichier
	parse_file("coefficients_p.c",poly, deg); //Ecriture des coefficients dans le tableau tab


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

	/* Methode 1-bis*/
	
	mpz_t den1bis;
	mpz_t num1bis;
	mpz_init(num1bis);
	mpz_init(den1bis);
	
	
	/* Methode 2 */
	mpz_t den;
	mpz_t num;
	mpz_init(num);
	mpz_init(den);


	/* Methode de Horner */
	mpz_t den3;
	mpz_t num3;
	mpz_init(num3);
	mpz_init(den3);

	/* Nous imprimons des resultats pour un degre qui varie entre 0 et le degre maximum saisi par l utilisateur, ainsi que le pas pour la partition */
	imprimer_resultats_deg("meth1num.txt", "meth1den.txt", "meth2num.txt", "meth2den.txt", "meth_hornernum.txt", "meth_hornerden.txt", poly, a, k, deg, partition, &num1bis, &den1bis, &num, &den, &num3, &den3);
	
	/* Liberation d espace memoire */
	mpz_clear(num1bis);
	mpz_clear(den1bis);
	mpz_clear(num);
	mpz_clear(den);
	mpz_clear(num3);
	mpz_clear(den3);
	
	free(poly);

	return 0;
}
