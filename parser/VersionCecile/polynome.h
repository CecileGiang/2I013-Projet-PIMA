#ifndef POLYNOME_
#define POLYNOME_

/* FONCTIONS DE LECTURE/ECRITURE DE POLYNÔMES - PARSER */

/* Fonction de lecture des coefficients d'un polynôme depuis un fichier c (parseur) */
void parse_file(char *nom_fichier, mpz_t *polynome, unsigned long int deg);

/* Fonction d'écriture de coefficients aléatoires d'un polynôme dans un fichier c */
void set_coefficients(char *nom_fichier, unsigned long int deg);

/* NOTE: 
Génération aléatoire de mpz_t. Ces mpz_t seront ensuite stockés dans un fichier qui pourra être lu par la fonction parse_file. Le nombre de mpz_t créés sera de n avec n e degré entré par l'utilisateur. 


Pour initialiser le noyau (seed) permettant de générer aléatoirement des nombres, il nous faut déclarer un paramètre (dit "random state parameter"), et l'initialiser par la fonction gmp_randinit_defaut.
Le noyau initial sera ensuite modifié avec la fonction gmp_randseed_ui.

La taille du noyau détermine le nombre de séquences de nombres aléatoires qu'il est possible de générer.

Une fois le noyau correctement paramétré, l'appel à mpz_rrandomb (mpz_t rop, gmp_randstate_t state, mp_bitcnt_t n) permet de générer un mpz_t aléatoire, compris entre 0 et 2^n-1 inclus, et de le stocker dans rop.
*/





/* FONCTIONS D'EVALUATION DES POLYNÔMES */
/* Par la suite nous ferons l'hypothèse que a et k sont compris entre (-2)⁶³ et 2⁶³-1, nous les stockerons donc dans des long long (int) */

/* METHODE 1: q(a/(2^k)) = c0*(a/(2^k))⁰ + c1*(a/(2^k))¹ +...+ cn*(a/(2^k))^n, où n est le degré du polynôme. */
void eval_poly_1(mpz_t *tab, long long a, long long k, unsigned long int deg, mpz_t res);

/*METHODE 2: q(a/(2^k)) = (c0*a⁰*2^(k*n) + c1*a¹*2^(k*(n-1)) +...+ cn*a^n*2^(k*(n-n)))/2^(n*k), où n est le degré du polynôme. */
void eval_poly_2(mpz_t *coeff, long long a, long long k, unsigned long int deg, mpz_t *res);

#endif
