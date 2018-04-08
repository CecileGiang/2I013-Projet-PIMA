
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

void eval_poly_1bis(mpz_t *coeff, int a, unsigned int k, unsigned long int deg, mpz_t *num, int *den){
	mpz_t tmp1;
	mpz_t tmp2;
	mpz_t num1;
	mpz_t num2;
	mpz_init(tmp1);
	mpz_init(tmp2);
	mpz_init(num1);
	mpz_init(num2);

	unsigned long int i;

	int ai = 1;
	int ai_decale = a;

	int den1 = 1;
	int den2 = 1<<k;
	mpz_mul_si(tmp1, coeff[0], ai);
	mpz_mul_si(tmp2, coeff[1], ai_decale);
	mpz_set(num1, tmp1);
	mpz_set(num2, tmp2);

	if (deg == 0){
		mpz_set(*num, num1);
		*den = den1;
		return;
	}

	mpz_mul_si(tmp1, num1, den2);
	mpz_mul_si(tmp2, num2, den1);
	mpz_add(tmp1, tmp1, tmp2);
	mpz_set(*num, tmp1);
	*den = den1*den2;

	for (i = 2; i <= deg; i++){
		mpz_set(num1, *num);
		den1 = *den;
		ai_decale = ai_decale * a;
		mpz_mul_si(tmp1, coeff[i], ai_decale);
		mpz_set(num2, tmp1);
		den2 = den1 * (1<<k);
		mpz_mul_si(tmp1, num1, den2);
		mpz_mul_si(tmp2, num2, den1);
		mpz_add(tmp1, tmp1, tmp2);
		mpz_set(*num, tmp1);
		*den = den1*den2;
	}
	mpz_clear(num1);
	mpz_clear(num2);
	mpz_clear(tmp1);
	mpz_clear(tmp2);
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
	for(i=0; i<=deg; i++) {
		gmp_printf("%Zd\n", poly[i]);
	}
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

	mpz_clear(res1);

	/* Méthode 1 avec des floats */
	//mpf_t res2;
	//mpf_init(res2);
	//eval_poly_1f(poly, a, k, deg, &res2);
	//gmp_printf("Resultat méthode 1 avec des floats: %Zd\n", res2);

	/* Méthode 1-bis*/
	int den1;
	mpz_t num1;
	mpz_init(num1);

	eval_poly_1bis(poly, a, k, deg, &num1, &den1);
	gmp_printf("Resultat méthode 1-bis: %Zd / %d\n", num1, den1);

	mpz_clear(num1);

	/* Méthode 2 */
	int den2;
	mpz_t num2;
	mpz_init(num2);

	eval_poly_2(poly, a, k, deg, &num2, &den2);
	gmp_printf("Resultat méthode 2: %Zd / %d\n", num2, den2);

	mpz_clear(num2);

	free(poly);
	return 0;
}
	
/* PROBLÈMES RENCONTRÉS LORS DE L'ÉVALUATION DE POLYNÔMES:
	PROBLÈME COMMUN AUX 2 FONCTIONS:
	La version flotant de la premiere methode possede un probleme - segfault apparement - 
*/
