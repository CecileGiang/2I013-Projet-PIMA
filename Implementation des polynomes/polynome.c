# include <stdio.h>
# include <gmp.h>
# include <stdint.h>
# include <stdlib.h>

/*Initialisation des coefficients du polynôme*/

void tab_polynome(mpz_t *tab, int n){ //n est le degré du polynôme. Donc la dernière case est d'indice n.

	int i;
	mpz_t k;
	mpz_init(k);
	printf("Saisie des coefficients\n");
	for(i = 0; i<=n; i++){
		gmp_scanf("%Zd",&k);
		mpz_set(tab[i], k);		
	}
	
	printf("Fin de la saisie.\n");
}

/*METHODE 1: q(a/(2^k)) = c0*(a/(2^k))⁰ + c1*(a/(2^k))¹ +...+ cn*(a/(2^k))^n, où n est le degré du polynôme. */

void eval_poly_1(mpz_t* coeffs, int n, mpz_t a, mpz_t k, mpz_t *res){
	

	mpz_t tmp1, tmp2, ak, m; //le k sert en fait pour la multiplication de la puissance de 2 dans la boucle for ("même valeur" que i)
	mpz_init(tmp1); mpz_init(ak), mpz_init(m); mpz_init(tmp2);
	gmp_printf("m = %Zd\n", m);
	int i;

	for(i = 0; i<=n; i++){
		printf("tour %d\n",i);
		gmp_printf("res=%Zd\n",*res);
		/* Question: division par 2^k: mpz_tdiv_q_2exp (mpz_t q, const mpz_t n, mp_bitcnt_t b) ou
		             mpz_tdiv_r_2exp (mpz_t r, const mpz_t n, mp_bitcnt_t b) ? */
		//Calcul de k*m
		mpz_mul(tmp2, k, m);
		gmp_printf("k = %Zd, m = %Zd, k*n = %Zd\n", k, m, tmp2);
		//tpm2 est la puissance de 2. C'est un mpz_t, on a besoin de le convertir en entier pour la suite.
		printf("haha\n");
		printf("%lu\n",mpz_get_ui(m));
		//Calcul de a^m
		mpz_pow_ui(ak, a, mpz_get_ui(m));
		//gmp_printf("a = %Zd, mpz_get_ui(m) = %lu, a^m = %Zd\n", a, mpz_get_ui(m), ak);
		gmp_printf("a^m=%Zd\n",ak);
		//Calcul de ci*a^m
		mpz_mul(ak,ak,coeffs[i]);
		gmp_printf("ci*(a^m)=%Zd\n",ak);
		//Calcul de (a^m)/2^k*m
		gmp_printf("tmp2=%lu\n",mpz_get_ui(tmp2));
		mpz_tdiv_q_2exp(tmp1, ak, mpz_get_ui(tmp2));
	/*	gmp_printf("a^m=%Zd\n",ak);
		gmp_printf("tmp2(int)=%lu\n",mpz_get_ui(tmp2));
		gmp_printf("tmp1=(a^m)/2^tmp2=%Zd\n",tmp1);
*/		gmp_printf("ci*(a^m)/2^k*m=%Zd\n",tmp1);
		
		//Calcul de ci*((a^m)/2^k*m)
	//	mpz_addmul(*res, coeffs[i], tmp1);
		gmp_printf("res=%Zd\n",*res);
		gmp_printf("tmp1=%Zd\n",tmp1);
		mpz_add(*res, *res, tmp1);
		gmp_printf("res+tmp1=%Zd\t",*res);
		mpz_add_ui(m,m,1);
		
	}
	gmp_printf("res fin de fonction: %Zd\n", *res);

	mpz_clears(tmp1, tmp2,ak, m);
}


/*METHODE 2: q(a/(2^k)) = (c0*a⁰*2^(k*n) + c1*a¹*2^(k*(n-1)) +...+ cn*a^n*2^(k*(n-n)))/2^(n*k), où n est le degré du polynôme. */

void eval_poly_2(mpz_t* coeffs, int n, mpz_t a, mpz_t k, mpz_t *res){


	mpz_t tmp1, tmp2, m, puiss2, ak; //le k sert en fait pour la multiplication de la puissance de 2 dans la boucle for ("même valeur" que i)
	mpz_init(tmp1); mpz_init(tmp2); mpz_init_set_ui(m,(unsigned long int)n); mpz_init(ak); mpz_init(puiss2);

	int i;

	for(i = 0; i<=n; i++){
		mpz_mul(puiss2, k, m);
		mpz_pow_ui(ak, a,i);
		mpz_mul(tmp1, coeffs[i], ak);
		// On convertit la puissance de 2 en entier sous peine d'erreur
		mpz_mul_2exp(tmp2,tmp1,(unsigned long int)puiss2);
		mpz_add(*res,*res,tmp2);
		mpz_set(m,m-1);
	}

	mpz_t puiss; mpz_init(puiss);
	mpz_mul(puiss, k, m);
	mpz_tdiv_q_2exp(*res, *res, (int)tmp2);
	mpz_clears(tmp1, tmp2, m, ak, puiss2, puiss);
}



int main(){
	mpz_t a, k, res;
	mpz_init(a); mpz_init(k); mpz_init(res);
	int n, i;
	printf("Degré du tableau:\n");
	scanf("%d",&n);

	// Création du tableau de coeffs
	mpz_t *tab = malloc((n+1)*sizeof(mpz_t));
	//Initialisation du tableau, sous peine d'erreur
	for(i = 0; i<=n; i++) mpz_init(tab[i]);
	tab_polynome(tab, n);

	printf("Calcul de q(a/(2^k), veuillez rentrer les valeurs de:\na:\n");
	gmp_scanf("%Zd", &a);
	printf("\nk:\n");
	gmp_scanf("%Zd", &k);

	// Test fonction d'évaluation de polynôme 1
	eval_poly_1(tab, n, a, k, &res);
	gmp_printf("Le résultat est: %Zd\n", res);
	printf("%d\n",*res);
/*
	// Test fonction d'évaluation de polynôme 2
	eval_poly_2(tab, n, a, k, &res);
	gmp_printf("Le résultat est: %Zd\n", res); */

	free(tab);
	mpz_clears(a,k,res);

	return 0;
}
