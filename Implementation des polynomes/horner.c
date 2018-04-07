void eval_poly_horner(mpz_t *coeff, int a, unsigned int k, unsigned long int deg, mpz_t *num, int *den){
    mpz_t tmp1;
    mpz_init(tmp1);
    mpz_t tmp2;
    mpz_init(tmp2);

    long int d = deg;
    long int i;
    *den = 1;
    mpz_set(*num, coeff[0]);
    
    for (i = d; i > 0; i--){
      mpz_mul_si(tmp1, *num, a);
      mpz_add(tmp2, coeff[i], tmp1);
      mpz_add(*num, *num, tmp2);
      *den = (*den)*(1<<k); 		// probleme: a partir d un certain rang, le denominateur vaut 0
	//printf("den: %d\n", *den);
  }
      mpz_clear(tmp1);
      mpz_clear(tmp2);
}
