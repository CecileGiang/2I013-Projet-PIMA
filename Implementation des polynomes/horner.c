void eval_poly_horner(mpz_t *coeff, int a, unsigned int k, unsigned long int deg, mpz_t *num, int *den){
    mpz_t tmp1;
    mpz_init(tmp1);
    mpz_t tmp2;
    mpz_init(tmp2);

    long int d = deg;
    long int i;
    *den = 1;
   // int den1 = (1<<k);

    for (i = d; i >= 0; i--){
	 if(i!=d) *den = *den*pow(2,k); // a partir d un certain rang, le denominateur vaut 0... pour quoi?
      mpz_mul_si(tmp1, *num, a);
	mpz_mul_si(tmp2, coeff[i], *den);
      mpz_add(*num, tmp2 , tmp1);
	//gmp_printf("tmp1: %Zd\n tmp3: %Zd\n tmp2: %Zd\n num:%Zd\n", tmp1, tmp3, tmp2, *num);
  		
	
  }
      mpz_clear(tmp1);
      mpz_clear(tmp2);
}
