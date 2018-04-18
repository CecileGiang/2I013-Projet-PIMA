//10/04/2018 Weijie YE
//header pour la version impelementer par cecile.
#ifndef POLYNOME_
#define POLYNOME_


/* Fonction de lecture des coefficients d'un polynôme depuis un fichier c (parseur) */
void parse_file(char *nom_fichier, mpz_t *polynome, unsigned long int deg);

/* Fonction d'écriture de coefficients aléatoires d'un polynôme dans un fichier c */
void set_coefficients(char *nom_fichier, unsigned long int deg,unsigned int max);



void eval_poly_1(mpz_t *coeff, long long a, long long k, unsigned long int deg, mpz_t *res_num, mpz_t *res_den);
void eval_poly_2(mpz_t *coeff, long long a, long long k, unsigned long int deg, mpz_t *res_num, mpz_t *res_den);
//void evaluateHorner(mpz_t *coeff, long long a, long long k, unsigned long int deg, mpz_t *res_num, mpz_t *res_den);
void evaluateHorner(mpz_t *coeff, long long a, long long k, unsigned long int deg, mpz_t *res_num, mpz_t *res_den);
//cette fonction efface les donnees qui sont dans le fichier ;
void effacer_fichier(char *nom_fichier);


// void enregistrement(char *nomfic, unsigned int i,float time,float time2);
void enregistrement_1(char *nomfic,unsigned long int n,float time);
void enregistrement_2(char *nomfic,unsigned long int n,float time, float time2);
void enregistrement_3(char *nomfic,unsigned long int n,float time, float time2,float time3);
//Cette fonction ecrit le resultat obtenu dans un fichier .
void ecrit_result(char *nom_fichier, mpz_t res);

#endif
