//Modif 10/04/2018 Weijie YE

# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>
# include <string.h>
# include <math.h>
# include <time.h>
# include <gmp.h>

//Cette fontion genere les randnums qui appartient a Z;
void set_coefficients(char *nom_fichier, unsigned long int deg,unsigned int max){

    FILE *f;
    f = fopen(nom_fichier, "w+");

    if(f==NULL){
        fprintf(stderr, "Impossible d'ouvrir le fichier\n");
        exit(1);
    }
    srand(time(NULL));
    mpz_t x;
    mpz_t randnum;

    mpz_init(randnum); //randnum prendra les valeurs des mpz_t générés aléatoirement
    mpz_init(x);
    mpz_set_si(x, -1);

    gmp_randstate_t seed;
    gmp_randinit_default(seed);
    gmp_randseed_ui(seed, rand()); //Paramétrage du noyau pour la génération aléatoire (cf. NOTE)

    unsigned long int i;
    for(i=0; i<=deg; i++){
        mpz_urandomb(randnum, seed, max);

        //pour generer des coeffs negatif avec une probabilite de 50%
        if(rand()%2==0){
            mpz_mul(randnum,randnum,x);
        }

        mpz_out_str (f, 10, randnum); //mpz_out_str (FILE *stream, int base, const mpz_t op) écrit dans le flux entré en paramètre le mpz_t op dans la base indiquée
        fputs("\n",f); //retour à la ligne
    }

    gmp_randclear(seed);
    mpz_clear(randnum);

    fclose(f);
}

void parse_file(char *nom_fichier, mpz_t *polynome, unsigned long int deg){

	FILE *f;
	f = fopen(nom_fichier, "r");

	unsigned long int i = 0;

	if(f==NULL){
		fprintf(stderr, "Impossible d'ouvrir le fichier\n");
		exit(1);
	}

	while(i<=deg){
		mpz_inp_str(polynome[i], f, 10);

		i++;
	}

	fclose(f);
}


void enregistrement_1(char *nomfic,unsigned long int n,float time, float time2){

    FILE *f = fopen(nomfic,"a");
    if (f==NULL) {
        printf("erreur d'ouverture du fichier \n");
        exit(1);
    }
    fprintf(f,"%lu ; %f\n",n,time);

    fclose(f);
}

void enregistrement_2(char *nomfic,unsigned long int n,float time, float time2){

    FILE *f = fopen(nomfic,"a");
    if (f==NULL) {
        printf("erreur d'ouverture du fichier \n");
        exit(1);
    }
    fprintf(f,"%lu ; %f ;%f\n",n,time,time2);

    fclose(f);
}


void enregistrement_3(char *nomfic,unsigned long int n,float time, float time2,float time3){

    FILE *f = fopen(nomfic,"a");
    if (f==NULL) {
        printf("erreur d'ouverture du fichier \n");
        exit(1);
    }
    fprintf(f,"%lu ; %f ; %f ; %f  \n",n,time,time2,time3);

    fclose(f);
}


void effacer_fichier(char *nom_fichier){
    FILE *f;
    f = fopen(nom_fichier,"w+");
    if (f== NULL)
       {
       fprintf(stderr,"Impossible d'ouvrir le fichier\n");
       exit(1);
       }
    //fputs("\n",f);

    fclose(f);
}


void ecrit_result(char *nom_fichier, mpz_t res){

    FILE *f;
    f = fopen(nom_fichier, "a");

    if(f==NULL){
        fprintf(stderr, "Impossible d'ouvrir le fichier\n");
        exit(1);
    }

    mpz_out_str (f, 10, res);
    fputs("\n",f); //deux fois retour à la ligne
    fputs("\n",f);



    fclose(f);
}
