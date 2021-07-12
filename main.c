//
//  main.c
//  newtonfrac
//
//  Created by Arnau Jutglar Puig on 10/07/2021.
//  Copyright © 2021 Arnau Jutglar Puig. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    double re;
    double im;
} C;

C suma(C, C);
C restar(C, C);
C multi(C, C);
C dividir(C, C);
C escalar(C, double);
C potencia(C, int);
double norma(C);
int base2(int [], int, int);
void polis(int **, int);
C evaluar(int [], int, C);
void derivar(int **, int **, int);
int newton(int *, int *, int, C, C *);
void imprimir(C z);
int omplearrels(int *, int *, C *, int);
void grauN(int n, FILE *);

double tol = 1e-4, h = 1, l1 = 10, l2 = -10;
int iter = 15;


int main(void) {
    int N=10, n;
    FILE *fp;
    
    fp = fopen("arrelsfrac.txt", "w");
    
    for (n=1; n <= N; n++) {
        grauN(n, fp);
    }
    
    return 0;
}


C suma(C z1, C z2) {
    C z;
    
    z.re = z1.re + z2.re;
    z.im = z1.im + z2.im;
    
    return z;
}

C restar(C z1, C z2) {
    
    return suma(z1, escalar(z2, -1));
}

C multi(C z1, C z2) {
    C z;
    
    z.re = z1.re * z2.re - z1.im * z2.im;
    z.im = z1.re * z2.im + z1.im * z2.re;
    
    return z;
}

C dividir(C z1, C z2) {
    C z, z2conj;
    
    z2conj.re = z2.re;
    z2conj.im = -z2.im;
    
    z = multi(z1, z2conj);
    
    z.re /= z2.re * z2.re + z2.im * z2.im;
    z.im /= z2.re * z2.re + z2.im * z2.im;
    
    return z;
}

C escalar(C z, double a) {
    
    z.re *= a;
    z.im *= a;
    
    return z;
}

C potencia(C z, int n) {
    C z0, r;
    int i;
    
    z0.re = 1;
    z0.im = 0;
    
    r.re = 1;
    r.im = 0;
    
    // No ho entenc
    /*
    if (n == 0) {
        return z0;
    } else if (n < 0) {
        printf("Exponent negatiu.\n");
        exit(1);
    } else if (n == 1) {
        printf("%d: ", n);
        imprimir(z);
        return z;
    }
    
    printf("%d: ", n);
    imprimir(z);
    return multi(z, potencia(z, n-1));*/
    
    if (n == 0) {
        return z0;
    }
    
    for (i=0; i < n; i++) {
        r = multi(r, z);
    }
    
    return r;
}

double norma(C z) {
    return sqrt(z.re * z.re + z.im * z.im);
}

int base2(int v[], int n, int y) {
    int i = 0, aux;
        
    for (i=0; i <= n; i++) {
        v[i] = -1;
    }
    
    i=0;
    
    while(y > 0) {
        aux = y % 2;
        
        if (aux == 1) {
            v[i] = 1;
        }
        y /= 2;
        i++;
    }
    
    // Retorna num xifres
    
    return i;
}

void polis(int **v, int n) {
    int i;
    
    /* n és el grau i indica el nombre de xifres màximes. Per tant, calcula 2^(n+1) i aquest és el numbru. */
    
    
    
    for (i=0; i < pow(2, n+1); i++) {
        base2(v[i], n, i);
        // Això omple cada vector amb el polinomi que toca.
    }
    
    return;
}

C evaluar(int v[], int n, C z) {
    int i;
    C r;
    
    r.re = 0;
    r.im = 0;
    
    // Activa això quan ho facis amb coeficients de veritat
    for (i=0; i <= n; i++) {
        r = suma(r, escalar(potencia(z, i), v[i]));
    }
    
    return r;
}

void derivar(int **v, int **dv, int n) {
    // 2^(n+1) vectors de mida n+1
    int i, j;
    
    for (j=0; j < pow(2, n+1); j++) {
        dv[j][n] = 0;
        for (i=n; i > 0; i--) {
            dv[j][i-1] = v[j][i] * i;
        }
    }
    
    return;
}

int newton(int *v, int *dv, int n, C z0, C *z) {
    int i=0;
    C f;
    C df;
    
    while (norma(evaluar(v, n, z0)) > tol && i < iter) {
        df = evaluar(dv, n, z0);
        
        if (norma(df) < tol) {
            // Divisió entre 0
            return 1;
        }
        
        f = evaluar(v, n, z0);
        z0 = restar(z0, dividir(f, df));
        i++;
    }
    
    if (i == iter) {
        *z = z0;
        return 1;
    }
    
    // Convergeix
    *z = z0;
    
    return 0;
}

void imprimir(C z) {
    printf("%lf %lfi\n", z.re, z.im);
    
    return;
}

int omplearrels(int *v, int *dv, C *arrels, int n) {
    int i, j, k, flag, conv, n_arrels=0;;
    C z0;
    C *z;
    
    z = malloc(sizeof(C));
    
    for (i = l2; i <= l1; i += h) {
        z0.re = i;
        for (j = l2; j <= l1; j += h) {
            z0.im = j;
            conv = newton(v, dv, n, z0, z);
            if (conv == 0) {
                // Convergeix
                flag = 0;
                for (k=0; k < n_arrels; k++) {
                    if (norma(restar(*z, arrels[k])) < tol) {
                        flag = 1;
                        break;
                    }
                }
                if (flag == 0 && n_arrels < n) {
                    // Llavors no hi és
                    arrels[n_arrels] = *z;
                    n_arrels++;
                } else if (flag == 0 && n_arrels == n) {
                    return 1;
                }
            }
        }
    }
    
    if (n_arrels == n) {
        return 0;
    } else {
        return -1;
    }
}

void grauN(int n, FILE *fp) {
    int i, j;
    int **v, **dv;
    C **arrels;
    C *z;
    
    z = calloc(1, sizeof(C));

    // Reservem memòria per la matriu de polinomis
    v = calloc(pow(2,n+1), sizeof(int *));
    if (v == NULL) {
        printf("Problemes memoria.\n");
        exit(1);
    }
    
    for (i=0; i < pow(2,n+1); i++) {
        v[i] = calloc(n+1, sizeof(int));
        if (v[i] == NULL) {
            printf("Problemes de memoria.\n");
            exit(1);
        }
    }
    
    dv = calloc(pow(2, n+1), sizeof(int *));
    if (dv == NULL) {
        printf("Problemes memòria.\n");
        exit(1);
    }
    
    for (i=0; i < pow(2, n+1); i++) {
        dv[i] = calloc(n+1, sizeof(int));
        if (dv[i] == NULL) {
            printf("Problemes memòria.\n");
            exit(1);
        }
    }
    
    arrels = calloc(pow(2, n+1), sizeof(C *));
    if (arrels == NULL) {
        printf("Problemes de memoria.\n");
        exit(1);
    }
    
    for (i=0; i < pow(2,n+1); i++) {
        arrels[i] = calloc(n, sizeof(C));
        if (arrels[i] == NULL) {
            printf("Problemes memòria.\n");
            exit(1);
        }
    }
    
    // Omplim la matriu amb tots els polinomis
    polis(v, n);
    
    for (j = 0; j < pow(2,n+1); j++) {
        printf("%d: ", j);
        for (i=0; i <= n; i++) {
            printf("%d ", v[j][i]);
        }
        printf("\n");
    }
    
    // Calculem les derivades
    derivar(v, dv, n);
    
    
    // Calculem les arrels de cada polinomi
    for (j=0; j < pow(2, n+1); j++) {
        omplearrels(v[j], dv[j], arrels[j], n);
    }
    
    // Escrivim al fitxer
    for (j=0; j < pow(2, n+1); j++) {
        for (i=0; i < n; i++) {
            fprintf(fp, "%lf %lf \n", arrels[j][i].re, arrels[j][i].im);
        }
    }
    
    // Alliberem la memòria
    for (i=0; i < pow(2, n+1); i++) {
        free(v[i]);
    }
    free(v);
    
    for (i=0; i < pow(2,n+1); i++) {
        free(dv[i]);
    }
    free(dv);
    
    for (i=0; i < pow(2,n+1); i++) {
        free(arrels[i]);
    }
    free(arrels);
    
    return;
}
