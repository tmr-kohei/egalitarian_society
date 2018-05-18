#include <stdio.h>
#include <math.h>

#define R 1
#define s_strong 2
#define s_weak 1
#define gamma 1
#define phi 0.2
#define w0 0
#define T 10

long double factorial(int a) {
    if (a > 0) {
        return a * factorial(a - 1);
    } else {
        return 1;
    }
}

// Function to calculate the multinomial distribution
long double multi(int k, int i, int j) {
	long double f1,f2,f3,f4;
	long double result;

	f1 = factorial(k+i+j);
	f2 = factorial(k);
	f3 = factorial(i);
	f4 = factorial(j);
	result = (long double)f1/f2;
	result = (long double)result/f3;
	result = (long double)result/f4;
	return(result);
}

// Function for evolutionary dynamics
int evo_dynamics(long double alpha, long double beta, long double c, int n) {

	FILE *fp;

	long double x00 = 0.97; //E0C0
	long double x01 = 0.01; //E0C1
	long double x10 = 0.01; //E1C0
	long double x11 = 0.01; //E1C1
	
	long double x00_prime, x01_prime, x10_prime, x11_prime;
	long double p, p_strong, p_weak;
	long double x00_strong, x01_strong, x10_strong, x11_strong;
	long double x00_weak, x01_weak, x10_weak, x11_weak;

	long double rho00, rho01, rho10, rho11;
	long double pi00, pi01, pi10, pi11;
	long double pi00_strong, pi01_strong, pi10_strong, pi11_strong;
	long double pi00_weak, pi01_weak, pi10_weak, pi11_weak;
	long double rho00_strong, rho01_strong, rho10_strong, rho11_strong;
	long double rho00_weak, rho01_weak, rho10_weak, rho11_weak;

	long double w00, w01, w10, w11, W;
	long double w00_strong, w01_strong, w10_strong, w11_strong;
	long double w00_weak, w01_weak, w10_weak, w11_weak;

	long double prob_win;
	long double sigma1 = 0.5;
	long double sigma2 = 0.5;
	long double s;
	long double q;

	int i, j, k, t;

	printf("alpha,beta,c,n,x00,x01,x10,x11\n");

	for (t=0; t<T; t++) {

		x00_strong = x00*phi;
		x00_weak = x00*(1-phi);
		x01_strong = x01*phi;
		x01_weak = x01*(1-phi);
		x10_strong = x10*phi;
		x10_weak = x10*(1-phi);
		x11_strong = x11*phi;
		x11_weak = x11*(1-phi);

		p_strong = x01_strong + x11_strong;
		p_weak = x01_weak + x11_weak;
		p = x01 + x11;

		// Payoff as an owner
		pi00_strong = 0.0;
		for (i=0; i<=(n-1); i++) {
			for (j=0; j<=(n-1-i); j++) {
				k = i + j;
				if (k == 0) {
					q = multi(n-1-k,i,j)*pow((1-p),(n-1-k))*pow(p_strong,i)*pow(p_weak,j);
					pi00_strong = pi00_strong + q*R;
				} else {
					q = multi(n-1-k,i,j)*pow((1-p),(n-1-k))*pow(p_strong,i)*pow(p_weak,j);
					s = (long double)(s_strong*i + s_weak*j)/k;
					prob_win = 1/( 1+exp(-gamma*(s_strong-pow(k,alpha)*s )) );
					pi00_strong = pi00_strong + q*( prob_win*R + (1-prob_win)*(-c) );
				}
			}			
		}

		pi00_weak = 0.0;
		for (i=0; i<=(n-1); i++) {
			for (j=0; j<=(n-1-i); j++) {
				k = i + j;
				if (k == 0) {
					q = multi(n-1-k,i,j)*pow((1-p),(n-1-k))*pow(p_strong,i)*pow(p_weak,j);
					pi00_weak = pi00_weak + q*R;
				} else {
					q = multi(n-1-k,i,j)*pow((1-p),(n-1-k))*pow(p_strong,i)*pow(p_weak,j);
					s = (long double)(s_strong*i + s_weak*j)/k;
					prob_win = 1/(1+exp(-gamma*(s_weak-pow(k,alpha)*s )));
					pi00_weak = pi00_weak + q*( prob_win*R + (1-prob_win)*(-c) );
				}
			}			
		}

		pi01_strong = pi00_strong;
		pi01_weak = pi00_weak;

		pi10 = (long double)R/n;
		pi11 = (long double)R/n;

		//Payoff as a peer
		rho00 = (x10+x11)*R/n + (x00+x01)*(0);

		rho01_strong = (x10+x11)*R/n;
		for (i=0; i<=(n-2); i++) {
			for (j=0; j<=(n-2-i); j++) {
				k = i + j + 1;
				q = multi(n-1-k,i,j)*pow((1-p),(n-1-k))*pow(p_strong,i)*pow(p_weak,j);
				s = (long double)(s_strong*(i+1) + s_weak*j)/k;
				prob_win = phi/(1+exp(gamma*(s_strong-pow(k,alpha)*s ))) + (1-phi)/(1+exp(gamma*(s_weak-pow(k,alpha)*s )));
				rho01_strong = rho01_strong + q*(x00+x01)*( prob_win*R/k + (1-prob_win)*(-c) );
			}			
		}

		rho01_weak = (x10+x11)*R/n;
		for (i=0; i<=(n-2); i++) {
			for (j=0; j<=(n-2-i); j++) {
				k = i + j + 1;
				q = multi(n-1-k,i,j)*pow((1-p),(n-1-k))*pow(p_strong,i)*pow(p_weak,j);
				s = (long double)(s_strong*i + s_weak*(j+1))/k;
				prob_win = phi/(1+exp(gamma*(s_strong-pow(k,alpha)*s ))) + (1-phi)/(1+exp(gamma*(s_weak-pow(k,alpha)*s )));
				rho01_weak = rho01_weak + q*(x00+x01)*( prob_win*R/k + (1-prob_win)*(-c) );
			}			
		}
		
		rho10 = (x10+x11)*R/n + (x00+x01)*(0);

		rho11_strong = (x10+x11)*R/n;
		for (i=0; i<=(n-2); i++) {
			for (j=0; j<=(n-2-i); j++) {
				k = i + j + 1;
				q = multi(n-1-k,i,j)*pow((1-p),(n-1-k))*pow(p_strong,i)*pow(p_weak,j);
				s = (long double)(s_strong*(i+1) + s_weak*j)/k;
				prob_win = phi/(1+exp(gamma*(s_strong-pow(k,alpha)*s ))) + (1-phi)/(1+exp(gamma*(s_weak-pow(k,alpha)*s )));
				rho11_strong = rho11_strong + q*(x00+x01)*( prob_win*R/k + (1-prob_win)*(-c) );
			}			
		}

		rho11_weak = (x10+x11)*R/n;
		for (i=0; i<=(n-2); i++) {
			for (j=0; j<=(n-2-i); j++) {
				k = i + j + 1;
				q = multi(n-1-k,i,j)*pow((1-p),(n-1-k))*pow(p_strong,i)*pow(p_weak,j);
				s = (long double)(s_strong*i + s_weak*(j+1))/k;
				prob_win = phi/(1+exp(gamma*(s_strong-pow(k,alpha)*s ))) + (1-phi)/(1+exp(gamma*(s_weak-pow(k,alpha)*s )));
				rho11_weak = rho11_weak + q*(x00+x01)*( prob_win*R/k + (1-prob_win)*(-c) );
			}			
		}

		// Fitness of each phenogenotypes
		w00_strong = w0 + sigma1*pi00_strong + sigma2*rho00;
		w00_weak = w0 + sigma1*pi00_weak + sigma2*rho00;
		w01_strong = w0 + sigma1*pi01_strong + sigma2*rho01_strong;
		w01_weak = w0 + sigma1*pi01_weak + sigma2*rho01_weak;
		w10_strong = w0 + sigma1*pi10 + sigma2*rho10;
		w10_weak = w0 + sigma1*pi10 + sigma2*rho10;
		w11_strong = w0 + sigma1*pi11 + sigma2*rho11_strong;
		w11_weak = w0 + sigma1*pi11 + sigma2*rho11_weak;

		// Mean fitness
		W = exp(beta*w00_strong)*x00_strong + exp(beta*w00_weak)*x00_weak + exp(beta*w01_strong)*x01_strong + exp(beta*w01_weak)*x01_weak + exp(beta*w10_strong)*x10_strong + exp(beta*w10_weak)*x10_weak + exp(beta*w11_strong)*x11_strong + exp(beta*w11_weak)*x11_weak;

		// The frequencies of four phenogenotypes in the next generation
		x00_prime = exp(beta*w00_strong)*x00_strong/W + exp(beta*w00_weak)*x00_weak/W;
		x01_prime = exp(beta*w01_strong)*x01_strong/W + exp(beta*w01_weak)*x01_weak/W;
		x10_prime = exp(beta*w10_strong)*x10_strong/W + exp(beta*w10_weak)*x10_weak/W;
		x11_prime = exp(beta*w11_strong)*x11_strong/W + exp(beta*w11_weak)*x11_weak/W; 

		x00 = x00_prime;
		x01 = x01_prime;
		x10 = x10_prime;
		x11 = x11_prime;

		printf("%Lf,%Lf,%Lf,%d,%Lf,%Lf,%Lf,%Lf\n", alpha, beta, c, n, x00, x01, x10, x11);
	}

	return(0);
}

int main(void) {
	long double alpha, beta, c;
	int n;
	long double a;

	alpha = 2;
	n = 20;
	beta = 1;
	c = 1.0;

	a = evo_dynamics(alpha, beta, c, n);
	
	return(0);
}
