#include<bits/stdc++.h>
// #include <boost/math/common_factor.hpp>
using namespace std;

typedef long long ll;
#define int long long
#define vi vector<int>
#define pb push_back
#define mp make_pair
#define ff first
#define ss second
#define inar(arr,n) for(int i=0; i<n; i++) cin>>arr[i];
#define suma(arr,n) for(int i=0; i<n; i++) sum+=arr[i];
#define for0(i,n) for(int i=0; i<n; i++)
#define for1(i,n) for(int i=1; i<=n; i++)
#define forr(i,n) for(int i=n-1; i>=0; i--)
#define rof(i,a,b) for(int i=a;i<=b;i++)
#define all(v) v.begin(), v.end()
#define rall(v) v.rbegin(), v.rend()
#define fast ios_base::sync_with_stdio(false); cin.tie(0); cout.tie(0)
const int mod = 1e9 + 7;
ll gcd (ll a, ll b) {return ( a ? gcd(b % a, a) : b );}
ll power(ll a, ll n) {ll p = 1; while (n > 0) {if (n % 2) {p = p * a;} n >>= 1; a *= a;} return p;}
ll power(ll a, ll n, ll mod) {ll p = 1; while (n > 0) {if (n % 2) {p = p * a; p %= mod;} n >>= 1; a *= a; a %= mod;} return p % mod;}

vi primes;
bool prime[100001];

// CRYPTOSYSTEMS
void RSA();
void Knapsack_Merkle_Hellman();
void ELGAMEL();
void Miller_Rabin();
void Mersenne_number();
void Fermats();


// NUMBER THEORY
void primes_sieve();
int extended_euclid_steps(int a, int b, int& x, int& y);
void GCD();
void GCD_steps(int a, int b);
void GCD_euler_steps();
void extended_euclid();
int phi(int n);
void PHI();
vector<int> prime_factors(int n);
void prime_factors();
void primitive_root();
void POWER_steps();
void ORDER_amodm();
int inv(int a, int m);
bool iscong(int x, int b, int m);
void findMinX(int num[], int rem[], int k);
void Chinese_rem();
void Chinese_rem(int n, int rem[], int num[]);
void Chinese_converted();

void INVERSE_euler();
int order(int a, int m);
void polypower();
void polypower_pow();
void PRIMITIVE_ROOT_steps();
int smallest_p_root(int m);
void INDEX_of_primitive_root();
void Linear_cong();
void Jacobi_symbol(int a, int q);
void Legendre_symbol();
int Legendre_symbol(int a, int n);

void Lucas_primality_test();
void Solovay_Strassen();

signed main() {
	fast;
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
#endif

	//////////////////////////////////////////////////////////////////////////////////////
	/* NUMBER THEORY  */

	// extended_euclid();		//a,b
	// GCD();					//t,a,b
	// GCD_euler_steps();		//a,b
	// INVERSE_euler();			//a,m
	// PHI();					//t,n
	// primes_sieve();
	prime_factors();			//n
	// primitive_root();		//n
	// POWER_steps();			//a,b,m 	a^b mod m
	// ORDER_amodm();			//a,m
	// polypower();				//y,b,m 	x^y cong bmodm
	// polypower_pow();				//y,b,m 	y^x cong bmodm
	// Fermats();	// a,n

	// int c = 209;
	// int d = (c - 1) / 2;
	// int n = 357;
	// cout << 31951492011 % 33 << "\n";
	// cout << power(21, 7, 33) << "\n";

	// cout << order(1, 7) << "\n";
	// cout << (power(157, phi(2668) - 1, 2668));
	// cout << phi(9);
	// cout << "\n\nbleah\n" << gcd(145 - 197, 221);
	// cout << iscong(17, 3, 7) << "\n";
	// cout << power((2396 * inv(power(435, 765, 2579), 2579)), 1, 2579) << "\n\n";
	// float f = 115 * ((1 - 1.0 / 5.0) * (1 - 1.0 / 23.0));
	// cout << f << phi(115);
	// cout << endl << 23455 % 751 << " " << 12121 % 9721 << "  " << 10346 % 12121;

	// Linear_cong();	//t, i=1 to t: {a, b, m}	ax cong bmodm
	// Chinese_rem();	//n, i=1 to n: {b[i], m[i]}
	// Chinese_converted();	//n, i=1 to n: {a[i], b[i], m[i]}

	//////////////////////////////////////////////////////////////////////////////////////
	/* CRYPTOSYSTEMS */

	// RSA();			//p,q, e, message(plaintext)
	// Knapsack_Merkle_Hellman();	//n,W,q,r, message m
	// ELGAMEL();					//A, a, p, k, messsage

	//////////////////////////////////////////////////////////////////////////////////////
	/* NUMBER THEORY REVISITED */
	// PRIMITIVE_ROOT_steps();	//m
	// INDEX_of_primitive_root();	//m

	/* AFTER MIDSEM */

	// Quad_res_find_p();	//a 	find 'p' ST 'a' is quad residue (mod p)
	// Legendre_symbol();	//a,p
	// Miller_Rabin();	//n  -> giving false results too
	// Mersenne_number();	//n
	// Lucas_primality_test();	//n, a
	// Solovay_Strassen();	//n, a

	return 0;
}

void Solovay_Strassen() {
	int n; cin >> n;
	cout << " \n Solovay Strassen probablistic prime test on " << n << "\n\n";
	int a;
	cout << "\n\nChoose an a (1<a<n)\n";
	cin >> a;

	if (gcd(a, n) > 1) {
		cout << n << " is composite!! because, gcd(a,n) is " << gcd(a, n) << "\n";
		return;
	}

	cout << "Finding LHS = power(a,(n-1)/2,n) and RHS = Legendre_symbol(a,n)\n\n";
	int lhs = power(a, (n - 1) / 2, n);
	int rhs = (Legendre_symbol(a, n) + n) % n;

	cout << " LHS = " << lhs << "\n RHS = " << rhs << "\n";

	if (!iscong(lhs, rhs, n)) {
		cout << n << " is composite!! because, lhs is not cong to rhs " << gcd(a, n) << "\n";
		return;
	}
	cout << n << " is probably a prime with given a = " << a << "\n";
}

void Lucas_primality_test() {
	int n;	cin >> n;
	cout << " \n Lucas Primality Test on " << n << "\n\n";
	cout << "Prime factors of " << n - 1 << " are ";
	vi factors = prime_factors(n - 1);

	int a;
	cout << "\n\nChoose an a (1<a<n)\n";
	cin >> a;

	if (iscong(power(a, n - 1, n), 1, n)) {
		cout << " a^(n-1) is cong to 1modn\n";
	}
	else {
		cout << "Choose another a!\n";
		return;
	}

	for (int q : factors) {
		if (iscong(power(a, (n - 1) / q, n), 1, n)) {
			cout << "Value of a failed for factor: " << q << " and power " << (n - 1) / q << "\nchoose another a!!\n";
			return;
		}
	}
	cout << endl << n << " is prime for the \"a\" value as " << a << "\n";

}

void primes_sieve() {
	// cout << "PRIME SEIVE\n";
	memset(prime, true, sizeof prime);
	int n = sizeof(prime) - 1;
	for (int p = 2; p * p <= n ; p++) {
		if (prime[p] == true) {
			for (int i = p * p; i <= n; i += p)
				prime[i] = false;
		}
	}
	for (int p = 2; p <= n; p++)
		if (prime[p])
			primes.pb(p);

	// for0(i, primes.size()) {
	// 	cout << primes[i] << " -> ";
	// }
}

int extended_euclid_steps(int a, int b, int& x, int& y)  {
	if (b == 0) {
		x = 1, y = 0;
		return a;
	}
	else {
		int x1, y1;
		int d = extended_euclid_steps(b, a % b, x1, y1);
		x = y1, y = x1 - y1 * (a / b);
		cout << d << " = " << a << " * " << x << "  +  " << b << " * " << y << "\n";
		return d;
	}
}

void GCD() {
	int t;	cin >> t;
	while (t--) {
		int a, b;	cin >> a >> b;
		cout << gcd(a, b) << "\n";
	}
}

void GCD_steps(int a, int b) {
	if (a) {
		cout << b << " = " << b / a << "*(" << a << ") + " << b % a << "\n";
		GCD_steps(b % a, a);
	}
	else {
		cout << "GCD = " << b << "\n";
	}
}

void GCD_euler_steps() {
	cout << "GCD WITH STEPS\n";
	int a, b;
	cin >> a >> b;
	cout << "\n GCD of (" << a << ", " << b << ") is :" << "\n";
	if (a > b) {
		swap(a, b);
	}

	GCD_steps(a, b);
}


void extended_euclid() {
	cout << "EXTENDED EUCLID\n";
	int a, b;	cin >> a >> b;
	if (a > b) {
		swap(a, b);
	}
	GCD_steps(a, b);
	cout << "	now, \n";
	int x = 0, y = 0;
	int d = extended_euclid_steps(a, b, x, y);
	cout << " Ans: " << d << " = " << a << " * " << x << "  +  " << b << " * " << y << "\n";
}

void Fermats() {
	cout << "#### FERMATS ####\n\n";
	int a, n;	cin >> a >> n;
	if (gcd(a, n) != 1)
	{	cout << " GCD of a and n is not 1, check!\n";
		return;
	}
	int lhs = power(a, n - 1, n);
	cout << "a^(n-1) mod n = " << lhs << "\n";
	if (iscong(lhs, 1, n))
		cout << " => Theorem says " << n << " is probablistic prime!\n";
	else
		cout << " => n is not a prime!\n";
}

int phi(int n) {
	primes_sieve();
	int res = n;
	if (prime[n]) {return n - 1;}
	for (int i = 2; i * i <= n; i++) {
		if (n % i == 0) {
			res /= i, res *= (i - 1);
			while (n % i == 0) n /= i;
		}
	}
	if (n > 1)	res /= n, res *= (n - 1);
	return res;
}

void PHI() {
	cout << "PHI - euler totient function\n";
	//Euler's Toitent function
	primes_sieve();
	//	IMPORTANT: https://primefan.tripod.com/Phi500.html
	int t;	cin >> t;
	while (t--) {
		int n;	cin >> n;
		int N = n;
		int res = n;

		if (prime[n]) {
			cout << "prime's totient -> " << n - 1 << "\n";
			continue;
		}

		cout << " Prime factors \n";
		// find all prime numbers in prime factorization of n
		for (int i = 2; i * i <= n; i++) {
			if (n % i == 0) {
				res /= i;
				res *= (i - 1);
				while (n % i == 0)
					n /= i;
			}
		}
		cout << "\n PHI of " << N << " : ";
		if (n > 1)
			res /= n, res *= (n - 1);
		cout << res << "\n";
	}
}

vector<int> prime_factors(int n) {
	primes_sieve();
	// finding the prime factors of n
	int k = 0;
	vi factors;
	for (auto i : primes) {
		if (i > n)	break;
		if (n % i == 0 && n > 1) {
			k = 0;
			while (n % i == 0 && n > 1) {
				n /= i;
				k++;
			}
			factors.pb(i);
			cout << i << "^" << k << " * ";
		}
	}
	return factors;
}

void prime_factors() {
	primes_sieve();
	int n;	cin >> n;
	// finding the prime factors of n
	int k = 0;
	for (auto i : primes) {
		if (i > n)	break;
		if (n % i == 0 && n > 1) {
			k = 0;
			while (n % i == 0 && n > 1) {
				n /= i;
				k++;
			}
			cout << i << "^" << k << " * ";
		}
	}
}

void primitive_root() {
	int n;	cin >> n;
	unordered_set<int> res;
	for (int a = 1; a < n; a++) {
		res.clear();
		for (int i = 1; i < n; i++) {
			int modulo = power(a, i, n);
			if (res.find(modulo) != res.end())	//already present and repeating
			{
				break;
			}
			res.insert(modulo);

		}
		cout << a << " -> ";
		for (auto i :  res) {
			cout << i << " ";
		}
		cout << "\n";

		if (res.size() == n - 1) {
			cout << a << " ans.\n";
			break;
		}
	}
}

void block(int x)
{
	vector<int> v;

	// Converting the decimal number
	// into its binary equivalent.
	cout << "Blocks for " << x << " : ";
	while (x > 0)
	{
		v.push_back(x % 2);
		x = x / 2;
	}

	// Displaying the output when
	// the bit is '1' in binary
	// equivalent of number.
	for (int i = 0; i < v.size(); i++)
	{
		if (v[i] == 1)
		{
			cout << "2^" << i << " ";
			if (i != v.size() - 1)
				cout << "+ ";
		}
	}
	cout << "\n\n";
}

void POWER_steps() {
	cout << "POWER:  a^b mod m steps\n\n";
	int a, b, m;	cin >> a >> b >> m;
	// writing b in binary form
	block(b);
	int p = 1, k = 1;
	int A = a;
	while (b > 0) {
		cout << A << " ^ " << k << " mod " << m << " = " << a << endl;
		if (b % 2) {
			p = p * a; p %= m;
			cout << "\tRes: p = " << a << "\n\n";
		}
		b >>= 1;
		// cout << a << " ^ " << b << "(mod " << m << ") * " << p << "(mod " << m << ")\n";
		a *= a;
		a %= m;
		k *= 2;
	}
	cout << "\n Ans: " << p % m << "\n";
}

void ORDER_amodm() {
	int a, m;
	cin >> a >> m;
	int h = -1;	//order
	for (int i = 1; i <= m * m; i++) {
		if (power(a, i, m) == 1) {
			h = i;
			break;
		}
	}
	cout << "Order of amodm ( " << a << "mod" << m << " ) is " << h << "\n";
	cout << "(min value of h>1 ST a^h congruent to 1modm => h is the order)\n";
}

int inv(int a, int m) {
	int m0 = m, t, q;
	int x0 = 0, x1 = 1;
	if (m == 1)
		return 0;
	// Apply extended Euclid Algorithm
	int num = 99999;
	while (a > 1 && num >= 0) {
		// q is quotient
		q = a / m;
		t = m;
		// m is remainder now, process same as
		// euclid's algo
		m = a % m, a = t;
		t = x0;
		x0 = x1 - q * x0;
		x1 = t;
		num--;
	}
	if (num == 0) {
		cout << "Infinite while loop\n";
		return 0;
	}
	// Make x1 positive
	if (x1 < 0)
		x1 += m0;
	return x1;
}

bool iscong(int x, int b, int m) {
	return ((b - x) % m == 0);
}

void findMinX(int num[], int rem[], int k) {
	// Compute product of all numbers
	int prod = 1;
	for (int i = 0; i < k; i++)
		prod *= num[i];
	cout << "Product of num's = " << prod << "\n\n";
	cout << "ni\tbi\t\txi\t\tNi\t\tbixiNi\n";
	// Initialize result
	int result = 0;
	// Apply above formula
	for (int i = 0; i < k; i++) {
		int pp = prod / num[i];
		int xi = inv(pp, num[i]);
		int xbN = rem[i] * xi * pp;
		result += xbN;
		cout << num[i] << "\t" << rem[i] << "\t\t" << xi << "\t\t" << pp << "\t\t" << xbN << "\n";
	}
	cout << "Summation (xibiNi) = " << result << "\n";
	int x = result % prod;
	// while (x <= 45) {
	cout << "\nx is: " << x << " mod " << prod << "\n\n";
	// x += 5;
	// }
	for0(i, k) {
		if (iscong(x, rem[i], num[i]) == 0)
			cout << "TEST FAILED at " << i << " !!\n";
		cout << "Passed: " << x << " " << rem[i] << " " << num[i] << "\n";
	}
}

void Chinese_rem() {
	cout << "CHINESE REMAINDER THEOREM\n\n";
	int n;	cin >> n;	//number of equations
	int num[n], rem[n];
	for0(i, n)	cin >> rem[i] >> num[i];
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			if (gcd(num[i], num[j]) != 1) {
				cout << "Chinese remainder theorem not valid in this case!\n As each par of modulo should be relatively prime for it to hold.\n";
				cout << num[i] << " and " << num[j] << " are not co - primes\n";
				return;
			}
		}
	}
	findMinX(num, rem, n);

	for0(i, n) {
		cout << "Phi " << num[i] << " -> " << phi(num[i]) << "\n";
	}
	/*
	Nixi = 1 mod ni
	xi = Ni^(phi(ni)-1) mod ni
	*/
}

void Chinese_rem(int n, int rem[], int num[]) {
	cout << "CHINESE REMAINDER THEOREM after conversions\n\n";
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			if (gcd(num[i], num[j]) != 1) {
				cout << "Chinese remainder theorem not valid in this case!\n As each par of modulo should be relatively prime for it to hold.\n";
				cout << num[i] << " and " << num[j] << " are not co - primes\n";
				return;
			}
		}
	}
	findMinX(num, rem, n);

	for0(i, n) {
		cout << "Phi " << num[i] << " -> " << phi(num[i]) << "\n";
	}
	/*
	Nixi = 1 mod ni
	xi = Ni^(phi(ni)-1) mod ni
	*/
}

void Chinese_converted() {
	cout << "Equations conversions - CHINESE REMAINDER\n\n";
	int n;	cin >> n;
	int a[n], b[n], m[n];
	for0(i, n)	cin >> a[i] >> b[i] >> m[i];
	int new_b;
	cout << "b's and m's for chinese remainder\n";
	cout << n << "\n";
	for0(i, n) {
		new_b = (power(a[i], phi(m[i]) - 1, m[i]) * b[i]) % m[i];
		b[i] = new_b;
		cout << b[i] << " " << m[i] << "\n";
	}
	Chinese_rem(n, b, m);
}

void INVERSE_euler() {
	cout << " INVERSE OF amodm\n";
	int a, m;	cin >> a >> m;
	cout << " Inverse is : a ^ (phi(m) - 1) (mod m) = \n\t";
	cout << power(a, phi(m) - 1, m) << "\n";
}

int order(int a, int m) {
	// least value of h for which:a^h cong 1modm
	int h = -1;	//order
	for (int i = 1; i <= m * m; i++) {
		if (power(a, i, m) == 1) {
			h = i;
			break;
		}
	}
	return h;
}

void polypower() {
	int x, y, b, m;	cin >> y >> b >> m;
	cout << "x^" << y << " cong " << b << " mod " << m << "\n";
	b %= m;
	for (int i = 1; i <= m * m; i++) {
		if (power(i, y, m) == b) {
			cout << "x is: " << i << "\n";
		}
	}
	cout << "No value found? \n";
}

void polypower_pow() {
	int y, x, b, m;	cin >> y >> b >> m;
	cout << y << "^x" << " cong " << b << " mod " << m << "\n";
	b %= m;
	for (int i = 1; i <= m; i++) {
		if (power(y, i, m) == b) {
			cout << "x is: " << i << "\n";
		}
	}
	cout << "No value found? \n";
}

void PRIMITIVE_ROOT_steps() {
	cout << "### PRIMITIVE ROOT ###\n\n";

	int m;	cin >> m;
	int phi_m = phi(m);
	cout << "Phi(" << m << ") is " << phi_m << "\n";
	cout << "NOTE : number of primitive roots of m \n\t = 0 or phi(phi(m)) = " << phi(phi(m)) << "\n";
	int a;
	bool found = 0;
	cout << "Orders of amodm for a = 2 to m - 1 : \n";
	for (a = 2; a < m; a++) {
		int o = order(a, m);
		cout << a << "mod" << m << " => " << o << "\n";
		if (o == phi_m) {
			cout << " We found a primitive root, a = " << a << "\n";
			// break;
			found = 1;
		}
	}
	if (found == 0)
		cout << " No primitive root exists for this n!\n";
}

int smallest_p_root(int m) {
	int phi_m = phi(m);
	int a;
	for (a = 2; a < m; a++) {
		int o = order(a, m);
		if (o == phi_m) {
			return a;
		}
	}
	return -1;
}

void INDEX_of_primitive_root() {
	// m should be prime
	cout << "### INDEX_of_primitive_root ###\n\n";
	int m;	cin >> m;
	int phi_m = phi(m);
	int g = smallest_p_root(m);
	cout << "Smallest primitive root of " << m << " is " << g << "\n";
	cout << "Now, a congruent g ^ k mod m \n\t = a cong " << g << "^k mod" << m << "\n";
	vi index;
	for (int a = 1; a < m; a++) {
		for (int k = 0; k < m; k++) {
			if (a == power(g, k, m)) {
				index.pb(k);
				cout << " a = " << a << ", k = " << k << "\n";
				break;
			}
		}
	}

	cout << "Values of 'a <- k' such that k and phi(m) are co prime\n";
	for (int i = 0; i < index.size(); i++) {
		int k = index[i];
		if (gcd(k, phi_m) == 1) {
			cout << i + 1 << " < - " << k << "\n";
		}
	}
}

void Linear_cong() {
	int t;	cin >> t;
	while (t--) {
		int a, b, m, M, A, B;
		// ax = bmodm
		cin >> a >> b >> m;
		M = m;
		A = a, B = b;
		// has a soluion only when gcd(a,m)|b
		int g = gcd(a, m);
		cout << " STEP-1 : \ngcd(a,m) is " << g << "\n";

		if (b % g) {
			cout << "GCD of a, m does not divide v, NO SOLUTION!\n";
			return;
		}
		a /= g, b /= g, m /= g;
		cout << "ax cong b mod m \n\t = > " << a << "x cong " << b << " mod " << m << "\n";
		cout << "\n STEP-1: \nx = (a^(phi(m)-1))*b mod m\n";
		int phi_m = phi(m);
		cout << "\n phi = " << phi_m << "\n";
		int apow = power(a, (phi_m - 1), m);
		int x = (apow) * b % m;
		cout << "x = " << a << "^" << phi_m - 1 << " * " << b << " mod " << m << "\n";
		cout << "x = " << apow << " * " << b << " mod " << m << "\n";
		cout << "x = " << apow * b << " mod " << m << "\n";

		cout << "\nSolutions: \n";
		while (x <= M) {
			cout << x << ", ";
			if (iscong((A * x), B, M) == 0)
				cout << "Test FAILED!!!\n";
			x += m;
		}
		if (iscong((A * x), B, M)) {
			cout << "Test passed\n\n";
		}
		else cout << "Test FAILED!!!\n\n";
	}
}

void Jacobi_symbol(int a, int q) {
	if (a > q) {
		cout << a << " / " << q << " = ";
		a %= q;
		cout << a << " / " << q << "\n";
	}
	else if (a < 0) {
		cout << a << " / " << q << " = ";
		a = (a + q) % q;
		cout << a << " / " << q << "\n";
	}

	// factorize q
	vector<pair<int, int>> factsq;	// factors of q
	vector<pair<int, int>> facts;	// factors of a
	int ANS = 1;
	if (!prime[q]) {
		// factorize a
		cout << "Denominator factorized into \n\t";
		cout << q << " = ";

		int Q = q;
		int k = 0;
		for (auto i : primes) {
			if (i > Q)	break;
			if (Q % i == 0 && Q > 1) {
				k = 0;
				while (Q % i == 0 && Q > 1) {
					Q /= i;
					k++;
				}
				cout << i << "^" << k << " * ";
				factsq.pb(make_pair(i, k));
			}
		}
		cout << "\n";
	}
	else factsq.pb(make_pair(q, 1));

	for (auto y : factsq) {
		facts.clear();
		int p = y.first;
		if (!prime[a]) {
			// factorize a
			cout << "Numerator factorized into \n\t";
			cout << a << " = ";

			int k = 0;
			for (auto i : primes) {
				if (i > a)	break;
				if (a % i == 0 && a > 1) {
					k = 0;
					while (a % i == 0 && a > 1) {
						a /= i;
						k++;
					}
				}
				cout << i << "^" << k << " * ";
				facts.pb(make_pair(i, k));
			}
			cout << "\n";
		}
		else facts.pb(make_pair(a, 1));
		int ans;
		ans = 1;
		for (auto x : facts) {
			int ai = x.first;
			int po = x.second;
			// p = p % ai;
			// ((ai^po)/p)
			int legx;
			if (po % 2 == 0) {
				legx = 1;
			}
			else if (ai == 2) {
				legx = power(-1, (p - 1) / 8);
			}
			else {
				legx = power(ai, (p - 1) / 2, p);
				if (legx == p - 1) {
					legx = -1;
				}
			}
			cout << " ai, po, p & (ai ^ po) / p are \n " << ai << ",  " << po << ",  " << p << ",  " << legx << "\n";
			ans *= legx;
		}
		cout << "Legendre symbol value of " << a << " / " << p << " is: " << ans << "\n\n";
		ANS *= ans;
		ANS = power(ANS, y.second);
	}
	cout << "\nJacobi symbol value of " << a << " / " << q << " is: " << ANS << "\n\n";
}

int Legendre_symbol(int a, int n) {
	primes_sieve();
	vector<pair<int, int>> facts;
	if (!prime[n]) {
		cout << "##### IMPORTANT: NOT A PRIME!! #####\n\n";
		// Jacobi_symbol(a, p);
		// return -5;
	}
	else
		cout << "*** A PRIME!! ***\n\n";

	cout << "##### Legendre Symbol of " << a << "/" << n << " #####\n";
	if (!a) {
		cout << "Legendre symbol value: " << 0 << "\n\n";
		return 0;// (0/n) = 0
	}

	int ans = 1;
	if (a < 0)
	{
		a = -a; // (a/n) = (-a/n)*(-1/n)
		if (n % 4 == 3)
			ans = -ans; // (-1/n) = -1 if n = 3 (mod 4)
	}

	if (a == 1) {
		cout << "Legendre symbol value: " << ans << "\n\n";
		return ans;// (1/n) = 1
	}

	while (a)
	{
		if (a < 0)
		{
			a = -a;// (a/n) = (-a/n)*(-1/n)
			if (n % 4 == 3)
				ans = -ans;// (-1/n) = -1 if n = 3 (mod 4)
		}

		while (a % 2 == 0)
		{
			a = a / 2;
			if (n % 8 == 3 || n % 8 == 5)
				ans = -ans;

		}

		swap(a, n);

		if (a % 4 == 3 && n % 4 == 3)
			ans = -ans;
		a = a % n;

		if (a > n / 2)
			a = a - n;

	}

	if (n == 1) {
		cout << "Legendre symbol value: " << ans << "\n\n";

		return ans;
	}
	cout << "Legendre symbol value: " << 0 << "\n\n";
	return 0;

}


void Legendre_symbol() {
	int a;	cin >> a;
	int p;	cin >> p;
	primes_sieve();
	vector<pair<int, int>> facts;
	if (!prime[p]) {
		cout << "##### Jacobi Symbol #####\n";
		Jacobi_symbol(a, p);
		return;
	}
	cout << "##### Legendre Symbol #####\n";
	if (a > p) {
		cout << a << " / " << p << " = ";
		a %= p;
		cout << a << " / " << p << "\n";
	}
	else if (a < 0) {
		cout << a << " / " << p << " = ";
		a = (a + p) % p;
		cout << a << " / " << p << "\n";
	}
	if (!prime[a]) {
		// factorize a
		cout << "Numerator factorized into\n\t";
		cout << a << " = ";

		int k = 0;
		for (auto i : primes) {
			if (i > a)	break;
			if (a % i == 0 && a > 1) {
				k = 0;
				while (a % i == 0 && a > 1) {
					a /= i;
					k++;
				}
			}
			cout << i << "^" << k << " * ";
			facts.pb(make_pair(i, k));
		}
		cout << "\n";
	}
	else facts.pb(make_pair(a, 1));
	int ans = 1;
	for (auto x : facts) {
		int ai = x.first;
		int po = x.second;
		// ((ai^po)/p)
		int legx;
		if (po % 2 == 0) {
			legx = 1;
		}
		else if (ai == 2) {
			legx = power(-1, (p - 1) / 8);
		}
		else {
			legx = power(ai, (p - 1) / 2, p);
			if (legx == p - 1) {
				legx = -1;
			}
		}
		cout << " ai, po, p & (ai ^ po) / p are \n " << ai << ",  " << po << ",  " << p << ",  " << legx << "\n";
		ans *= legx;
	}
	cout << "Legendre symbol value: " << ans << "\n\n";

}

void ELGAMEL() {
	cout << "------------ELGAMEL cryptosystem------------\n\n";
	// a prime p, it's orimitive root A (alpha), secret key a, calculated key B
	// Public keys: p, A, B
	// Private keys: a (random int chosen)
	// Range allowed: 1 <= A, a <= p-1

	// BOB
	int A, a, p;	cin >> A >> a >> p;
	// B = A^a mod p
	int B = power(A, a, p);
	cout << "Public keys of Bob: p, A, B: " << p << ", " << A << ", " << B << "\n";
	cout << "Private key of Bob: a = " << a << "\n";

	// Encryption - ALICE
	// using Bob's public key
	int k;	cin >> k;	// random number < p
	int message;	cin >> message;
	int x = message;	//alias for message
	// Eb(x,k) = (y1, y2)
	//	 where y1 = A^k mod p and y2 = (x*power(B,k)) mod p
	int y1, y2;
	y1 = power(A, k, p);
	y2 = power(x * power(B, k, p), 1, p);	// called masking the plaintext!
	cout << "Encryption:\n y1 = " << y1 << "\n y2 = " << y2 << "\n";

	// Decryption - BOB
	// gets (y1,y2)
	// D(y1,y2) = y2 * ( (y1^a)^(-1) ) mod p

	int dec_message = power(y2 * inv(power(y1, a, p), p), 1, p);
	cout << "Decrypted message = " << dec_message << "\n";


}

void Knapsack_Merkle_Hellman() {
	cout << "------------Knapsack cryptosystem------------\n\n";
	int n;	cin >> n;
	vi W(n);	for0(i, n)	cin >> W[i];
	int W_sum = accumulate(W.begin(), W.end(), 0);
	cout << "sum of Wi is : " << W_sum << "\n";
	int q;	cin >> q;
	if (q < W_sum) {
		cout << "INVALID q\n";
		return;
	}
	int r;	cin >> r;
	if (gcd(q, r) != 1) {
		cout << "INVALID r\n";
		return;
	}
	cout << "Private keys accepted! (W, q, r)\n Generating public key B (bi = Wi * r(mod q))\n B is : ";
	vi B(n);
	for0(i, n) {
		B[i] = (W[i] * r) % 881;
		cout << B[i] << "  ";
	}
	cout << "\n\n";

	cout << "### ENCRYPTION in knapsack ###\n";
	cout << " input " << n << " digit message m\n";
	string m;	cin >> m;
	cout << "ciphertext : c = summation(mi * bi)\n";
	int c = 0;
	for0(i, n) {
		c += ((m[i] - '0') * B[i]);
	}
	cout << c << "  < - ciphertext\n\n";

	cout << "### DECRYPTION in knapsack ###\n";
	// we have c, and the private keys
	int r_dash, c_dash;
	r_dash = power(r, phi(q) - 1, q);
	c_dash = (c * r_dash) % q;
	cout << " r` : " << r_dash << "\n c` : " << c_dash << "\n";
	// APPLY Greedy to get those indices in W
	// represent c1 as sum of Wi using freedy
	vi indices;
	forr(i, n) {
		if (c_dash == 0)	break;
		if (W[i] <= c_dash) {
			indices.push_back(i + 1);
			c_dash -= W[i];
			// cout << "c` : " << c_dash << "\n";
		}
	}
	if (c_dash) {
		cout << " % % % % % % Something wrong in knapsack subset finsing step!\n";
		return;
	}
	int decrypted_m = 0;
	cout << "printing indices\n";
	for0(i, indices.size()) {
		cout << indices[i] << " ";
	}
	cout << "\n";
	for0(i, indices.size()) {
		decrypted_m += power(2, n - indices[i]);
	}
	cout << "Decrypted message is : " << decrypted_m << "\n\n";

	/*
	8
	2 7 11 21 42 89 180 354
	881
	588
	01100001
	*/
}

void RSA() {
	cout << "---------- -RSA cryptosystem------------\n\n";
	int p, q;	cin >> p >> q;
	int n = p * q;
	int phi_n = phi(n);
	cout << " n : " << n << "\n phi of n : " << phi_n << "\n";

	// given e
	cout << "choose e\n";
	int d, e;
	cin >> e;
	d = power(e, phi(phi_n) - 1, phi_n);
	if (gcd(e, phi_n) != 1)	{
		cout << "INVALID e\n";
		return;
	}
	if (gcd(d, phi_n) != 1)	{
		cout << "INVALID d\n";
		return;
	}
	if (iscong(d * e, 1, phi_n) == 0)	{
		cout << "INVALI e * d not cong 1 mod phi(n)\n";
		return;
	}
	cout << "e and d : " << e << " " << d << "\n";

	cout << " ed = 1 mod(phi_n) : " << " e" << " * d = 1 * mod( " << phi_n << " )\n";
	cout << " ed = 1 mod(phi_n) : " << e << " * " << d << " = 1 * mod( " << phi_n << " )\n";

	cout << "Public Keys : n, e : " << n << ", " << e << "\n";
	cout << "Private Keys : p, q, d : " << p << ", " << q << ", " << d << "\n";

	// Encryption:
	int ciphertext, plaintext;
	cin >> plaintext;
	cout << "Plaintext is : " << plaintext << "\n";
	cout << "E(m) = m ^ e modn\n" << "E(m) = " << plaintext << "^" << e << " mod" << n << "\n";
	ciphertext = power(plaintext, e, n);
	cout << "E(m)	= " << ciphertext << " \tis the ciphertext\n";

	// Decryption:
	cout << "D(m) = c ^ d modn\n" << "D(m) = " << ciphertext << "^" << d << " mod" << n << "\n";
	cout << "D(m) = " << power(ciphertext, d, n) << "\n";
}

void Miller_Rabin() {
	cout << " MILLER RABIN PRILIM TEST\n\n";
	int n;	cin >> n;

	primes_sieve();
	if (!prime[n]) {
		cout << "##### IMPORTANT: NOT A PRIME!! #####\n\n";
		// Jacobi_symbol(a, p);
		// return -5;
	}
	else
		cout << "*** A PRIME!! ***\n\n";


	// Step-1. Represent n-1 in the form of 2^s * t
	int s = 0, t = n - 1;
	while (t % 2 == 0) {
		t /= 2;
		s++;
	}
	cout << " n-1 = 2^s * t, s and t are: " << s << " " << t << "\n";
	cout << " => " << n << " = " << "2^" << s << " * " << t << " + 1\n";
	// Step-2. choose a random x from 2 to N-2
	int x = 2 + rand() % (n - 4);
	// x = 2;
	cout << "Chosen x is: " << x << "\n";

	int xpownminus1 = power(x, n - 1, n);
	if (!iscong(xpownminus1, 1, n)) {
		cout << " Ooops, x^(n-1) is not cong to 1modn\n";
		cout << " \n Composite !!! \n";
		return;
	}
	cout << " Yes, x^(n-1) is cong to 1modn, now the loop: \n\n";
	cout << " for any i=1 to s, if: \n Condition1 && Condition2 && Condition3 \n then, number is composite!\n";
	cout << " Condition 1 : x^(2^i) * t cong 1modn\n";
	cout << " Condition 2 : x^(2^(i-1)) * t not cong 1modn\n";
	cout << " Condition 2 : x^(2^(i-1)) * t not cong -1modn\n";


	int thisPow, B;
	for (int i = 1; i <= s; i++) {
		int tmp = 2 ^ i;
		thisPow = (power(x, tmp * t, n)) % n;
		B = (power(x, tmp / 2, n) * t) % n;
		bool c1, c2, c3;
		c1 = iscong(thisPow, 1, n);
		c2 = !iscong(B, 1, n) ;
		c3 = !iscong(B, - 1, n);
		cout << "i: c1 c2 c3 => " << i << " : " << c1 << " " << c2 << " " << c3 << " \n";
		if (c1 && c2 && c3) {
			cout << " Ooops, condition failed for i = " << i << "\n";
			cout << " \n Composite !!! \n";
			return;
		}
	}
	cout << "\n Prime !!! \n";
}

void Mersenne_number() {
	cout << "\n this is Mersenne_number function \n";
	cout << "\n Input: a prime n \n";
	int n;	cin >> n;
	cout << "n is : " << n << "\n";
	// Mn = 2^n - 1, n is a known prime
	int Mn = pow(2, n) - 1;
	cout << " Mn is " << Mn << "\n\n";
	// Mn is a prime if and only if mn divides S(n-2)
	int Si = 4;
	for (int i = 0; i < n - 2; i++) {
		cout << "S" << i << " is " << Si << "\n";
		Si = pow(Si, 2) - 2;
	}
	int Sn_minus_2 = Si;
	int rem = Sn_minus_2 % Mn;
	cout << "\n Sn_minus_2 = " << Sn_minus_2 << "\n";
	if (!rem) {
		cout << " Mn is a prime!!\n";
	}
	else cout << " Mn is not a prime!! Rem is : " << rem << "\n";
}
