#include<bits/stdc++.h>
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
bool prime[100000007];

// CRYPTOSYSTEMS
void RSA();

// NUMBER THEORY
void primes_sieve();
int phi(int n);
void PHI();
vector<int> prime_factors(int n);
void prime_factors();
bool iscong(int x, int b, int m);

signed main() {
	fast;
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
#endif

	//////////////////////////////////////////////////////////////////////////////////////
	/* NUMBER THEORY  */

	// PHI();					//t,n
	// primes_sieve();
	// prime_factors();			//n

	// cout<<prime_factors(70)<<"\n";
	// int d = (c - 1) / 2;
	// cout << gcd(383, 491389) << "\n";
	// cout << power(2 , 420, 491389) << "\n";

	// cout << (power(157, phi(2668) - 1, 2668));
	// cout << phi(96);// cout << iscong(17, 3, 7) << "\n";
	// cout << phi(55465219);
	// cout << endl << 23455 % 751 << " " << 12121 % 9721 << "  " << 10346 % 12121;

	RSA();			//p,q, e, message(plaintext)

	return 0;
}


void RSA() {
	cout << "---------- RSA cryptosystem------------\n\n";
	int p, q;	cin >> p >> q;
	int n = p * q;
	int phi_n = (p - 1) * (q - 1);
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
	string s, t;
	cin >> s;
	for (int i = 0; i < s.length(); i++) {
		// cout << s[i];
		int a = (s[i] - 'A' + 1) % 26;
		cout << a << " ";
		t += a;
	}
	cout << "t = " << t << "\n";
	// cin >> plaintext;
	plaintext = 31951492011;
	cout << "Plaintext is : " << plaintext << "\n";
	cout << "E(m) = m ^ e modn\n" << "E(m) = " << plaintext << "^" << e << " mod" << n << "\n";
	ciphertext = power(plaintext, e, n);
	cout << "E(m) = " << ciphertext << " \tis the ciphertext\n";

	// ciphertext = 123;
	// Decryption:
	cout << "D(m) = c ^ d modn\n" << "D(m) = " << ciphertext << "^" << d << " mod" << n << "\n";
	cout << "D(m) = " << power(ciphertext, d, n) << "\t is the plaintext\n";
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
		int res = n;

		if (prime[n]) {
			cout << "prime's totient -> " << n - 1 << "\n";
			continue;
		}

		// find all prime numbers in prime factorization of n
		for (int i = 2; i * i <= n; i++) {
			if (n % i == 0) {
				res /= i;
				res *= (i - 1);

				while (n % i == 0)
					n /= i;
			}
		}
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

bool iscong(int x, int b, int m) {
	return ((b - x) % m == 0);
}